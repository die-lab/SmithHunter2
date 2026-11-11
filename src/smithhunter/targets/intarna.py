"""IntaRNA wrapper for module B: score seed sites.

Each site is scored in a local window of the target region (the site plus ``flank`` bases
each side, so the accessibility is computed in context: IntaRNA's own window is 150 nt).
The interaction is confined to the site (``--tRegion``) and the seed is pinned to the
guide seed range and to the seed match on the target (``--seedQRange``/``--seedTRange``),
so the energy is that of the interaction anchored on the site, not of the best duplex
anywhere in the window.

Coordinates: ``Site.start``/``end`` are 0-based half-open in the region; IntaRNA uses
1-based inclusive positions in the sequence it is given (here the window).

Two entry points give the same numbers: ``score_site`` (one process per site, the
reference) and ``score_pairs`` (several windows or several small RNAs in one process; see
there for what is shared and why it is exact).
"""

from __future__ import annotations

import os
import subprocess
import tempfile

from .sites import Site

# guide positions paired by the seed of each canonical class
SEED_OF_CLASS = {"8mer": (2, 8), "7mer-m8": (2, 8), "7mer-A1": (2, 7), "6mer": (2, 7),
                 "offset-6mer": (3, 8)}

CSV_COLUMNS = ["id1", "id2", "start1", "end1", "start2", "end2", "E", "ED1", "ED2",
               "Pu1", "Pu2", "E_hybrid", "hybridDP"]


def seed_range(site_class: str) -> tuple[int, int]:
    """Guide positions (1-based, inclusive) of the seed of a site class."""
    if site_class in SEED_OF_CLASS:
        return SEED_OF_CLASS[site_class]
    if site_class.startswith("seed"):  # seed<start>-<end>
        a, b = site_class[4:].split("-")
        return int(a), int(b)
    raise ValueError(f"unknown site class {site_class!r}")


def rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def window_coords(n: int, site: Site, guide_len: int, flank: int = 150, slack: int = 8):
    """Geometry of the window of a region of length ``n`` (no sequence needed):
    (lo, hi, t_region, seed_t). ``lo``/``hi`` are 0-based half-open in the region;
    ``t_region`` and ``seed_t`` are 1-based inclusive (from, to) in the window."""
    a, b = seed_range(site.site_class)
    t1 = site.t1
    reg_lo = max(0, t1 - (guide_len - 1) - slack)
    reg_hi = min(n, t1 + 1 + 3)
    lo = max(0, reg_lo - flank)
    hi = min(n, reg_hi + flank)
    seed_lo, seed_hi = t1 - b + 1, t1 - a + 2  # 0-based half-open
    return lo, hi, (reg_lo - lo + 1, reg_hi - lo), (seed_lo - lo + 1, seed_hi - lo)


def window(target: str, site: Site, guide_len: int, flank: int = 150, slack: int = 8):
    """Window of the target and the site coordinates inside it.

    Returns (sequence, lo, t_region, seed_t): ``lo`` is the window start in the region
    (0-based); ``t_region`` and ``seed_t`` are 1-based inclusive (from, to) in the window.
    The interaction may extend from the guide 3' end (``guide_len - 1`` bases upstream of
    t1 plus ``slack`` for bulges) to just past t1.
    """
    lo, hi, t_region, seed_t = window_coords(len(target), site, guide_len, flank, slack)
    return target[lo:hi], lo, t_region, seed_t


def build_command(query: str, target: str, site_class: str, t_region, seed_t,
                  temperature: float = 37.0, binary: str = "IntaRNA",
                  extra: list[str] | None = None) -> list[str]:
    a, b = seed_range(site_class)
    cmd = [binary, "-q", rna(query), "-t", rna(target), "--outMode", "C",
           "--outCsvCols", ",".join(CSV_COLUMNS),
           "--seedBP", str(b - a + 1), "--seedQRange", f"{a}-{b}",
           "--seedTRange", f"{seed_t[0]}-{seed_t[1]}",
           "--tRegion", f"{t_region[0]}-{t_region[1]}",
           "--temperature", str(temperature), "--threads", "1"]
    return cmd + list(extra or [])


def _result(row: dict) -> dict:
    def num(key):
        try:
            return float(row[key])
        except (KeyError, ValueError):
            return None

    return {
        "intarna_E": num("E"), "intarna_ED_target": num("ED1"), "intarna_ED_query": num("ED2"),
        "intarna_Pu_target": num("Pu1"), "intarna_Pu_query": num("Pu2"),
        "intarna_E_hybrid": num("E_hybrid"),
        "intarna_t_start": int(row["start1"]), "intarna_t_end": int(row["end1"]),
        "intarna_q_start": int(row["start2"]), "intarna_q_end": int(row["end2"]),
        "intarna_duplex": row.get("hybridDP", ""),
    }


def parse_csv(text: str) -> dict | None:
    """The best interaction of IntaRNA's CSV output (``--outMode C``), or None when no
    interaction satisfies the constraints."""
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    if len(lines) < 2:
        return None
    header = lines[0].split(";")
    return _result(dict(zip(header, lines[1].split(";"))))


def parse_csv_batch(text: str) -> dict[tuple[str, str], dict]:
    """All rows of a multi-target / multi-query call: {(target id, query id): result}. A
    pair with no row has no interaction under the constraints (as in ``parse_csv``). The
    first row of a pair wins (``--outNumber`` is 1, so there is one)."""
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    if len(lines) < 2:
        return {}
    header = lines[0].split(";")
    out: dict[tuple[str, str], dict] = {}
    for ln in lines[1:]:
        row = dict(zip(header, ln.split(";")))
        key = (row.get("id1", ""), row.get("id2", ""))
        if key not in out:
            out[key] = _result(row)
    return out


def run_tool(cmd: list[str], timeout: float | None = None):
    """subprocess.run with a timeout; a hang becomes an ordinary error of that site."""
    try:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"{cmd[0]} timed out after {timeout} s") from None


def score_site(guide: str, target: str, site: Site, flank: int = 150,
               temperature: float = 37.0, binary: str = "IntaRNA",
               extra: list[str] | None = None, timeout: float | None = 300) -> dict | None:
    seq, lo, t_region, seed_t = window(target, site, len(guide), flank)
    cmd = build_command(guide, seq, site.site_class, t_region, seed_t, temperature, binary,
                        extra)
    proc = run_tool(cmd, timeout)
    if proc.returncode != 0:
        raise RuntimeError(f"IntaRNA failed ({proc.returncode}): {proc.stderr.strip()[:300]}")
    result = parse_csv(proc.stdout)
    if result:  # report the interaction in region coordinates (1-based)
        result["intarna_t_start"] += lo
        result["intarna_t_end"] += lo
    return result


# ---- batch path ----------------------------------------------------------------------------
# IntaRNA spends ~90% of its time on the accessibility of the target window (a partition
# function over the 150 nt windows of the sequence it is given) and ~10% on start-up, the
# query and the interaction. Handing several targets to one process saves the start-up (and
# the query accessibility); handing several queries for ONE target saves the target
# accessibility as well, because it is computed once per target and reused for every query
# (measured: 8 queries on one window cost 0.1 s instead of 0.85 s, same numbers). The second
# case occurs when small RNAs share a window: same region, same t1 and the same layout
# (the same guide length and seed range), e.g. a candidate and its keep_seed decoys.
# ``--tRegion`` and ``--seedTRange`` apply to every target of a call, so the windows of one
# call must share the layout (``t_region``, ``seed_t``) and the seed range; windows of
# other layouts (sites near a region end) go into other calls and are never padded.

def _pair_single(query, target, site_class, t_region, seed_t, temperature, binary, extra,
                 timeout):
    """One pair through the reference command line; an exception instance on failure."""
    try:
        cmd = build_command(query, target, site_class, t_region, seed_t, temperature, binary,
                            extra)
        proc = run_tool(cmd, timeout)
        if proc.returncode != 0:
            raise RuntimeError(f"IntaRNA failed ({proc.returncode}): {proc.stderr.strip()[:300]}")
        return parse_csv(proc.stdout)
    except Exception as exc:  # noqa: BLE001  (reported per pair)
        return exc


def score_pairs(queries: list[str], windows: list[str], site_class: str, t_region, seed_t,
                temperature: float = 37.0, binary: str = "IntaRNA",
                extra: list[str] | None = None, timeout: float | None = 300,
                tmpdir: str | None = None) -> dict:
    """Every query against every window in ONE IntaRNA process (the caller makes sure that
    is the set of pairs it wants: one window with several queries, or one query with
    several windows). The windows share ``t_region`` / ``seed_t`` (window coordinates, as
    from ``window``) and the seed range of ``site_class``. Returns {(window index, query
    index): result in WINDOW coordinates, or None (no interaction), or an Exception}.

    If the call fails (non-zero exit, timeout), the larger side is halved and retried down to
    single pairs, which use the reference command line: a bad window costs its own pair only."""
    out: dict = {}
    nq, nw = len(queries), len(windows)
    if nq * nw == 0:
        return out
    if nq * nw == 1:
        out[(0, 0)] = _pair_single(queries[0], windows[0], site_class, t_region, seed_t,
                                   temperature, binary, extra, timeout)
        return out
    qpath = tpath = None
    try:
        fd, qpath = tempfile.mkstemp(suffix=".fa", prefix="intarna_q_", dir=tmpdir)
        with os.fdopen(fd, "w") as fh:
            fh.write("".join(f">q{i}\n{rna(q)}\n" for i, q in enumerate(queries)))
        fd, tpath = tempfile.mkstemp(suffix=".fa", prefix="intarna_t_", dir=tmpdir)
        with os.fdopen(fd, "w") as fh:
            fh.write("".join(f">t{i}\n{rna(w)}\n" for i, w in enumerate(windows)))
        cmd = build_command("ACGU", "ACGU", site_class, t_region, seed_t, temperature, binary,
                            extra)
        cmd[cmd.index("-q") + 1] = qpath
        cmd[cmd.index("-t") + 1] = tpath
        proc = run_tool(cmd, timeout)
        if proc.returncode != 0:
            raise RuntimeError(f"IntaRNA failed ({proc.returncode}): {proc.stderr.strip()[:300]}")
        rows = parse_csv_batch(proc.stdout)
    except Exception:  # noqa: BLE001  (halve the larger side; single pairs at the bottom)
        rows = None
    finally:
        for p in (qpath, tpath):
            if p:
                try:
                    os.unlink(p)
                except OSError:
                    pass
    if rows is not None:
        for wi in range(nw):
            for qi in range(nq):
                out[(wi, qi)] = rows.get((f"t{wi}", f"q{qi}"))
        return out
    args = (site_class, t_region, seed_t, temperature, binary, extra, timeout, tmpdir)
    if nw >= nq:
        half = nw // 2
        for off, part in ((0, windows[:half]), (half, windows[half:])):
            for (wi, qi), r in score_pairs(queries, part, *args).items():
                out[(wi + off, qi)] = r
    else:
        half = nq // 2
        for off, part in ((0, queries[:half]), (half, queries[half:])):
            for (wi, qi), r in score_pairs(part, windows, *args).items():
                out[(wi, qi + off)] = r
    return out
