"""IntaRNA wrapper for module B: score one seed site.

Each site is scored in a local window of the target region (the site plus ``flank`` bases
each side, so the accessibility is computed in context: IntaRNA's own window is 150 nt).
The interaction is confined to the site (``--tRegion``) and the seed is pinned to the
guide seed range and to the seed match on the target (``--seedQRange``/``--seedTRange``),
so the energy is that of the interaction anchored on the site, not of the best duplex
anywhere in the window.

Coordinates: ``Site.start``/``end`` are 0-based half-open in the region; IntaRNA uses
1-based inclusive positions in the sequence it is given (here the window).
"""

from __future__ import annotations

import subprocess

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


def window(target: str, site: Site, guide_len: int, flank: int = 150, slack: int = 8):
    """Window of the target and the site coordinates inside it.

    Returns (sequence, lo, t_region, seed_t): ``lo`` is the window start in the region
    (0-based); ``t_region`` and ``seed_t`` are 1-based inclusive (from, to) in the window.
    The interaction may extend from the guide 3' end (``guide_len - 1`` bases upstream of
    t1 plus ``slack`` for bulges) to just past t1.
    """
    a, b = seed_range(site.site_class)
    t1 = site.t1
    reg_lo = max(0, t1 - (guide_len - 1) - slack)
    reg_hi = min(len(target), t1 + 1 + 3)
    lo = max(0, reg_lo - flank)
    hi = min(len(target), reg_hi + flank)
    seed_lo, seed_hi = t1 - b + 1, t1 - a + 2  # 0-based half-open
    return (target[lo:hi], lo, (reg_lo - lo + 1, reg_hi - lo), (seed_lo - lo + 1, seed_hi - lo))


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


def parse_csv(text: str) -> dict | None:
    """The best interaction of IntaRNA's CSV output (``--outMode C``), or None when no
    interaction satisfies the constraints."""
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    if len(lines) < 2:
        return None
    header = lines[0].split(";")
    row = dict(zip(header, lines[1].split(";")))

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
