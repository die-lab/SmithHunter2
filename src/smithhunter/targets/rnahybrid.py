"""RNAhybrid wrapper for module B: score seed sites, calibrate p-values.

The window is the site with ``flank`` bases each side (the guide pairs upstream of the
seed on the target). The helix constraint ``-f`` is the seed range of the site class,
counted from the 5' end of the guide, so the best duplex must contain the seed helix.
p-values use the extreme value parameters from ``RNAcalibrate`` run on the actual target
regions (``-d xi,theta``) instead of the fly/worm/human presets; without calibration the
``3utr_human`` preset is used and the p-value is only indicative.

Two entry points give the same numbers: ``score_site`` (one process per site, the reference)
and ``score_sites`` (many windows of one small RNA in one process, the production path).
"""

from __future__ import annotations

import math
import os
import subprocess  # noqa: F401  (tests patch rnahybrid.subprocess.run)
import tempfile

from . import intarna
from .intarna import rna, seed_range
from .sites import Site

PRESET = "3utr_human"


def build_command(query: str, target: str, seed: tuple[int, int], xi_theta=None,
                  binary: str = "RNAhybrid") -> list[str]:
    cmd = [binary, "-c", "-b", "1", "-f", f"{seed[0]},{seed[1]}",
           "-m", str(max(len(target), 30) + 1), "-n", str(max(len(query), 30) + 1)]
    cmd += ["-d", f"{xi_theta[0]},{xi_theta[1]}"] if xi_theta else ["-s", PRESET]
    return cmd + [rna(target), rna(query)]


def parse_compact(text: str) -> dict | None:
    """Best duplex of RNAhybrid's compact output (``-c``):
    ``target:length:query:length:mfe:p-value:position:alignment...`` (colon separated;
    ``position`` is the 1-based start in the target)."""
    for line in text.splitlines():
        parts = line.strip().split(":")
        if len(parts) >= 7:
            try:
                mfe = float(parts[4])
                if mfe >= 0:  # no duplex satisfies the helix constraint (mfe 0.0, p 1.0)
                    return None
                out = {"rnahybrid_mfe": mfe, "rnahybrid_p": float(parts[5]),
                       "rnahybrid_pos": int(parts[6])}
                if len(parts) >= 9:  # target bases in the alignment: unpaired + paired rows
                    out["rnahybrid_span"] = sum(len(r.replace(" ", "")) for r in parts[7:9])
                return out
            except ValueError:
                continue
    return None


def window(target: str, site: Site, guide_len: int, flank: int = 15):
    """(sequence, lo): site window; the guide 3' end pairs upstream of the seed."""
    lo = max(0, site.t1 - (guide_len - 1) - flank)
    hi = min(len(target), site.t1 + 1 + flank)
    return target[lo:hi], lo


def _place(result: dict | None, lo: int, site: Site, require_site: bool) -> dict | None:
    """Region coordinates and the at-site test for a parsed duplex of the window that starts
    at ``lo`` (shared by the per-site and the batch path, so both give the same numbers)."""
    if result:
        result["rnahybrid_pos"] += lo  # region coordinates (1-based)
        # the duplex must cover the seed of this site, not another seed match in the window
        # (a seed-only site leaves the 3' half of the guide unpaired, so the start alone says
        # nothing). Seed on the target, 1-based inclusive: [t1 - b + 2, t1 - a + 2].
        a, b = seed_range(site.site_class)
        start = result["rnahybrid_pos"]
        end = start + result.get("rnahybrid_span", 0) - 1
        result["rnahybrid_at_site"] = int(start <= site.t1 - b + 3 and end >= site.t1 - a + 1)
        if require_site and not result["rnahybrid_at_site"]:
            return None
    return result


def score_site(guide: str, target: str, site: Site, xi_theta=None, flank: int = 15,
               binary: str = "RNAhybrid", require_site: bool = True,
               timeout: float | None = 300) -> dict | None:
    """With ``require_site`` a duplex that does not cover the seed of this site (the best
    helix of the window lies on another seed match) counts as no score, for real and decoy
    sites alike."""
    seq, lo = window(target, site, len(guide), flank)
    cmd = build_command(guide, seq, seed_range(site.site_class), xi_theta, binary)
    proc = intarna.run_tool(cmd, timeout)
    if proc.returncode != 0:
        raise RuntimeError(f"RNAhybrid failed ({proc.returncode}): {proc.stderr.strip()[:300]}")
    return _place(parse_compact(proc.stdout), lo, site, require_site)


# ---- batch path: many windows of one small RNA in ONE RNAhybrid process -------------------
# A call costs ~10 ms of process start-up for ~1 ms of work, so the windows of a chunk of
# sites go into one multi-FASTA (-t) and the per-target lines of the compact output are
# read back. All windows of a call share the small RNA, the helix constraint (the seed range
# of the site class) and xi/theta; -m only has to be at least the longest window.

def build_batch_command(query: str, seed: tuple[int, int], max_target: int, targets_path: str,
                        xi_theta=None, binary: str = "RNAhybrid") -> list[str]:
    cmd = [binary, "-c", "-b", "1", "-f", f"{seed[0]},{seed[1]}",
           "-m", str(max(max_target, 30) + 1), "-n", str(max(len(query), 30) + 1)]
    cmd += ["-d", f"{xi_theta[0]},{xi_theta[1]}"] if xi_theta else ["-s", PRESET]
    return cmd + ["-t", targets_path, rna(query)]


def _no_duplex(line: str) -> bool:
    """A well-formed result line (>= 7 fields, numeric energy >= 0): the helix constraint
    cannot be satisfied in that target."""
    parts = line.strip().split(":")
    if len(parts) < 7:
        return False
    try:
        return float(parts[4]) >= 0
    except ValueError:
        return False


def parse_compact_batch(text: str) -> dict[str, dict | None]:
    """Per target of a multi-FASTA call: the parsed best duplex, or None when the line says
    that no duplex satisfies the helix constraint (mfe >= 0). A target without a usable
    line is absent from the result: the caller must not read that as 'no duplex'. Lines are
    read with the rules of ``parse_compact``; the first usable line of a target wins."""
    out: dict[str, dict | None] = {}
    for line in text.splitlines():
        name = line.split(":", 1)[0]
        if name in out:
            continue
        parsed = parse_compact(line)
        if parsed is not None or _no_duplex(line):
            out[name] = parsed
    return out


def _batch_call(guide, windows, seed, xi_theta, binary, timeout, tmpdir):
    """windows: list of sequences. Returns {index: parsed duplex or None} for the targets
    the tool reported; raises RuntimeError when the call fails."""
    fd, path = tempfile.mkstemp(suffix=".fa", prefix="rnahybrid_batch_", dir=tmpdir)
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write("".join(f">w{i}\n{rna(w)}\n" for i, w in enumerate(windows)))
        cmd = build_batch_command(guide, seed, max(map(len, windows)), path, xi_theta, binary)
        proc = intarna.run_tool(cmd, timeout)
    finally:
        try:
            os.unlink(path)
        except OSError:
            pass
    if proc.returncode != 0:
        raise RuntimeError(f"RNAhybrid failed ({proc.returncode}): {proc.stderr.strip()[:300]}")
    parsed = parse_compact_batch(proc.stdout)
    return {int(k[1:]): v for k, v in parsed.items() if k[:1] == "w" and k[1:].isdigit()}


def score_sites(guide: str, items, xi_theta=None, flank: int = 15, binary: str = "RNAhybrid",
                require_site: bool = True, timeout: float | None = 300,
                tmpdir: str | None = None, min_split: int = 8) -> list:
    """Batch version of ``score_site``: ``items`` is a list of (target, Site) with the same
    seed range. Returns a list aligned with ``items``; each element is the result dict, None
    (no score) or, for a site that could not be scored, an Exception instance (the caller
    decides what an error means). Same numbers as ``score_site`` for every site.

    A chunk whose call fails (non-zero exit, timeout) is halved and retried, down to
    ``min_split`` sites, which are then scored one by one: a bad window costs its own site,
    not the chunk. A target missing from a successful output is re-scored alone rather than
    read as 'no duplex'. ``tmpdir`` receives the temporary FASTA (the work directory)."""
    if not items:
        return []
    seeds = {seed_range(site.site_class) for _, site in items}
    if len(seeds) != 1:
        raise ValueError("score_sites: all sites of a call must share the seed range")
    seed = seeds.pop()
    wins = [window(target, site, len(guide), flank) for target, site in items]
    out: list = [None] * len(items)

    def one(i):
        target, site = items[i]
        try:
            out[i] = score_site(guide, target, site, xi_theta, flank, binary, require_site,
                                timeout)
        except Exception as exc:  # noqa: BLE001  (reported per site)
            out[i] = exc

    def run(idx):
        try:
            got = _batch_call(guide, [wins[i][0] for i in idx], seed, xi_theta, binary,
                              timeout, tmpdir)
        except Exception:  # noqa: BLE001  (RuntimeError, OSError: halve, then single calls)
            if len(idx) <= max(1, min_split):
                for i in idx:
                    one(i)
            else:
                half = len(idx) // 2
                run(idx[:half])
                run(idx[half:])
            return
        for j, i in enumerate(idx):
            if j in got:
                out[i] = _place(got[j], wins[i][1], items[i][1], require_site)
            else:
                one(i)

    run(list(range(len(items))))
    return out


def parse_calibration(text: str) -> tuple[float, float] | None:
    """(xi, theta) from RNAcalibrate: ``query  <n>  xi  theta`` per query."""
    for line in text.splitlines():
        parts = line.split()
        if len(parts) >= 3:
            try:
                xi, theta = float(parts[-2]), float(parts[-1])
            except ValueError:
                continue
            # tiny samples can fit nonsense (inf, negative theta): refuse and use the preset
            if math.isfinite(xi) and math.isfinite(theta) and xi > 0 and theta > 0:
                return xi, theta
            return None
    return None


def calibrate(guide: str, targets_fasta: str, seed: tuple[int, int], mean_len: int,
              sd_len: int, samples: int = 2000, binary: str = "RNAcalibrate"):
    """Extreme value parameters of ``guide`` against dinucleotide-shuffled copies of the
    target regions in ``targets_fasta`` (RNA alphabet, FASTA)."""
    import os
    import tempfile
    # query as a FASTA file next to the targets (never in a system temp directory)
    with tempfile.NamedTemporaryFile("w", suffix=".fa", dir=os.path.dirname(targets_fasta),
                                     delete=False) as fh:
        fh.write(f">q\n{rna(guide)}\n")
        qpath = fh.name
    try:
        cmd = [binary, "-k", str(samples), "-l", f"{mean_len},{sd_len}",
               "-f", f"{seed[0]},{seed[1]}", "-m", str(mean_len + 6 * sd_len + 100),
               "-n", str(max(len(guide), 30) + 1), "-t", targets_fasta, "-s", "-q", qpath]
        proc = intarna.run_tool(cmd, 3600)
    finally:
        os.unlink(qpath)
    if proc.returncode != 0:
        raise RuntimeError(f"RNAcalibrate failed ({proc.returncode}): {proc.stderr.strip()[:300]}")
    return parse_calibration(proc.stdout)
