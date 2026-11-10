"""RNAhybrid wrapper for module B: score one seed site, calibrate p-values.

The window is the site with ``flank`` bases each side (the guide pairs upstream of the
seed on the target). The helix constraint ``-f`` is the seed range of the site class,
counted from the 5' end of the guide, so the best duplex must contain the seed helix.
p-values use the extreme value parameters from ``RNAcalibrate`` run on the actual target
regions (``-d xi,theta``) instead of the fly/worm/human presets; without calibration the
``3utr_human`` preset is used and the p-value is only indicative.
"""

from __future__ import annotations

import math
import subprocess  # noqa: F401  (tests patch rnahybrid.subprocess.run)

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
    result = parse_compact(proc.stdout)
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
