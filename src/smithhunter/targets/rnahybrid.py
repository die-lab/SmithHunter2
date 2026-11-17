"""RNAhybrid on one seed site, and RNAcalibrate: command builders and parsers.

RNAhybrid scores a window around the site: from (guide length + flank) nt upstream of
the seed match, where the 3' part of the guide can pair, to ``flank`` nt downstream.
The helix constraint ``-f`` is the paired stretch of the site class, so the seed used
by the scan is the seed RNAhybrid must form (v0 did not pass it). RNAhybrid cannot be
told where on the target the helix goes; ``rnahybrid_on_site`` says whether the duplex
it returned pairs the seed with the site of the scan (it does unless the window holds a
second, better seed match).

p-values need the parameters of an extreme value distribution. They are estimated per
guide by RNAcalibrate on windows sampled from the target regions themselves (``-d
xi,theta``), not taken from the fly / worm / human presets. The calibration runs without
the helix constraint: with it most random windows only form a forced, near-zero duplex
and the fit fails for about a quarter of the small RNAs. The p-value is therefore the
chance that a random window of the target regions pairs with the guide, anywhere and
with any structure, as well as the seed-anchored duplex does (conservative).
"""

from __future__ import annotations

import math
import random

from .intarna import seed_pairing
from .sites import Site

FLANK = 15


def window(target: str, site: Site, guide_length: int, seed: tuple[int, int] = (2, 8),
           flank: int = FLANK) -> tuple[str, int]:
    """Window of the target around the site and its offset in the region."""
    _, _, t_from, t_to = seed_pairing(site, seed)
    lo = max(0, t_from - guide_length - flank)
    hi = min(len(target), t_to + 2 + flank)
    return target[lo:hi], lo


def command(guide: str, win: str, site: Site, seed: tuple[int, int] = (2, 8),
            params: tuple[float, float] | None = None,
            binary: str = "RNAhybrid") -> list[str]:
    """RNAhybrid command for one window (compact output, best hit only).

    Without calibrated ``params`` a preset is passed only because RNAhybrid requires one;
    the p-value is then discarded by the caller.
    """
    g_from, g_to, _, _ = seed_pairing(site, seed)
    cmd = [binary, "-c", "-b", "1", "-f", f"{g_from},{g_to}",
           "-m", str(max(len(win), 1)), "-n", str(max(len(guide), 1))]
    cmd += ["-d", f"{params[0]},{params[1]}"] if params else ["-s", "3utr_human"]
    return cmd + [win, guide]


def duplex_pairs(lines: list[str], position: int, query_length: int) -> list[tuple[int, int]]:
    """(target position in the window, 1-based; guide position, 1-based) of each base
    pair, from the four alignment lines of the compact output: target unpaired, target
    paired, query paired, query unpaired. The target runs 5'->3', the query 3'->5'.
    ``position`` is the window position of the first target base shown."""
    width = max(len(x) for x in lines)
    t_un, t_pair, q_pair, q_un = (x.ljust(width) for x in lines)
    pairs = []
    t, q = position - 1, query_length + 1
    for c in range(width):
        has_t = t_un[c] != " " or t_pair[c] != " "
        has_q = q_un[c] != " " or q_pair[c] != " "
        t += has_t
        q -= has_q
        if t_pair[c] != " " and q_pair[c] != " ":
            pairs.append((t, q))
    return pairs


def parse_compact(text: str, site: Site, offset: int, seed: tuple[int, int] = (2, 8),
                  with_pvalue: bool = True) -> dict:
    """Scores from one line of ``RNAhybrid -c``; {} when there is no hit.

    Fields: target, target length, query, query length, mfe, p-value, position and the
    four alignment lines, separated by ':'.
    """
    line = next((ln for ln in text.splitlines() if ln.count(":") >= 10), "")
    if not line:
        return {}
    f = line.split(":")
    mfe, pvalue, position = float(f[4]), float(f[5]), int(f[6])
    pairs = duplex_pairs(f[7:11], position, int(f[3]))
    g_from, _, _, t_to = seed_pairing(site, seed)
    partner = {q: t for t, q in pairs}.get(g_from)
    t_span = [t for t, _ in pairs]
    return {
        "rnahybrid_mfe": mfe,
        "rnahybrid_p": pvalue if with_pvalue and not math.isnan(pvalue) else "",
        "rnahybrid_target": f"{min(t_span) + offset}-{max(t_span) + offset}" if pairs else "",
        "rnahybrid_pairs": len(pairs),
        "rnahybrid_on_site": int(partner is not None and partner + offset == t_to + 1),
    }


def calibration_windows(sequences: list[str], length: int, n: int, random_seed: int = 1
                        ) -> list[str]:
    """``n`` windows of ``length`` nt drawn from the target sequences, each sequence
    weighted by its length (shorter sequences are used whole)."""
    usable = [s for s in sequences if len(s) >= 20]
    if not usable:
        return []
    rng = random.Random(random_seed)
    picks = rng.choices(usable, weights=[len(s) for s in usable], k=n)
    out = []
    for s in picks:
        if len(s) <= length:
            out.append(s)
        else:
            p = rng.randrange(len(s) - length + 1)
            out.append(s[p:p + length])
    return out


def calibrate_command(target_fasta: str, query_fasta: str, max_target: int, max_query: int,
                      binary: str = "RNAcalibrate") -> list[str]:
    return [binary, "-m", str(max_target), "-n", str(max_query),
            "-t", target_fasta, "-q", query_fasta]


def parse_calibrate(text: str) -> dict[str, tuple[float, float]]:
    """{query: (xi, theta)} from RNAcalibrate (query, data points, xi, theta). Queries
    whose fit failed (nan, or too few data points) are left out."""
    out = {}
    for line in text.splitlines():
        f = line.split()
        if len(f) != 4:
            continue
        try:
            xi, theta = float(f[2]), float(f[3])
        except ValueError:
            continue
        if math.isfinite(xi) and math.isfinite(theta) and theta > 0:
            out[f[0]] = (xi, theta)
    return out
