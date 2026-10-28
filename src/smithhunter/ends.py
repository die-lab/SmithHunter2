"""End-precision score of a locus, after sharp_smith.R in SmithHunter v0.

Positions are ranked by read count; the smallest set of top positions reaching
``n_thre`` of all reads is kept. With k kept positions spread over a span that
contains g positions not kept, the score is ``1 / (k + penalty * g)``:
1.0 for a single dominant end, 0.5 for two adjacent ends, lower for scattered ends.

Differences from v0, on purpose:
- a single dominant position scores 1.0 (v0 gave 1 / (1 - penalty), above 1);
- the gap counts positions inside the span that were not kept (v0 used span - 2
  after dropping positions with tied coverage values, which depended on ties).

The score takes few values (1, 0.5, 0.476, 0.455, 0.333, ...): ``score > 0.5`` holds
exactly when the top position alone reaches ``n_thre``. The 5' filter of a locus can
therefore use the share of reads at the top position directly (``dominant``), the share
within a few bases of it (``jitter``, tolerating 5' ends shifted by 1-2 nt), or the score
(``score``, as before); see ``FIVE_FILTERS``.
"""

from __future__ import annotations

FIVE_FILTERS = ("dominant", "jitter", "score")


def end_score(profile: dict[int, float], n_thre: float = 0.5, penalty: float = 0.1):
    """Return (score, fraction of reads at the dominant position)."""
    total = sum(profile.values())
    if total <= 0:
        return 0.0, 0.0
    ranked = sorted(profile.items(), key=lambda kv: (-kv[1], kv[0]))
    chosen = []
    accumulated = 0.0
    for position, count in ranked:
        chosen.append(position)
        accumulated += count
        if accumulated >= n_thre * total - 1e-9:
            break
    gaps = (max(chosen) - min(chosen) + 1) - len(chosen)
    return 1.0 / (len(chosen) + penalty * gaps), ranked[0][1] / total


def jitter_fraction(profile: dict[int, float], jitter: int = 1) -> float:
    """Share of reads within ``jitter`` positions of the top position (ties: leftmost)."""
    total = sum(profile.values())
    if total <= 0:
        return 0.0
    top = min(profile, key=lambda p: (-profile[p], p))
    return sum(c for p, c in profile.items() if abs(p - top) <= jitter) / total
