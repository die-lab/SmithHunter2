"""Empirical false discovery rates from decoy small RNAs.

Decoys go through the same scan as the real candidates. With ``n`` decoys per
candidate, the expected number of false sites is the decoy count divided by ``n``:

    FDR = (decoy sites / n) / real sites

``site_count_fdr`` does this per site class, for each candidate and pooled.
``qvalues`` does it along a score (used once energy predictors score each site).
"""

from __future__ import annotations

from bisect import bisect_right
from statistics import mean, pstdev


def site_count_fdr(real: dict[str, int], decoys: list[dict[str, int]], classes) -> list[dict]:
    """Rows of class, real sites, decoy mean and sd, enrichment and FDR."""
    rows = []
    for cls in list(classes) + ["any"]:
        pick = (lambda c: sum(c.values())) if cls == "any" else (lambda c, k=cls: c.get(k, 0))
        n_real = pick(real)
        counts = [pick(d) for d in decoys]
        m = mean(counts) if counts else 0.0
        rows.append({
            "site_class": cls, "real_sites": n_real,
            "decoy_mean": round(m, 3), "decoy_sd": round(pstdev(counts), 3) if counts else 0.0,
            "enrichment": round(n_real / m, 3) if m else "",
            "fdr": round(min(1.0, m / n_real), 4) if n_real else "",
        })
    return rows


def qvalues(real: list[float], decoy: list[float], n_decoys: int,
            higher_is_better: bool = False) -> list[float]:
    """q-value of each real score: the lowest FDR among thresholds it passes.

    By default lower scores are better (energies)."""
    if not real:
        return []
    sign = -1.0 if higher_is_better else 1.0
    r = sorted(sign * x for x in real)
    d = sorted(sign * x for x in decoy)
    n = max(1, n_decoys)
    # FDR at each distinct real score, from the best (lowest) upwards
    fdr_at = {}
    for x in sorted(set(r)):
        passed_real = bisect_right(r, x)
        passed_decoy = bisect_right(d, x)
        fdr_at[x] = min(1.0, passed_decoy / n / passed_real)
    q, best = {}, 1.0
    for x in sorted(fdr_at, reverse=True):
        best = min(best, fdr_at[x])
        q[x] = best
    return [round(q[sign * x], 6) for x in real]

