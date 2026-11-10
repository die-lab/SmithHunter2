"""Empirical false discovery rates from decoy small RNAs.

Decoys go through the same scan as the real candidates. With ``n`` decoys per
candidate, the expected number of false sites is the decoy count divided by ``n``:

    FDR = (decoy sites / n) / real sites

``site_count_fdr`` does this per site class, for each candidate and pooled.
``qvalues`` does it along a score (used once energy predictors score each site).
A q-value of 0 means "better than every decoy site"; with ``n`` decoys per candidate
its resolution is about 1/n of a site, so it is not a probability of exactly zero.
"""

from __future__ import annotations

from bisect import bisect_right
from itertools import accumulate
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
            higher_is_better: bool = False,
            decoy_weights: list[float] | None = None,
            plus_one: float = 0.0) -> list[float]:
    """q-value of each real score: the lowest FDR among thresholds it passes.

    By default lower scores are better (energies). ``decoy_weights`` gives each decoy
    site its own weight (1 / number of decoys of its candidate, when candidates have
    different numbers of decoys); without it every decoy site weighs 1 / ``n_decoys``.
    ``plus_one`` is a pseudo decoy site (its weight, normally 1 / decoys of the candidate)
    added to the numerator at every threshold: FDR = (decoy weight + plus_one) / real, the
    (D + 1) / T estimate that keeps a site from reaching q = 0 just because no decoy of the
    few that were drawn scored better (Barber & Candes 2015; Phipson & Smyth 2010).
    Without decoy sites the FDR is 0 only if decoys were scored and none of their sites
    reached the threshold (``n_decoys`` >= 1 or ``decoy_weights`` given); with no decoys
    at all it is not estimable and raises ValueError instead of returning zeros."""
    if not real:
        return []
    if not decoy and n_decoys < 1 and decoy_weights is None:
        raise ValueError("q-values need decoys to estimate the false discovery rate")
    sign = -1.0 if higher_is_better else 1.0
    r = sorted(sign * x for x in real)
    pairs = sorted(zip((sign * x for x in decoy),
                       decoy_weights if decoy_weights is not None
                       else [1.0 / max(1, n_decoys)] * len(decoy)))
    d = [x for x, _ in pairs]
    cum = list(accumulate(w for _, w in pairs))
    fdr_at = {}
    for x in sorted(set(r)):
        passed_real = bisect_right(r, x)
        k = bisect_right(d, x)
        false_sites = cum[k - 1] if k else 0.0
        fdr_at[x] = min(1.0, (false_sites + plus_one) / passed_real)
    q, best = {}, 1.0
    for x in sorted(fdr_at, reverse=True):
        best = min(best, fdr_at[x])
        q[x] = best
    return [round(q[sign * x], 6) for x in real]
