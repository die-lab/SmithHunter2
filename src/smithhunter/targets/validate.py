"""Check scored targets against experimentally validated small RNA -> gene pairs.

With only a handful of validated pairs, a threshold cannot be tuned and the numbers are
not a precision/recall estimate. What this reports, per validated pair and per predictor:

- the rank (and percentile) of the validated gene among all genes of the target set for
  that small RNA, by its best site (lowest IntaRNA energy / RNAhybrid MFE) and by
  ``support``;
- the AUC (probability that a validated gene outranks a random gene of the same small
  RNA), pooled over all pairs.

Genes without any scored site (no strong seed site) rank last: missing the validated site
is a failure of the predictor and has to count. Percentile 0 is best.
"""

from __future__ import annotations

from collections import defaultdict

PREDICTORS = {"intarna": ("best_intarna_E", min), "rnahybrid": ("best_rnahybrid_mfe", min)}


def gene_table(target_rows):
    """(candidate, gene) -> best energies and support over its transcripts."""
    out: dict[tuple, dict] = {}
    for r in target_rows:
        key = (r["candidate"], r["gene"])
        g = out.setdefault(key, {"intarna": None, "rnahybrid": None, "support": 0, "q": None})
        for p, (col, fn) in PREDICTORS.items():
            v = r.get(col)
            if v not in (None, ""):
                v = float(v)
                g[p] = v if g[p] is None else fn(g[p], v)
        g["support"] = max(g["support"], int(r.get("support") or 0))
        for col in ("best_intarna_q", "best_rnahybrid_q"):
            q = r.get(col)
            if q not in (None, ""):
                g["q"] = float(q) if g["q"] is None else min(g["q"], float(q))
    return out


def ranks(values: dict[str, float | None], n_total: int):
    """gene -> (rank, percentile); lower value is better, missing is last (ties share the
    mean rank)."""
    scored = sorted((v, g) for g, v in values.items() if v is not None)
    out = {}
    i = 0
    while i < len(scored):
        j = i
        while j < len(scored) and scored[j][0] == scored[i][0]:
            j += 1
        mean_rank = (i + 1 + j) / 2
        for k in range(i, j):
            out[scored[k][1]] = mean_rank
        i = j
    # unscored genes share the last ranks: mean of len(scored)+1 .. n_total
    last = (len(scored) + 1 + n_total) / 2 if n_total > len(scored) else None
    return out, last


def validate(target_rows, truth, genes_per_candidate: dict[str, set] | None = None):
    """``truth``: rows with ``candidate`` and ``gene``. ``genes_per_candidate`` gives the
    full gene universe (from the regions); default is the genes that have a scored site.
    Returns (rows, summary)."""
    table = gene_table(target_rows)
    by_cand: dict[str, dict[str, dict]] = defaultdict(dict)
    for (cid, gene), g in table.items():
        by_cand[cid][gene] = g
    rows, auc_parts = [], {p: [0.0, 0] for p in list(PREDICTORS) + ["support"]}
    pcts: dict[str, list[float]] = {p: [] for p in auc_parts}
    for t in truth:
        cid, gene = t["candidate"], t["gene"]
        genes = by_cand.get(cid, {})
        universe = set(genes_per_candidate.get(cid, ())) if genes_per_candidate else set()
        universe |= set(genes) | {gene}
        n = len(universe)
        row = {"candidate": cid, "gene": gene, "n_genes": n, "scored": int(gene in genes),
               "note": t.get("note", "")}
        for p in list(PREDICTORS) + ["support"]:
            if p == "support":  # higher is better: rank on the negative
                vals = {g: -v["support"] if v["support"] else None for g, v in genes.items()}
            else:
                vals = {g: v[p] for g, v in genes.items()}
            rk, last = ranks(vals, n)
            rank = rk.get(gene, last if last is not None else n)
            pct = (rank - 1) / max(1, n - 1)
            row[f"{p}_rank"] = round(rank, 1)
            row[f"{p}_pct"] = round(pct, 4)
            # AUC share of this pair: fraction of other genes ranked below the validated one;
            # a universe of one gene says nothing, so such a pair is left out of the AUC
            if n >= 2:
                auc_parts[p][0] += 1 - pct
                auc_parts[p][1] += 1
                pcts[p].append(pct)
        q = genes.get(gene, {}).get("q")
        row["best_q"] = "" if q is None else q
        row["support"] = genes.get(gene, {}).get("support", 0)
        rows.append(row)
    summary = {p: (round(a / k, 4) if k else None) for p, (a, k) in auc_parts.items()}
    from statistics import median
    detail = {p: {"median_pct": round(median(v), 4), "top5": round(sum(x <= 0.05 for x in v) / len(v), 4),
                  "top10": round(sum(x <= 0.10 for x in v) / len(v), 4)}
              for p, v in pcts.items() if v}
    return rows, {"auc": summary, "detail": detail, "pairs": len(rows),
                  "found": sum(r["scored"] for r in rows)}


ROW_COLUMNS = ["candidate", "gene", "n_genes", "scored", "intarna_rank", "intarna_pct",
               "rnahybrid_rank", "rnahybrid_pct", "support_rank", "support_pct", "support",
               "best_q", "note"]
