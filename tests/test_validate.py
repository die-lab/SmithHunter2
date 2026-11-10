from smithhunter.targets.validate import gene_table, ranks, validate


def row(cand, gene, e, mfe, support=0, q=""):
    return {"candidate": cand, "gene": gene, "transcript": gene + ".1",
            "best_intarna_E": e, "best_rnahybrid_mfe": mfe, "support": support,
            "best_intarna_q": q, "best_rnahybrid_q": q}


def test_ranks_ties_and_missing():
    rk, last = ranks({"a": -10.0, "b": -10.0, "c": -5.0, "d": None}, 5)
    assert rk == {"a": 1.5, "b": 1.5, "c": 3.0}
    assert last == 4.5            # d and one gene without any site share ranks 4 and 5


def test_gene_table_takes_best_over_transcripts():
    rows = [row("c1", "g1", "-10", "-20", 1, "0.2"), row("c1", "g1", "-12", "-18", 2, "0.01")]
    g = gene_table(rows)[("c1", "g1")]
    assert g["intarna"] == -12.0 and g["rnahybrid"] == -20.0 and g["support"] == 2 and g["q"] == 0.01


def test_validate_ranks_validated_gene_and_penalises_missing():
    rows = [row("c1", "g1", "-30", "-35", 2), row("c1", "g2", "-20", "-25", 0),
            row("c1", "g3", "-15", "-22", 0), row("c1", "g4", "-12", "-21", 0)]
    truth = [{"candidate": "c1", "gene": "g1"}, {"candidate": "c1", "gene": "g9", "note": "no site"}]
    out, summary = validate(rows, truth, {"c1": {"g1", "g2", "g3", "g4", "g9"}})
    first, second = out
    assert first["intarna_rank"] == 1 and first["intarna_pct"] == 0.0 and first["support"] == 2
    assert second["scored"] == 0 and second["intarna_pct"] == 1.0   # unscored = last
    assert summary["pairs"] == 2 and summary["found"] == 1
    assert summary["auc"]["intarna"] == 0.5       # (1 + 0) / 2


def test_validate_universe_defaults_to_scored_genes():
    rows = [row("c1", "g1", "-30", "-35"), row("c1", "g2", "-20", "-25")]
    out, _ = validate(rows, [{"candidate": "c1", "gene": "g2"}])
    assert out[0]["n_genes"] == 2 and out[0]["intarna_rank"] == 2.0
