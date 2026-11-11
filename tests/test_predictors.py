"""Predictor policy: which predictors run, which count, prefilter safety, suitability report."""
import random
import shutil
from collections import Counter

import pytest

from smithhunter import cli
from smithhunter.targets import predictors as pred
from smithhunter.targets import score
from smithhunter.targets.regions import Region, revcomp, write_regions
from smithhunter.targets.shuffle import make_decoys
from smithhunter.seqio import read_tsv

GUIDE = "TAGCTTATCAGACTGATGTTGA"
M8 = revcomp(GUIDE[1:8])
QUIET = lambda *_: None  # noqa: E731


def rnd(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def setup_regions(decoys):
    shared = rnd(120, 11) + M8 + "A" + rnd(120, 12)
    planted = rnd(100, 13)
    for _, _, seq in decoys[:3]:
        planted += revcomp(seq[1:8]) + "A" + rnd(40, 14)
    return [
        Region("tx1.3utr", "tx1", "geneA", "three_prime_UTR", shared, "chr1", "+", [(0, len(shared))]),
        Region("tx2.3utr", "tx2", "geneB", "three_prime_UTR", shared, "chr1", "+", [(0, len(shared))]),
        Region("tx3.3utr", "tx3", "geneC", "three_prime_UTR", planted, "chr2", "+", [(0, len(planted))]),
    ]


def faked(intarna="separates", rnahybrid="separates", calls=None):
    """A score_one in which each predictor either separates real from decoy sites (real is
    better than every decoy site) or is unsuited (real is worse than every decoy site)."""
    def one(guide, target, site, xi, cfg):
        real = guide == GUIDE
        jitter = (sum(map(ord, guide)) % 7) / 10
        out = {}
        if "intarna" in cfg.predictors:
            if calls is not None:
                calls["intarna"] += 1
            good = intarna == "separates"
            out.update(intarna_E=(-25.0 if real else -12.0 + jitter) if good
                       else (-10.0 if real else -30.0 + jitter))
        if "rnahybrid" in cfg.predictors:
            if calls is not None:
                calls["rnahybrid"] += 1
            good = rnahybrid == "separates"
            out.update(rnahybrid_mfe=(-30.0 if real else -15.0 + jitter) if good
                       else (-10.0 if real else -35.0 + jitter),
                       rnahybrid_p=0.01, rnahybrid_pos=10, rnahybrid_at_site=1)
        return out
    return one


def run(monkeypatch, one, cfg, log=QUIET):
    monkeypatch.setattr(score, "score_one", one)
    monkeypatch.setattr(score, "check_tools", lambda c: None)
    decoys = make_decoys([("c1", GUIDE)], 20, (2, 8), 5)
    return score.run_scoring([("c1", GUIDE)], decoys, setup_regions(decoys), (2, 8), None, cfg,
                             ".", log=log)


@pytest.fixture(autouse=True)
def _work_in_tmp(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)


# ---- policy: which predictors are computed ----------------------------------------------------

def test_mode_defaults():
    assert pred.resolve_predictors("ago") == ("intarna", "rnahybrid")
    assert pred.resolve_predictors("bacterial") == ("intarna",)
    assert pred.resolve_predictors("cross-kingdom") == ("intarna", "rnahybrid")


def test_target_set_overrides_global_overrides_mode():
    assert pred.resolve_predictors("bacterial", None, ["intarna", "rnahybrid"]) == ("intarna", "rnahybrid")
    assert pred.resolve_predictors("ago", ["rnahybrid"], ["intarna"]) == ("rnahybrid",)
    assert pred.resolve_predictors("ago", "intarna, rnahybrid") == ("intarna", "rnahybrid")
    assert pred.resolve_predictors("ago", [], ["intarna"]) == ()      # explicit: seed sites only


def test_unknown_predictor_or_mode_is_an_error():
    with pytest.raises(ValueError, match="unknown predictor"):
        pred.resolve_predictors("ago", ["intarna", "miranda"])
    with pytest.raises(ValueError, match="unknown mode"):
        pred.resolve_predictors("plant")


def test_counted_predictors():
    both = ("intarna", "rnahybrid")
    assert pred.counted_predictors(both) == both
    assert pred.counted_predictors(both, ["intarna"]) == ("intarna",)
    assert pred.counted_predictors(("intarna",), ["intarna", "rnahybrid"]) == ("intarna",)
    with pytest.raises(ValueError, match="none of the scored"):
        pred.counted_predictors(("intarna",), ["rnahybrid"])
    with pytest.raises(ValueError, match="unknown predictor"):
        pred.counted_predictors(both, ["miranda"])


def test_a_skipped_predictor_is_not_run_and_its_columns_are_empty(monkeypatch, tmp_path):
    calls = Counter()
    cfg = score.ScoreConfig(predictors=("rnahybrid",), calibrate_samples=0)
    _, sites, targets, curves = run(monkeypatch, faked(calls=calls), cfg)
    assert calls["intarna"] == 0 and calls["rnahybrid"] > 0
    assert all("intarna_E" not in r and "intarna_q" not in r for r in sites)
    assert {c["predictor"] for c in curves} == {"rnahybrid"}
    score.write_scores(sites, targets, curves, "s.tsv", "t.tsv", "c.tsv")
    for row in read_tsv("s.tsv"):
        assert row["intarna_E"] == "" and row["intarna_q"] == "" and row["intarna_q_cand"] == ""
        assert "intarna" not in row["called_by"]
    for row in read_tsv("t.tsv"):
        assert row["best_intarna_E"] == "" and row["best_intarna_q"] == ""


# ---- calls: union, support_from, called_by, best_q --------------------------------------------

def test_calls_of_is_the_union_of_counted_predictors():
    both = ("intarna", "rnahybrid")
    assert pred.calls_of({"intarna": 0.01, "rnahybrid": 0.9}, both, 0.05) == (1, "intarna", 0.01)
    assert pred.calls_of({"intarna": 0.2, "rnahybrid": 0.03}, both, 0.05) == (1, "rnahybrid", 0.03)
    assert pred.calls_of({"intarna": 0.01, "rnahybrid": 0.02}, both, 0.05) == (2, "intarna,rnahybrid", 0.01)
    assert pred.calls_of({"intarna": 0.2, "rnahybrid": 0.3}, both, 0.05) == (0, "", 0.2)
    assert pred.calls_of({"intarna": None, "rnahybrid": None}, both, 0.05) == (0, "", "")
    assert pred.calls_of({"intarna": 0.04}, both, 0.05) == (1, "intarna", 0.04)   # one q missing
    assert pred.calls_of({"intarna": 0.05}, ("intarna",), 0.05)[0] == 1            # q <= fdr


def test_union_default_a_disagreeing_tool_never_removes_a_call(monkeypatch):
    cfg = score.ScoreConfig(calibrate_samples=0, resume=False)   # the fakes change between loops
    for good, bad in (("intarna", "rnahybrid"), ("rnahybrid", "intarna")):
        one = faked(**{good: "separates", bad: "unsuited"})
        _, sites, targets, _ = run(monkeypatch, one, cfg)
        for r in sites:
            assert r["support"] == 1 and r["called_by"] == good
            assert r["best_q"] == r[f"{good}_q"] <= 0.05 < r[f"{bad}_q"]
        assert all(t["called_by"] == good and t["support"] == 1 for t in targets)


def test_both_agree_gives_support_two_and_both_called_by(monkeypatch):
    _, sites, targets, _ = run(monkeypatch, faked(), score.ScoreConfig(calibrate_samples=0))
    assert all(r["support"] == 2 and r["called_by"] == "intarna,rnahybrid" for r in sites)
    assert all(t["called_by"] == "intarna,rnahybrid" and t["best_q"] <= 0.05 for t in targets)


def test_support_from_ignores_a_predictor_but_still_reports_it(monkeypatch):
    # RNAhybrid separates, IntaRNA does not: counting IntaRNA only leaves no call, yet the
    # RNAhybrid score and q-value stay in the table
    one = faked(intarna="unsuited", rnahybrid="separates")
    cfg = score.ScoreConfig(calibrate_samples=0, support_from=("intarna",))
    messages = []
    _, sites, targets, _ = run(monkeypatch, one, cfg, log=messages.append)
    for r in sites:
        assert r["support"] == 0 and r["called_by"] == ""
        assert r["best_q"] == r["intarna_q"] > 0.05
        assert r["rnahybrid_q"] <= 0.05 and r["rnahybrid_mfe"] == -30.0
    assert all(t["support"] == 0 and t["called_by"] == "" and t["sites_q_pass"] == 0 for t in targets)
    assert any("IGNORED when forming calls" in m for m in messages)
    # the other way round, the same data and the call is made
    cfg = score.ScoreConfig(calibrate_samples=0, support_from=("rnahybrid",), resume=False)
    _, sites, _, _ = run(monkeypatch, one, cfg)
    assert all(r["support"] == 1 and r["called_by"] == "rnahybrid" for r in sites)


def test_support_from_does_not_change_the_scores_or_the_checkpoint(monkeypatch):
    base = score.ScoreConfig(calibrate_samples=0)
    other = score.ScoreConfig(calibrate_samples=0, support_from=("intarna",))
    args = ({"c1": GUIDE}, ["x"], (2, 8), ("8mer",))
    assert score.run_signature("final", *args, base) == score.run_signature("final", *args, other)


def test_targets_aggregate_called_by_as_union_and_best_q_as_minimum():
    cfg = score.ScoreConfig()
    base = {"candidate": "c", "transcript": "t", "gene": "g", "site_class": "8mer",
            "intarna_E": -10.0, "rnahybrid_mfe": -20.0}
    rows = [dict(base, support=1, called_by="intarna", best_q=0.04, intarna_q=0.04, rnahybrid_q=0.5),
            dict(base, support=1, called_by="rnahybrid", best_q=0.01, intarna_q=0.6, rnahybrid_q=0.01),
            dict(base, support=0, called_by="", best_q=0.7, intarna_q=0.7, rnahybrid_q=0.9)]
    (t,) = score.summarise_targets(rows, cfg)
    assert t["called_by"] == "intarna,rnahybrid" and t["best_q"] == 0.01 and t["support"] == 1
    assert t["sites_q_pass"] == 2


def test_written_tables_have_the_new_columns(monkeypatch):
    _, sites, targets, curves = run(monkeypatch, faked(), score.ScoreConfig(calibrate_samples=0))
    score.write_scores(sites, targets, curves, "s.tsv", "t.tsv", "c.tsv")
    for path in ("s.tsv", "t.tsv"):
        header = open(path).readline().split()
        assert {"support", "called_by", "best_q"} <= set(header)
    assert read_tsv("s.tsv")[0]["called_by"] == "intarna,rnahybrid"


# ---- the suitability diagnostic ---------------------------------------------------------------

def synthetic(n=200):
    """Real and decoy sites (one decoy per candidate): IntaRNA separates, RNAhybrid has the same
    scores for real and decoy sites, so it cannot tell them apart."""
    scores, qpool = {}, {}
    for i in range(n):
        v = -10.0 - i * 0.1
        real, decoy = ("c1", 0, i, i + 8, i + 7, "8mer"), ("c1__decoy1", 0, i, i + 8, i + 7, "8mer")
        scores[real] = {"intarna_E": v - 20.0, "rnahybrid_mfe": v}
        scores[decoy] = {"intarna_E": v, "rnahybrid_mfe": v}
        qpool[("intarna", real)] = 0.01
    return scores, {"c1__decoy1": "c1"}, Counter({"c1": 1}), qpool


def test_auc_basics():
    assert pred.weighted_auc([-5.0], [-1.0], [1.0]) == 1.0               # real is better (lower)
    assert pred.weighted_auc([-1.0], [-5.0], [1.0]) == 0.0
    assert pred.weighted_auc([-3.0], [-3.0], [1.0]) == 0.5               # ties count half
    assert pred.weighted_auc([1.0], [5.0], [1.0], higher_is_better=True) == 0.0
    assert pred.weighted_auc([-2.0], [-5.0, -1.0], [3.0, 1.0]) == 0.25   # weights: 1 of 4 is worse
    assert pred.weighted_auc([], [1.0], [1.0]) is None and pred.weighted_auc([1.0], [], []) is None
    assert pred.weighted_median([3.0, 1.0, 2.0]) == 2.0
    assert pred.weighted_median([1.0, 2.0], [1.0, 3.0]) == 2.0


def test_diagnostic_tells_a_separating_predictor_from_one_that_adds_nothing():
    scores, decoy_of, n_decoys, qpool = synthetic()
    rows = pred.diagnose(scores, decoy_of, n_decoys, ("intarna", "rnahybrid"), ("8mer",), qpool, 0.05)
    by = {r["predictor"]: r for r in rows}
    i, h = by["intarna"], by["rnahybrid"]
    assert i["auc"] == 1.0 and i["real_q_le_fdr"] == 200 and i["verdict"] == "separates"
    assert i["median_real"] < i["median_decoy"] and i["real_scored"] == 200 and i["decoy_scored"] == 200
    assert h["auc"] == 0.5 and h["real_q_le_fdr"] == 0 and h["verdict"] == "no_separation"
    assert "no separation from decoys: this predictor adds nothing for this data set" in h["note"]
    assert h["counted"] == "yes"


def test_diagnostic_marks_a_predictor_that_is_reported_only():
    scores, decoy_of, n_decoys, qpool = synthetic()
    rows = pred.diagnose(scores, decoy_of, n_decoys, ("intarna", "rnahybrid"), ("8mer",), qpool,
                         0.05, counted=("intarna",))
    h = [r for r in rows if r["predictor"] == "rnahybrid"][0]
    assert h["counted"] == "no" and "reported only" in h["note"]


def test_diagnostic_verdicts_for_the_other_cases():
    f = pred.advise
    assert f(0, 10, None, 0, 0.0, 0.05)[0] == "no_scores"
    assert f(10, 0, None, 0, 1.0, 0.05)[0] == "no_scores"
    assert f(10, 100, 0.9, 5, 1.0, 0.05)[0] == "too_few"
    assert f(500, 500, 0.30, 0, 1.0, 0.05)[0] == "worse_than_decoys"
    assert f(500, 500, 0.52, 0, 1.0, 0.05)[0] == "no_separation"        # inside 0.5 +- 0.02
    assert f(500, 500, 0.4801, 0, 1.0, 0.05)[0] == "no_separation"
    assert f(500, 500, 0.52, 3, 1.0, 0.05)[0] == "calls_without_separation"
    assert f(500, 500, 0.75, 0, 1.0, 0.05)[0] == "separates_no_calls"
    assert f(500, 500, 0.75, 9, 1.0, 0.05)[0] == "separates"
    note = f(500, 500, 0.75, 9, 0.15, 0.05)[1]
    assert "15% of the real sites" in note and "prefilter" in note


def test_diagnostic_uses_decoy_weights_and_skips_sites_without_decoys():
    scores = {("c1", 0, 0, 8, 7, "8mer"): {"intarna_E": -10.0},
              ("c2", 0, 0, 8, 7, "8mer"): {"intarna_E": -10.0},        # c2 has no decoys
              ("d1", 0, 0, 8, 7, "8mer"): {"intarna_E": -20.0},        # decoys better than real
              ("d2", 0, 0, 8, 7, "8mer"): {"intarna_E": -5.0}}
    rows = pred.diagnose(scores, {"d1": "c1", "d2": "c1"}, Counter({"c1": 2}), ("intarna",),
                         ("8mer",), {}, 0.05)
    assert rows[0]["real_sites"] == 1 and rows[0]["auc"] == 0.5        # one decoy better, one worse


def test_run_writes_the_report_and_logs_three_lines(monkeypatch, tmp_path):
    messages = []
    cfg = score.ScoreConfig(calibrate_samples=0, report_dir=str(tmp_path / "out"))
    run(monkeypatch, faked(), cfg, log=messages.append)
    rows = read_tsv(tmp_path / "out" / "predictor_report.tsv")
    assert {r["predictor"] for r in rows} == {"intarna", "rnahybrid"}
    assert set(rows[0]) >= {"site_class", "real_scored", "decoy_scored", "median_real",
                            "median_decoy", "auc", "real_q_le_fdr", "verdict", "note"}
    block = [m for m in messages if "predictor suitability" in m or m.startswith("  intarna:")
             or m.startswith("  rnahybrid:")]
    assert len(block) == 3
    pf = read_tsv(tmp_path / "out" / "prefilter_summary.tsv")
    assert "no prefilter" in pf[0]["note"]


def test_the_diagnostic_never_changes_the_calls(monkeypatch):
    # RNAhybrid is unsuited here (real sites worse than decoys) and says so in the report ...
    messages = []
    one = faked(intarna="separates", rnahybrid="unsuited")
    _, sites, _, _ = run(monkeypatch, one, score.ScoreConfig(calibrate_samples=0), log=messages.append)
    assert any(m.startswith("  rnahybrid:") for m in messages)      # it is reported ...
    # ... but the calls stay the union: IntaRNA still calls the site
    assert all(r["called_by"] == "intarna" for r in sites)


# ---- prefilter safety -------------------------------------------------------------------------

def test_prefilter_keep_classes_rule():
    scores = {("a", 0, 0, 8, 7, "8mer"): {"rnahybrid_mfe": -10.0},        # poor MFE
              ("b", 0, 0, 8, 7, "8mer"): {},                               # no RNAhybrid duplex
              ("c", 0, 0, 8, 7, "7mer-m8"): {"rnahybrid_mfe": -10.0},
              ("d", 0, 0, 8, 7, "7mer-m8"): {"rnahybrid_mfe": -30.0}}
    jobs = list(scores)
    cutoff = {"8mer": -20.0, "7mer-m8": -20.0}
    assert pred.prefilter_keep(jobs, scores, cutoff) == [jobs[3]]
    assert pred.prefilter_keep(jobs, scores, cutoff, ("8mer",)) == [jobs[0], jobs[1], jobs[3]]
    assert pred.prefilter_keep(jobs, scores, cutoff, ["8mer", "7mer-m8"]) == jobs


def test_prefilter_summary_counts_what_intarna_never_saw():
    scores = {("a", 0, 0, 8, 7, "8mer"): {"rnahybrid_mfe": -10.0},
              ("b", 0, 0, 8, 7, "8mer"): {},
              ("d1", 0, 0, 8, 7, "8mer"): {"rnahybrid_mfe": -30.0},
              ("e", 0, 0, 8, 7, "7mer-m8"): {"rnahybrid_mfe": -10.0},
              ("f", 0, 0, 8, 7, "7mer-m8"): {"rnahybrid_mfe": -30.0}}
    decoy_of = {"d1": "a"}
    cutoff = {"8mer": -20.0, "7mer-m8": -20.0}
    rows = pred.prefilter_summary(scores, decoy_of, cutoff, (), ("8mer", "7mer-m8"))
    real8 = [r for r in rows if r["site_class"] == "8mer" and r["group"] == "real"][0]
    assert (real8["sites"], real8["no_rnahybrid_score"], real8["scored_by_intarna"],
            real8["never_scored_by_intarna"], real8["share_never_scored"]) == (2, 1, 0, 2, 1.0)
    rows = pred.prefilter_summary(scores, decoy_of, cutoff, ("8mer",), ("8mer", "7mer-m8"))
    real8 = [r for r in rows if r["site_class"] == "8mer" and r["group"] == "real"][0]
    assert (real8["never_scored_by_intarna"], real8["kept_by_class_rule"]) == (0, 2)
    assert "always scored" in real8["note"]
    m7 = [r for r in rows if r["site_class"] == "7mer-m8"][0]
    assert (m7["sites"], m7["never_scored_by_intarna"], m7["kept_by_class_rule"]) == (2, 1, 0)


def two_stage_run(monkeypatch, keep, tmp_path):
    calls = Counter()
    messages = []
    cfg = score.ScoreConfig(calibrate_samples=0, intarna_prefilter=0.2,
                            prefilter_keep_classes=keep, report_dir=str(tmp_path / "rep"))
    _, sites, _, _ = run(monkeypatch, faked(calls=calls), cfg, log=messages.append)
    return calls, messages, sites, read_tsv(tmp_path / "rep" / "prefilter_summary.tsv")


def test_two_stage_warns_and_reports_what_intarna_never_scored(monkeypatch, tmp_path):
    calls, messages, sites, pf = two_stage_run(monkeypatch, (), tmp_path)
    assert any("NEVER scored by IntaRNA" in m for m in messages)
    assert any(m.startswith("prefilter cost:") and "never scored by IntaRNA" in m for m in messages)
    never = sum(int(r["never_scored_by_intarna"]) for r in pf)
    assert never > 0 and calls["intarna"] < calls["rnahybrid"]
    assert calls["rnahybrid"] - calls["intarna"] == never           # the file adds up with the calls
    assert sites[0]["support"] == 2                                  # the real site still passes


def test_keep_classes_make_intarna_score_every_site_of_the_class(monkeypatch, tmp_path):
    calls, messages, sites, pf = two_stage_run(monkeypatch, ("8mer",), tmp_path)
    assert calls["intarna"] == calls["rnahybrid"]                    # all sites here are 8mers
    assert all(int(r["never_scored_by_intarna"]) == 0 for r in pf)
    assert any("always scored: ['8mer']" in m for m in messages)


def test_keep_classes_change_the_checkpoint_signature_only_when_they_matter():
    args = ({"c1": GUIDE}, ["x"], (2, 8), ("8mer",))
    sig = lambda **kw: score.run_signature("final", *args, score.ScoreConfig(**kw))  # noqa: E731
    assert sig(intarna_prefilter=0.2) != sig(intarna_prefilter=0.2, prefilter_keep_classes=("8mer",))
    assert sig() == sig(prefilter_keep_classes=("8mer",))            # no prefilter: no effect
    assert score.run_signature("stage1", *args, score.ScoreConfig()) == \
        score.run_signature("stage1", *args, score.ScoreConfig(prefilter_keep_classes=("8mer",)))


def test_check_config_messages():
    msgs = []
    pred.check_config(score.ScoreConfig(intarna_prefilter=0.15, support_from=("intarna",)),
                      msgs.append)
    text = "\n".join(msgs)
    assert "NEVER scored by IntaRNA" in text and "better-suited tool" in text
    assert "rnahybrid is computed and reported but IGNORED" in text
    assert "still decides which sites IntaRNA sees" in text
    msgs = []
    pred.check_config(score.ScoreConfig(mode="bacterial", intarna_prefilter=0.15), msgs.append)
    assert any("not meant for bacterial" in m for m in msgs)
    msgs = []
    pred.check_config(score.ScoreConfig(predictors=("intarna",), intarna_prefilter=0.15), msgs.append)
    assert any("ignored" in m and "needs both predictors" in m for m in msgs)
    with pytest.raises(ValueError, match="prefilter_by"):
        pred.check_config(score.ScoreConfig(prefilter_by="intarna"), QUIET)
    with pytest.raises(ValueError, match="unknown mode"):
        pred.check_config(score.ScoreConfig(mode="plant"), QUIET)


# ---- command line -----------------------------------------------------------------------------

def cli_inputs(tmp_path):
    decoys = make_decoys([("c1", GUIDE)], 20, (2, 8), 5)
    (tmp_path / "cand.fa").write_text(f">c1\n{GUIDE}\n")
    (tmp_path / "decoys.fa").write_text("".join(f">{d}\n{s}\n" for d, _, s in decoys))
    write_regions(str(tmp_path / "regions.tsv"), setup_regions(decoys))
    return ["score-sites", "--candidates", str(tmp_path / "cand.fa"),
            "--decoys", str(tmp_path / "decoys.fa"), "--regions", str(tmp_path / "regions.tsv"),
            "--temperature", "25", "--calibrate-samples", "0", "--workdir", str(tmp_path / "w"),
            "--sites", str(tmp_path / "out" / "sites.tsv"),
            "--targets", str(tmp_path / "out" / "targets.tsv"),
            "--curves", str(tmp_path / "out" / "curves.tsv")]


def test_cli_bacterial_mode_runs_intarna_only_and_writes_the_reports(monkeypatch, tmp_path, capsys):
    calls = Counter()
    monkeypatch.setattr(score, "score_one", faked(calls=calls))
    monkeypatch.setattr(score, "check_tools", lambda c: None)
    (tmp_path / "out").mkdir()
    cli.main(cli_inputs(tmp_path) + ["--mode", "bacterial"])
    assert calls["intarna"] > 0 and calls["rnahybrid"] == 0
    rows = read_tsv(tmp_path / "out" / "sites.tsv")
    assert rows and all(r["rnahybrid_mfe"] == "" and r["rnahybrid_q"] == "" for r in rows)
    assert all(r["called_by"] == "intarna" for r in rows)
    report = read_tsv(tmp_path / "out" / "predictor_report.tsv")
    assert {r["predictor"] for r in report} == {"intarna"}
    assert (tmp_path / "out" / "prefilter_summary.tsv").exists()


def test_cli_explicit_predictors_and_support_from_override_the_mode(monkeypatch, tmp_path):
    calls = Counter()
    monkeypatch.setattr(score, "score_one", faked(intarna="unsuited", calls=calls))
    monkeypatch.setattr(score, "check_tools", lambda c: None)
    (tmp_path / "out").mkdir()
    cli.main(cli_inputs(tmp_path) + ["--mode", "bacterial", "--predictors", "intarna,rnahybrid",
                                     "--support-from", "rnahybrid"])
    assert calls["rnahybrid"] > 0
    rows = read_tsv(tmp_path / "out" / "sites.tsv")
    assert all(r["called_by"] == "rnahybrid" and r["intarna_q"] != "" for r in rows)
    report = {r["predictor"]: r["counted"] for r in read_tsv(tmp_path / "out" / "predictor_report.tsv")}
    assert report == {"intarna": "no", "rnahybrid": "yes"}


def test_cli_rejects_an_unknown_predictor(tmp_path):
    (tmp_path / "out").mkdir()
    with pytest.raises(SystemExit, match="unknown predictor"):
        cli.main(cli_inputs(tmp_path) + ["--predictors", "intarna,miranda"])


# ---- with the real tools (skipped where they are not installed) --------------------------------

needs_tools = pytest.mark.skipif(
    not (shutil.which("IntaRNA") and shutil.which("RNAhybrid")),
    reason="IntaRNA and RNAhybrid are not on PATH")


@needs_tools
def test_real_tools_intarna_only_skips_rnahybrid_and_both_report_the_diagnostic(tmp_path):
    target = rnd(150, 51) + revcomp(GUIDE[8:22]) + M8 + "A" + rnd(150, 52)
    regions = [Region("r", "t", "geneA", "three_prime_UTR", target, "c", "+", [(0, len(target))])]
    decoys = make_decoys([("c1", GUIDE)], 5, (2, 8), 3)
    only = score.ScoreConfig(predictors=("intarna",), calibrate_samples=0, resume=False,
                             report_dir=str(tmp_path / "only"))
    _, sites, targets, _ = score.run_scoring([("c1", GUIDE)], decoys, regions, (2, 8), None, only,
                                             str(tmp_path / "w1"), log=QUIET)
    assert sites and all(r["intarna_E"] < -15 for r in sites)
    assert all("rnahybrid_mfe" not in r and "rnahybrid_q" not in r for r in sites)
    assert {r["predictor"] for r in read_tsv(tmp_path / "only" / "predictor_report.tsv")} == {"intarna"}
    both = score.ScoreConfig(calibrate_samples=0, resume=False, support_from=("intarna",),
                             report_dir=str(tmp_path / "both"))
    _, sites, targets, _ = score.run_scoring([("c1", GUIDE)], decoys, regions, (2, 8), None, both,
                                             str(tmp_path / "w2"), log=QUIET)
    assert sites[0]["rnahybrid_mfe"] < -25 and sites[0]["rnahybrid_q"] != ""
    assert sites[0]["called_by"] in ("", "intarna")                 # never "rnahybrid": ignored
    report = {r["predictor"]: r for r in read_tsv(tmp_path / "both" / "predictor_report.tsv")}
    assert report["rnahybrid"]["counted"] == "no" and report["intarna"]["counted"] == "yes"


# ---- defaults keep the old results ------------------------------------------------------------

def test_default_output_is_what_it_was_before(monkeypatch):
    """No policy option set: both predictors computed, both counted, ``support`` is the old
    count of predictors with q <= fdr, the rows are sorted as before and every old column keeps
    its value; the new columns only add information."""
    cfg = score.ScoreConfig(calibrate_samples=0)
    assert cfg.predictors == ("intarna", "rnahybrid") and cfg.support_from is None
    assert cfg.prefilter_keep_classes == () and cfg.mode == "ago"
    one = faked(intarna="separates", rnahybrid="unsuited")
    _, sites, targets, _ = run(monkeypatch, one, cfg)
    for r in sites:
        legacy = sum(1 for p in cfg.predictors if r[f"{p}_q"] != "" and r[f"{p}_q"] <= cfg.fdr)
        assert r["support"] == legacy
    keys = [(r["candidate"], -r["support"]) for r in sites]
    assert keys == sorted(keys)
    for t in targets:
        assert t["support"] == max(r["support"] for r in sites
                                   if (r["candidate"], r["transcript"], r["gene"]) ==
                                   (t["candidate"], t["transcript"], t["gene"]))
        assert t["best_intarna_q"] != "" and t["best_rnahybrid_q"] != ""
    # an explicit policy that spells out the defaults gives the same numbers
    explicit = score.ScoreConfig(calibrate_samples=0, resume=False, predictors=("intarna", "rnahybrid"),
                                 support_from=("intarna", "rnahybrid"), prefilter_by="rnahybrid",
                                 prefilter_keep_classes=())
    _, sites2, targets2, _ = run(monkeypatch, one, explicit)
    assert sites2 == sites and targets2 == targets


def test_batch_options_reach_the_config(monkeypatch, tmp_path):
    seen = []
    real = score.run_scoring

    def spy(cands, decoys, regions, seed, classes, cfg, *rest, **kw):
        seen.append(cfg)
        return real(cands, decoys, regions, seed, classes, cfg, *rest, **kw)

    monkeypatch.setattr(score, "score_one", faked())
    monkeypatch.setattr(score, "check_tools", lambda c: None)
    monkeypatch.setattr(score, "run_scoring", spy)
    (tmp_path / "out").mkdir()
    cli.main(cli_inputs(tmp_path))
    assert (seen[-1].rnahybrid_batch, seen[-1].intarna_batch) == (2000, 50)
    cli.main(cli_inputs(tmp_path) + ["--rnahybrid-batch", "7", "--intarna-batch", "3"])
    assert (seen[-1].rnahybrid_batch, seen[-1].intarna_batch) == (7, 3)
