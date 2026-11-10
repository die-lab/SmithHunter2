"""Regression tests for the issues found in the review of module B step 2."""
import random
import subprocess

import pytest

from smithhunter.targets import intarna, rnahybrid, score
from smithhunter.targets.fdr import qvalues
from smithhunter.targets.regions import Region, revcomp
from smithhunter.targets.shuffle import make_decoys, seed_of
from smithhunter.targets.sites import Site, scan
from smithhunter.targets.validate import validate

G = "TGAGGTAGTAGGTTGTATAGTT"  # let-7a, 22 nt


def rnd(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


# ---- t1 of a non-canonical seed ---------------------------------------------------------

def test_t1_of_non_canonical_seed_points_at_the_guide_1_position():
    a, b = 4, 10
    target = rnd(30, 1) + revcomp(G[a - 1:b]) + rnd(30, 2)
    site = [s for s in scan(G, target, (a, b)) if s.start == 30][0]
    assert site.site_class == "seed4-10"
    assert site.t1 == site.start + b - 1            # guide position k pairs with t1 - (k - 1)
    seq, lo, t_region, seed_t = intarna.window(target, site, len(G), flank=40)
    assert seq[seed_t[0] - 1:seed_t[1]] == revcomp(G[a - 1:b])   # the pinned seed is the match


@pytest.mark.skipif(not __import__("shutil").which("IntaRNA"), reason="IntaRNA not on PATH")
def test_non_canonical_seed_is_scored_by_intarna():
    a, b = 4, 10
    target = rnd(60, 3) + revcomp(G[a - 1:b + 6]) + rnd(60, 4)
    site = [s for s in scan(G, target, (a, b)) if s.start == 66][0]   # seed = end of the planted block
    assert intarna.score_site(G, target, site) is not None


# ---- q-values ------------------------------------------------------------------------------

def test_qvalues_refuse_to_run_without_decoys():
    with pytest.raises(ValueError):
        qvalues([-5.0, -3.0], [], 0)
    assert qvalues([-5.0, -3.0], [], 3) == [0.0, 0.0]   # decoys existed, none reached these


def test_decoy_weights_for_unequal_numbers_of_decoys():
    # candidate A has 2 decoys (weight 0.5), candidate B has 4 (weight 0.25)
    q = qvalues([-10.0, -2.0], [-9.0, -9.0, -1.0], 0, decoy_weights=[0.5, 0.25, 0.25])
    # at -2: decoy sites <= -2: -9 (0.5) + -9 (0.25) = 0.75 expected false; 2 real -> 0.375
    assert q == [0.0, 0.375]


def make_regions():
    shared = rnd(120, 11) + revcomp(G[1:8]) + "A" + rnd(120, 12)
    return [Region("t.3utr", "t", "g", "three_prime_UTR", shared, "c", "+", [(0, len(shared))])]


def fake(guide, target, site, xi, cfg):
    real = guide == G
    return {"intarna_E": -20.0 if real else -8.0, "rnahybrid_mfe": -30.0 if real else -12.0,
            "rnahybrid_p": 0.1, "rnahybrid_pos": 1, "rnahybrid_at_site": 1}


def test_candidate_without_decoys_gets_no_q_value(monkeypatch):
    monkeypatch.setattr(score, "score_one", fake)
    monkeypatch.setattr(score, "check_tools", lambda cfg: None)
    decoys = make_decoys([("c1", G)], 40, (2, 8), 5)   # 40: with the +1 a lone site needs q floor 1/40
    messages = []
    cfg = score.ScoreConfig(calibrate_samples=0)
    _, sites, _, _ = score.run_scoring([("c1", G), ("c2", G)], decoys, make_regions(), (2, 8),
                                       None, cfg, ".", log=messages.append)
    by = {r["candidate"]: r for r in sites}
    assert by["c1"]["intarna_q"] != "" and by["c1"]["support"] >= 1
    assert by["c2"]["intarna_q"] == "" and by["c2"]["support"] == 0     # not "q = 0"
    assert any("no decoys" in m for m in messages)


def test_no_decoys_at_all_is_an_error(monkeypatch):
    monkeypatch.setattr(score, "score_one", fake)
    with pytest.raises(ValueError, match="needs decoys"):
        score.run_scoring([("c1", G)], [], make_regions(), (2, 8), None,
                          score.ScoreConfig(calibrate_samples=0), ".", log=lambda *_: None)


def test_unequal_decoy_counts_are_weighted(monkeypatch):
    monkeypatch.setattr(score, "score_one", fake)
    monkeypatch.setattr(score, "check_tools", lambda cfg: None)
    d1 = make_decoys([("c1", G)], 20, (2, 8), 5)
    d2 = make_decoys([("c2", G)], 2, (2, 8), 6)
    cfg = score.ScoreConfig(calibrate_samples=0)
    _, sites, _, curves = score.run_scoring([("c1", G), ("c2", G)], d1 + d2, make_regions(),
                                            (2, 8), None, cfg, ".", log=lambda *_: None)
    assert all(r["intarna_q"] != "" for r in sites)
    assert all(0 <= c["fdr"] <= 1 for c in curves)


def test_stage_two_error_keeps_the_rnahybrid_score():
    scores = {("c", 0, 0, 8, 7, "8mer"): {"rnahybrid_mfe": -20.0, "_error": "boom"}}
    for i in range(1, 40):
        scores[("c", 0, i, i + 8, i + 7, "8mer")] = {"rnahybrid_mfe": -10.0}
    score.summarise_errors(scores, lambda *_: None, 0.5)
    assert scores[("c", 0, 0, 8, 7, "8mer")] == {"rnahybrid_mfe": -20.0}


def test_best_site_class_is_the_strongest_among_equal_support():
    rows = [{"candidate": "c", "transcript": "t", "gene": "g", "site_class": cls, "support": 1,
             "intarna_E": e, "rnahybrid_mfe": -20.0, "intarna_q": 0.01, "rnahybrid_q": 0.01}
            for cls, e in (("7mer-A1", -30.0), ("8mer", -10.0), ("7mer-m8", -15.0))]
    out = score.summarise_targets(rows, score.ScoreConfig())
    assert out[0]["best_site_class"] == "8mer"


# ---- tool robustness -------------------------------------------------------------------

def test_a_hanging_tool_is_a_site_error_not_a_hang(monkeypatch):
    def hang(cmd, **kw):
        raise subprocess.TimeoutExpired(cmd, kw.get("timeout"))
    monkeypatch.setattr(intarna.subprocess, "run", hang)
    site = Site(214, 222, 221, "8mer")
    with pytest.raises(RuntimeError, match="timed out"):
        intarna.score_site(G, rnd(500, 8), site, timeout=1)


def test_calibration_rejects_non_finite_parameters():
    assert rnahybrid.parse_calibration("q 6 inf -inf") is None
    assert rnahybrid.parse_calibration("q 6 1.5 -0.2") is None
    assert rnahybrid.parse_calibration("q 6 1.5 0.2") == (1.5, 0.2)


def test_require_site_drops_a_duplex_elsewhere_in_the_window(monkeypatch):
    site = Site(100, 107, 107, "7mer-m8")   # seed on the target: 1-based 101..107
    target = "A" * 300

    def fake_run(start, span):
        out = f"t:300:q:22:-20.0:0.01:{start}:{'A' * span}: {' ' * span}:{' ' * span}:   " + chr(10)
        return type("P", (), {"returncode": 0, "stdout": out, "stderr": ""})()

    monkeypatch.setattr(rnahybrid.subprocess, "run",
                        lambda *a, **k: fake_run(150, 10))      # duplex far from the seed
    assert rnahybrid.score_site(G, target, site, flank=300) is None
    kept = rnahybrid.score_site(G, target, site, flank=300, require_site=False)
    assert kept["rnahybrid_at_site"] == 0


# ---- decoys that keep the seed ---------------------------------------------------------

def test_keep_seed_decoys_share_the_seed_and_differ_beyond_it():
    cands = [("a", G), ("b", "TAGCTTATCAGACTGATGTTGA")]
    decoys = make_decoys(cands, 8, (2, 8), 3, keep_seed=True)
    assert len(decoys) == 16
    for did, cid, seq in decoys:
        real = dict(cands)[cid]
        assert seq[:8] == real[:8] and len(seq) == len(real) and seq != real
        assert sorted(seq[8:]) == sorted(real[8:])           # same composition, shuffled
    assert len({s for _, c, s in decoys if c == "a"}) == 8   # distinct
    plain = make_decoys(cands, 8, (2, 8), 3)
    assert any(seed_of(s, (2, 8)) != seed_of(dict(cands)[c], (2, 8)) for _, c, s in plain)


# ---- validate ---------------------------------------------------------------------------

def test_a_single_gene_universe_does_not_inflate_the_auc():
    out, summary = validate([], [{"candidate": "c", "gene": "g"}])
    assert out[0]["scored"] == 0
    assert summary["auc"]["intarna"] is None


def test_validate_reports_percentile_detail():
    rows = [{"candidate": "c", "gene": f"g{i}", "transcript": f"g{i}.1", "best_intarna_E": str(-i),
             "best_rnahybrid_mfe": str(-i), "support": "0", "best_intarna_q": "", "best_rnahybrid_q": ""}
            for i in range(1, 101)]
    truth = [{"candidate": "c", "gene": "g100"}, {"candidate": "c", "gene": "g50"}]
    _, summary = validate(rows, truth)
    d = summary["detail"]["intarna"]
    assert d["top5"] == 0.5 and d["top10"] == 0.5          # g100 is rank 1, g50 is rank 51
    assert 0.2 < d["median_pct"] < 0.3


@pytest.fixture(autouse=True)
def _work_in_tmp(tmp_path, monkeypatch):
    """run_scoring writes checkpoints into its work directory; keep tests apart."""
    monkeypatch.chdir(tmp_path)


def _count_calls(monkeypatch):
    calls = {"n": 0}

    def counting(guide, target, site, xi, cfg):
        calls["n"] += len(cfg.predictors)
        return fake(guide, target, site, xi, cfg)

    monkeypatch.setattr(score, "score_one", counting)
    monkeypatch.setattr(score, "check_tools", lambda cfg: None)
    return calls


def test_scores_are_resumed_from_a_matching_checkpoint(monkeypatch):
    calls = _count_calls(monkeypatch)
    decoys = make_decoys([("c1", G)], 10, (2, 8), 5)
    cfg = score.ScoreConfig(calibrate_samples=0)
    msgs = []
    first = score.run_scoring([("c1", G)], decoys, make_regions(), (2, 8), None, cfg, "w",
                              log=msgs.append)
    n_first = calls["n"]
    assert n_first > 0
    msgs.clear()
    second = score.run_scoring([("c1", G)], decoys, make_regions(), (2, 8), None, cfg, "w",
                               log=msgs.append)
    assert calls["n"] == n_first                     # nothing was scored again
    assert any("resumed all" in m for m in msgs)
    assert [r["intarna_q"] for r in first[1]] == [r["intarna_q"] for r in second[1]]


def test_a_changed_setting_or_sequence_invalidates_the_checkpoint(monkeypatch):
    calls = _count_calls(monkeypatch)
    decoys = make_decoys([("c1", G)], 10, (2, 8), 5)
    base = score.ScoreConfig(calibrate_samples=0)
    score.run_scoring([("c1", G)], decoys, make_regions(), (2, 8), None, base, "w", log=lambda *_: None)
    n = calls["n"]
    score.run_scoring([("c1", G)], decoys, make_regions(), (2, 8), None,
                      score.ScoreConfig(calibrate_samples=0, temperature=25.0), "w", log=lambda *_: None)
    assert calls["n"] > n                            # new temperature: scored again
    n = calls["n"]
    other = make_decoys([("c1", G)], 10, (2, 8), 99)  # other decoys
    score.run_scoring([("c1", G)], other, make_regions(), (2, 8), None, base, "w", log=lambda *_: None)
    assert calls["n"] > n
    n = calls["n"]
    score.run_scoring([("c1", G)], other, make_regions(), (2, 8), None,
                      score.ScoreConfig(calibrate_samples=0, resume=False), "w", log=lambda *_: None)
    assert calls["n"] > n                            # resume off: always scored


def test_stage_one_checkpoint_survives_a_failure_in_stage_two(monkeypatch):
    state = {"calls": 0, "fail": True}

    def flaky(guide, target, site, xi, cfg):
        state["calls"] += 1
        if "intarna" in cfg.predictors and state["fail"]:
            raise KeyboardInterrupt            # simulates the run being killed in stage 2
        return fake(guide, target, site, xi, cfg)

    monkeypatch.setattr(score, "score_one", flaky)
    monkeypatch.setattr(score, "check_tools", lambda cfg: None)
    decoys = make_decoys([("c1", G)], 10, (2, 8), 5)
    cfg = score.ScoreConfig(calibrate_samples=0, intarna_prefilter=0.5)
    with pytest.raises(KeyboardInterrupt):
        score.run_scoring([("c1", G)], decoys, make_regions(), (2, 8), None, cfg, "w", log=lambda *_: None)
    stage1_calls = state["calls"]
    state["fail"] = False
    msgs = []
    score.run_scoring([("c1", G)], decoys, make_regions(), (2, 8), None, cfg, "w", log=msgs.append)
    assert any("resumed the RNAhybrid scores" in m for m in msgs)
    # stage 1 (one call per site) was not repeated: only IntaRNA calls were added
    assert state["calls"] - stage1_calls <= stage1_calls


# ---- pseudocount, GC-matched decoys, calibration window ----------------------------------

def test_pseudocount_keeps_q_above_zero():
    # 20 decoys per candidate -> pseudo weight 0.05; one real site can reach at best 0.05
    assert qvalues([-9.0], [-1.0], 20, plus_one=0.05) == [0.05]
    # with 10 real sites better than every decoy the floor is 0.005
    assert qvalues([-9.0] * 10, [-1.0], 20, plus_one=0.05) == [0.005] * 10
    # without the pseudocount the same sites get q = 0
    assert qvalues([-9.0] * 10, [-1.0], 20) == [0.0] * 10


def test_pseudocount_in_run_scoring(monkeypatch):
    monkeypatch.setattr(score, "score_one", fake)
    monkeypatch.setattr(score, "check_tools", lambda cfg: None)
    decoys = make_decoys([("c1", G)], 20, (2, 8), 5)
    with_pc = score.run_scoring([("c1", G)], decoys, make_regions(), (2, 8), None,
                                score.ScoreConfig(calibrate_samples=0), "a", log=lambda *_: None)[1]
    without = score.run_scoring([("c1", G)], decoys, make_regions(), (2, 8), None,
                                score.ScoreConfig(calibrate_samples=0, fdr_pseudocount=False),
                                "b", log=lambda *_: None)[1]
    assert without[0]["intarna_q"] == 0.0
    assert with_pc[0]["intarna_q"] == pytest.approx(1 / 20)      # one real site, 20 decoys


def test_match_seed_gc_decoys_have_the_real_seed_gc():
    cands = [("a", "TGAGGTAGTAGGTTGTATAGTT"), ("b", "CGCGGGATGGAGCAGTCTGGTAGC")]
    from smithhunter.targets.shuffle import gc_count
    decoys = make_decoys(cands, 15, (2, 8), 7, mode="match_seed_gc")
    assert len(decoys) == 30
    real = dict(cands)
    matched = 0
    for _, cid, seq in decoys:
        assert sorted(seq) == sorted(real[cid])                  # same composition
        matched += gc_count(seed_of(seq, (2, 8))) == gc_count(seed_of(real[cid], (2, 8)))
        assert seed_of(seq, (2, 8)) not in {seed_of(s, (2, 8)) for s in real.values()}
    assert matched >= 27                                         # almost all match exactly
    plain = make_decoys(cands, 15, (2, 8), 7)
    from collections import Counter
    spread = lambda ds: len({gc_count(seed_of(s, (2, 8))) for _, c, s in ds if c == "b"})
    assert spread(decoys) == 1 and spread(plain) > 1


def test_unknown_decoy_mode_is_rejected():
    with pytest.raises(ValueError, match="decoy mode"):
        make_decoys([("a", G)], 2, (2, 8), 1, mode="nope")


def test_calibration_uses_long_random_sequences_not_the_short_window(monkeypatch):
    seen = {}

    def fake_calibrate(guide, targets_fasta, seed, mean_len, sd_len, samples, binary):
        seen["mean"], seen["sd"] = mean_len, sd_len
        return (1.5, 0.2)

    monkeypatch.setattr(score.rnahybrid, "calibrate", fake_calibrate)
    monkeypatch.setattr(score, "score_one", fake)
    monkeypatch.setattr(score, "check_tools", lambda cfg: None)
    decoys = make_decoys([("c1", G)], 20, (2, 8), 5)
    score.run_scoring([("c1", G)], decoys, make_regions(), (2, 8), None,
                      score.ScoreConfig(calibrate_samples=50), "w", log=lambda *_: None)
    # short windows make the extreme-value fit unstable (negative xi), see score.py
    assert (seen["mean"], seen["sd"]) == score.CALIBRATION_LENGTH == (500, 300)
