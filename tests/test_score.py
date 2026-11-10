import random
import shutil

import pytest

from smithhunter.targets import intarna, rnahybrid, score
from smithhunter.targets.regions import Region, revcomp
from smithhunter.targets.shuffle import make_decoys
from smithhunter.targets.sites import Site, scan

GUIDE = "TAGCTTATCAGACTGATGTTGA"  # 22 nt, seed 2-8 AGCTTAT
M8 = revcomp(GUIDE[1:8])


def rnd(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


# ---- parsers, on output captured from IntaRNA 3.4.1 and RNAhybrid 2.1.2 ----------------

INTARNA_OUT = (
    "id1;id2;start1;end1;start2;end2;E;ED1;ED2;Pu1;Pu2;E_hybrid;hybridDP\n"
    "target;query;159;180;1;22;-28.9;5.34;2.07;0.000172637;0.0347825;-36.31;"
    "((((((((((((((((((((((&))))))))))))))))))))))\n")
RNAHYBRID_OUT = ("command_line:52:command_line:22:-41.5:0.000000:15:G                      C: "
                 "UCAACAUCAGUCUGAUAAGCUA : AGUUGUAGUCAGACUAUUCGAU :                        \n")
RNAHYBRID_NONE = "command_line:52:command_line:22:0.0:1.000000:71:G  C: AA :  UU :\n"


def test_parse_intarna_csv():
    r = intarna.parse_csv(INTARNA_OUT)
    assert r["intarna_E"] == -28.9 and r["intarna_ED_target"] == 5.34
    assert (r["intarna_t_start"], r["intarna_t_end"]) == (159, 180)
    assert (r["intarna_q_start"], r["intarna_q_end"]) == (1, 22)
    assert intarna.parse_csv("id1;id2;start1;end1;start2;end2;E\n") is None
    assert intarna.parse_csv("") is None


def test_parse_rnahybrid_compact_and_no_duplex():
    r = rnahybrid.parse_compact(RNAHYBRID_OUT)
    assert r == {"rnahybrid_mfe": -41.5, "rnahybrid_p": 0.0, "rnahybrid_pos": 15,
                 "rnahybrid_span": 24}
    assert rnahybrid.parse_compact(RNAHYBRID_NONE) is None
    assert rnahybrid.parse_compact("") is None


def test_parse_rnacalibrate():
    assert rnahybrid.parse_calibration("g 130 1.634881 0.296835\n") == (1.634881, 0.296835)
    assert rnahybrid.parse_calibration("") is None


# ---- geometry ----------------------------------------------------------------------------

@pytest.mark.parametrize("cls,extra,seedlen", [("8mer", "A", 7), ("7mer-m8", "C", 7)])
def test_window_pins_seed_on_the_site(cls, extra, seedlen):
    target = rnd(300, 1) + M8 + extra + rnd(300, 2)
    site = [s for s in scan(GUIDE, target) if s.site_class == cls][0]
    seq, lo, t_region, seed_t = intarna.window(target, site, len(GUIDE), flank=50)
    assert seq[seed_t[0] - 1:seed_t[1]] == M8            # seed on the target, 1-based inclusive
    assert seed_t[1] - seed_t[0] + 1 == seedlen
    assert t_region[0] <= seed_t[0] and seed_t[1] <= t_region[1]


def test_window_seed_for_six_mer_and_offset():
    target = rnd(200, 3) + "G" + revcomp(GUIDE[1:7]) + "C" + rnd(200, 4)
    site = [s for s in scan(GUIDE, target) if s.site_class == "6mer"][0]
    seq, lo, t_region, seed_t = intarna.window(target, site, len(GUIDE))
    assert seq[seed_t[0] - 1:seed_t[1]] == revcomp(GUIDE[1:7])
    off = rnd(200, 5) + revcomp(GUIDE[2:8]) + "CC" + rnd(200, 6)
    site = [s for s in scan(GUIDE, off) if s.site_class == "offset-6mer"][0]
    seq, lo, t_region, seed_t = intarna.window(off, site, len(GUIDE))
    assert seq[seed_t[0] - 1:seed_t[1]] == revcomp(GUIDE[2:8])
    assert intarna.seed_range("offset-6mer") == (3, 8)
    assert intarna.seed_range("seed4-10") == (4, 10)
    with pytest.raises(ValueError):
        intarna.seed_range("nonsense")


def test_window_clamped_at_region_ends():
    target = M8 + "A" + rnd(20, 7)   # site at the very start: less than guide length upstream
    site = scan(GUIDE, target)[0]
    seq, lo, t_region, seed_t = intarna.window(target, site, len(GUIDE))
    assert lo == 0 and t_region[0] == 1 and t_region[1] <= len(seq)


def test_command_lines():
    site = Site(214, 222, 221, "8mer")
    seq, lo, treg, seed_t = intarna.window(rnd(500, 8), site, len(GUIDE))
    cmd = intarna.build_command(GUIDE, seq, "8mer", treg, seed_t, 25.0)
    assert cmd[0] == "IntaRNA" and "U" in cmd[cmd.index("-q") + 1] and "T" not in cmd[cmd.index("-q") + 1]
    assert cmd[cmd.index("--seedQRange") + 1] == "2-8" and cmd[cmd.index("--seedBP") + 1] == "7"
    assert cmd[cmd.index("--temperature") + 1] == "25.0" and cmd[cmd.index("--threads") + 1] == "1"
    cmd = intarna.build_command(GUIDE, seq, "6mer", treg, seed_t)
    assert cmd[cmd.index("--seedQRange") + 1] == "2-7" and cmd[cmd.index("--seedBP") + 1] == "6"
    rh = rnahybrid.build_command(GUIDE, "ACGT" * 20, (2, 8), (1.5, 0.25))
    assert rh[rh.index("-f") + 1] == "2,8" and rh[rh.index("-d") + 1] == "1.5,0.25"
    assert "-s" not in rh
    assert "-s" in rnahybrid.build_command(GUIDE, "ACGT" * 20, (2, 7))


# ---- scoring and q-values, with the external tools replaced ----------------------------

def fake_score_one(guide, target, site, xi_theta, cfg):
    """Real guide binds well (-25 / -30); decoys bind worse (-12 / -15), with some spread."""
    real = guide == GUIDE
    jitter = (sum(map(ord, guide)) % 7) / 10
    out = {}
    if "intarna" in cfg.predictors:
        out.update(intarna_E=(-25.0 if real else -12.0 + jitter), intarna_ED_target=1.0)
    if "rnahybrid" in cfg.predictors:
        out.update(rnahybrid_mfe=(-30.0 if real else -15.0 + jitter), rnahybrid_p=0.01,
                   rnahybrid_pos=10, rnahybrid_at_site=1)
    return out


def setup_regions(decoys):
    """Two isoforms sharing a 3'UTR with the real 8mer site, plus a UTR carrying 8mer sites
    of three decoys (so decoys also have sites to score)."""
    shared = rnd(120, 11) + M8 + "A" + rnd(120, 12)
    planted = rnd(100, 13)
    for _, _, seq in decoys[:3]:
        planted += revcomp(seq[1:8]) + "A" + rnd(40, 14)
    return [
        Region("tx1.3utr", "tx1", "geneA", "three_prime_UTR", shared, "chr1", "+", [(0, len(shared))]),
        Region("tx2.3utr", "tx2", "geneA", "three_prime_UTR", shared, "chr1", "+", [(0, len(shared))]),
        Region("tx3.3utr", "tx3", "geneB", "three_prime_UTR", planted, "chr2", "+", [(0, len(planted))]),
    ]


@pytest.fixture
def scored(monkeypatch):
    monkeypatch.setattr(score, "score_one", fake_score_one)
    monkeypatch.setattr(score, "check_tools", lambda cfg: None)
    decoys = make_decoys([("c1", GUIDE)], 20, (2, 8), 5)
    regions = setup_regions(decoys)
    cfg = score.ScoreConfig(calibrate_samples=0, workers=1)
    result = score.run_scoring([("c1", GUIDE)], decoys, regions, (2, 8), None, cfg, ".",
                               log=lambda *_: None)
    return result, decoys, regions


def test_real_sites_scored_once_per_sequence_and_listed_per_isoform(scored):
    (classes, sites, targets, curves), _, _ = scored
    assert classes == ("8mer", "7mer-m8", "7mer-A1")          # weak classes not scored by default
    assert {r["transcript"] for r in sites} == {"tx1", "tx2"}  # decoy-only region has no real site
    assert len(sites) == 2 and all(r["site_class"] == "8mer" for r in sites)
    assert {r["gene"] for r in sites} == {"geneA"}
    assert targets[0]["gene"] == "geneA" or targets[0]["scored_sites"] == 1


def test_qvalues_use_decoys_and_support_counts_predictors(scored):
    (_, sites, _, curves), decoys, _ = scored
    r = sites[0]
    assert r["intarna_E"] == -25.0 and r["rnahybrid_mfe"] == -30.0
    assert r["intarna_q"] <= 0.05 and r["rnahybrid_q"] <= 0.05 and r["support"] == 2
    assert r["intarna_q_cand"] != ""
    assert {c["predictor"] for c in curves} == {"intarna", "rnahybrid"}
    assert all(0 <= c["fdr"] <= 1 for c in curves)


def test_decoys_that_score_better_give_high_q(monkeypatch):
    def reversed_scores(guide, target, site, xi, cfg):
        real = guide == GUIDE
        return {"intarna_E": -10.0 if real else -30.0, "rnahybrid_mfe": -10.0 if real else -30.0,
                "rnahybrid_p": 0.5, "rnahybrid_pos": 1, "rnahybrid_at_site": 1}
    monkeypatch.setattr(score, "score_one", reversed_scores)
    monkeypatch.setattr(score, "check_tools", lambda cfg: None)
    decoys = make_decoys([("c1", GUIDE)], 20, (2, 8), 5)
    cfg = score.ScoreConfig(calibrate_samples=0)
    _, sites, _, _ = score.run_scoring([("c1", GUIDE)], decoys, setup_regions(decoys), (2, 8),
                                       None, cfg, ".", log=lambda *_: None)
    assert all(r["intarna_q"] > 0.05 and r["support"] == 0 for r in sites)


def test_weak_sites_only_with_flag(monkeypatch):
    monkeypatch.setattr(score, "score_one", fake_score_one)
    monkeypatch.setattr(score, "check_tools", lambda cfg: None)
    target = rnd(100, 20) + "G" + revcomp(GUIDE[1:7]) + "C" + rnd(100, 21)   # a 6mer only
    regions = [Region("r", "t", "g", "three_prime_UTR", target, "c", "+", [(0, len(target))])]
    decoys = make_decoys([("c1", GUIDE)], 5, (2, 8), 1)
    off = score.ScoreConfig(calibrate_samples=0)
    on = score.ScoreConfig(calibrate_samples=0, score_weak_sites=True)
    quiet = lambda *_: None
    assert score.run_scoring([("c1", GUIDE)], decoys, regions, (2, 8), None, off, ".", log=quiet)[1] == []
    rows = score.run_scoring([("c1", GUIDE)], decoys, regions, (2, 8), None, on, ".", log=quiet)[1]
    assert [r["site_class"] for r in rows] == ["6mer"]


def test_single_predictor_and_dry_run(monkeypatch, capsys):
    monkeypatch.setattr(score, "score_one", fake_score_one)
    monkeypatch.setattr(score, "check_tools", lambda cfg: None)
    decoys = make_decoys([("c1", GUIDE)], 20, (2, 8), 5)
    regions = setup_regions(decoys)
    cfg = score.ScoreConfig(predictors=("intarna",), calibrate_samples=0)
    _, sites, _, curves = score.run_scoring([("c1", GUIDE)], decoys, regions, (2, 8), None, cfg,
                                            ".", log=lambda *_: None)
    assert all(r.get("rnahybrid_mfe") in (None, "") for r in sites) and sites[0]["support"] == 1
    assert {c["predictor"] for c in curves} == {"intarna"}
    out = score.run_scoring([("c1", GUIDE)], decoys, regions, (2, 8), None, cfg, ".", dry_run=True)
    assert out is None
    assert "sites to score" in capsys.readouterr().out


def test_write_scores_roundtrip(scored, tmp_path):
    (_, sites, targets, curves), _, _ = scored
    score.write_scores(sites, targets, curves, tmp_path / "s.tsv", tmp_path / "t.tsv",
                       tmp_path / "c.tsv")
    header = (tmp_path / "s.tsv").read_text().splitlines()[0].split("\t")
    assert "intarna_E" in header and "rnahybrid_q" in header and "support" in header
    assert len((tmp_path / "s.tsv").read_text().splitlines()) == len(sites) + 1


# ---- with the real tools (skipped where they are not installed) --------------------------

needs_tools = pytest.mark.skipif(
    not (shutil.which("IntaRNA") and shutil.which("RNAhybrid") and shutil.which("RNAcalibrate")),
    reason="IntaRNA, RNAhybrid and RNAcalibrate are not on PATH")


@needs_tools
def test_real_tools_score_planted_site_better_than_nothing(tmp_path):
    # a guide with a planted 8mer plus supplementary pairing, against random sequence
    target = rnd(200, 31) + revcomp(GUIDE[8:22]) + M8 + "A" + rnd(200, 32)
    site = [s for s in scan(GUIDE, target) if s.site_class == "8mer"][0]
    i = intarna.score_site(GUIDE, target, site)
    h = rnahybrid.score_site(GUIDE, target, site)
    assert i["intarna_E"] < -15 and h["rnahybrid_mfe"] < -25 and h["rnahybrid_at_site"] == 1
    assert i["intarna_t_end"] == site.t1 + 1                      # interaction ends at t1
    # a seed pinned where there is no complementarity finds nothing or something very weak
    fake = Site(100, 107, 107, "7mer-m8")
    f = intarna.score_site(GUIDE, target, fake)
    assert f is None or f["intarna_E"] > i["intarna_E"] + 10
    assert rnahybrid.score_site(GUIDE, target, fake) is None


@needs_tools
def test_real_calibration_and_scoring_pipeline(tmp_path):
    targets = tmp_path / "t.fa"
    targets.write_text("".join(f">r{i}\n{intarna.rna(rnd(300, 100 + i))}\n" for i in range(20)))
    xt = rnahybrid.calibrate(GUIDE, str(targets), (2, 8), 300, 50, samples=200)
    assert xt is not None and xt[0] > 0 and xt[1] > 0
    target = rnd(200, 41) + revcomp(GUIDE[8:22]) + M8 + "A" + rnd(200, 42)
    site = [s for s in scan(GUIDE, target) if s.site_class == "8mer"][0]   # the planted one
    h = rnahybrid.score_site(GUIDE, target, site, xt)
    assert h["rnahybrid_at_site"] == 1 and h["rnahybrid_p"] < 0.01


def test_two_stage_prefilter_runs_intarna_only_on_best_sites(monkeypatch):
    calls = {"intarna": 0, "rnahybrid": 0}

    def counting(guide, target, site, xi, cfg):
        for p in cfg.predictors:
            calls[p] += 1
        return fake_score_one(guide, target, site, xi, cfg)

    monkeypatch.setattr(score, "score_one", counting)
    monkeypatch.setattr(score, "check_tools", lambda cfg: None)
    decoys = make_decoys([("c1", GUIDE)], 20, (2, 8), 5)
    regions = setup_regions(decoys)
    full = score.ScoreConfig(calibrate_samples=0)
    score.run_scoring([("c1", GUIDE)], decoys, regions, (2, 8), None, full, ".", log=lambda *_: None)
    assert calls["intarna"] == calls["rnahybrid"]
    n_all = calls["rnahybrid"]
    calls.update(intarna=0, rnahybrid=0)
    two = score.ScoreConfig(calibrate_samples=0, intarna_prefilter=0.2)
    _, sites, _, _ = score.run_scoring([("c1", GUIDE)], decoys, regions, (2, 8), None, two, ".",
                                       log=lambda *_: None)
    assert calls["rnahybrid"] == n_all and 0 < calls["intarna"] < n_all
    # the real site has the best MFE of its class, so it passes and gets both predictors
    assert all(r["intarna_E"] == -25.0 and r["rnahybrid_mfe"] == -30.0 for r in sites)
    assert sites[0]["support"] == 2


def test_prefilter_cutoff_is_a_pooled_quantile_per_class():
    scores = {("a", 0, 0, 8, 7, "8mer"): {"rnahybrid_mfe": -30.0},
              ("b", 0, 0, 8, 7, "8mer"): {"rnahybrid_mfe": -10.0},
              ("c", 0, 0, 8, 7, "8mer"): {"rnahybrid_mfe": -20.0},
              ("d", 0, 0, 8, 7, "8mer"): {"rnahybrid_mfe": None},
              ("e", 0, 0, 8, 7, "6mer"): {"rnahybrid_mfe": -5.0}}
    cut = score.prefilter_cutoffs(scores, 0.34)
    assert cut == {"8mer": -30.0, "6mer": -5.0}


def test_a_failing_site_is_left_unscored_not_fatal(monkeypatch):
    state = {"n": 0}

    def flaky(guide, target, site, xi, cfg):
        state["n"] += 1
        if guide == GUIDE and state["n"] == 1:
            raise RuntimeError("IntaRNA failed (1): boom")
        return fake_score_one(guide, target, site, xi, cfg)

    monkeypatch.setattr(score, "score_one", flaky)
    monkeypatch.setattr(score, "check_tools", lambda cfg: None)
    decoys = make_decoys([("c1", GUIDE)], 20, (2, 8), 5)
    messages = []
    cfg = score.ScoreConfig(calibrate_samples=0, max_error_share=0.5)
    _, sites, _, _ = score.run_scoring([("c1", GUIDE)], decoys, setup_regions(decoys), (2, 8),
                                       None, cfg, ".", log=messages.append)
    assert any("failed and are left unscored" in m for m in messages)
    assert all(r.get("intarna_E") in (None, "") for r in sites)   # the only real site failed


def test_too_many_failures_stop_the_run(monkeypatch):
    def broken(guide, target, site, xi, cfg):
        raise RuntimeError("tool crashed")

    monkeypatch.setattr(score, "score_one", broken)
    monkeypatch.setattr(score, "check_tools", lambda cfg: None)
    decoys = make_decoys([("c1", GUIDE)], 20, (2, 8), 5)
    with pytest.raises(RuntimeError, match="sites failed"):
        score.run_scoring([("c1", GUIDE)], decoys, setup_regions(decoys), (2, 8), None,
                          score.ScoreConfig(calibrate_samples=0), ".", log=lambda *_: None)


def test_missing_binary_is_reported_before_scoring():
    cfg = score.ScoreConfig(intarna_bin="IntaRNA_not_installed")
    with pytest.raises(SystemExit, match="IntaRNA_not_installed"):
        score.check_tools(cfg)


def test_at_site_means_the_duplex_covers_the_seed(monkeypatch):
    """A seed-only site (3' half of the guide unpaired) is still on the site; a duplex that
    lies elsewhere in the window is not."""
    site = Site(100, 107, 107, "7mer-m8")   # seed on the target: 1-based 102..108
    target = "A" * 300

    def fake(start, span):
        out = f"t:300:q:22:-20.0:0.01:{start}:{'A' * span}: {' ' * span}:{' ' * span}:{' ' * 3}" + chr(10)
        return type("P", (), {"returncode": 0, "stdout": out, "stderr": ""})()

    for start, span, expected in ((95, 14, 1), (102, 7, 1), (115, 10, 0), (70, 20, 0)):
        monkeypatch.setattr(rnahybrid.subprocess, "run", lambda *a, _s=start, _p=span, **k: fake(_s, _p))
        # start is reported in the window; window starts at lo, so pass lo = 0 via a long flank
        res = rnahybrid.score_site(GUIDE, target, site, flank=300, require_site=False)
        assert res["rnahybrid_at_site"] == expected, (start, span)


@pytest.fixture(autouse=True)
def _work_in_tmp(tmp_path, monkeypatch):
    """run_scoring writes checkpoints into its work directory; keep tests apart."""
    monkeypatch.chdir(tmp_path)
