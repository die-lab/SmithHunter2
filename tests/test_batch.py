"""Batch scoring (RNAhybrid multi-FASTA, IntaRNA multi-target / multi-query, chunked workers):
parsers, grouping, fallbacks, and exact identity with the one-process-per-site code.

Tests that need the real programs run only when they are on PATH.
"""
import pickle
import random
import shutil
import subprocess
from types import SimpleNamespace

import pytest

from smithhunter.targets import intarna, rnahybrid, score
from smithhunter.targets.regions import Region, revcomp
from smithhunter.targets.sites import Site, scan

GUIDE = "TAGCTTATCAGACTGATGTTGA"  # 22 nt, seed 2-8 AGCTTAT
M8 = revcomp(GUIDE[1:8])
HAVE_RNAHYBRID = shutil.which("RNAhybrid") is not None
HAVE_INTARNA = shutil.which("IntaRNA") is not None


def rnd(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


# ---- RNAhybrid: parser of a multi-target output ------------------------------------------

LINE = "{n}:52:q:22:{mfe}:{p}:{pos}:G                      C: UCAACAUCAGUCUGAUAAGCUA : AGUUGUAGUCAGACUAUUCGAU :  "


def test_parse_compact_batch_per_target():
    text = "\n".join([LINE.format(n="w0", mfe="-41.5", p="0.000000", pos=15),
                      "w1:30:q:22:0.0:1.000000:0::::",          # no duplex under the constraint
                      LINE.format(n="w3", mfe="-20.1", p="0.5", pos=3),
                      "garbage line without fields",
                      ""])
    got = rnahybrid.parse_compact_batch(text)
    assert got["w0"] == rnahybrid.parse_compact(LINE.format(n="w0", mfe="-41.5", p="0.000000", pos=15))
    assert got["w0"]["rnahybrid_mfe"] == -41.5 and got["w0"]["rnahybrid_span"] == 24
    assert got["w1"] is None and "w1" in got                       # known, but no score
    assert got["w3"]["rnahybrid_pos"] == 3
    assert "w2" not in got                                          # absent: NOT 'no duplex'
    assert rnahybrid.parse_compact_batch("") == {}


def test_parse_compact_batch_keeps_the_first_line_of_a_target():
    text = LINE.format(n="w0", mfe="-10.0", p="0.1", pos=1) + "\n" + \
        LINE.format(n="w0", mfe="-30.0", p="0.1", pos=2)
    assert rnahybrid.parse_compact_batch(text)["w0"]["rnahybrid_mfe"] == -10.0


def test_batch_command_line():
    cmd = rnahybrid.build_batch_command(GUIDE, (2, 8), 120, "/x/t.fa", (1.5, 0.25))
    assert cmd[cmd.index("-t") + 1] == "/x/t.fa" and cmd[cmd.index("-f") + 1] == "2,8"
    assert cmd[cmd.index("-m") + 1] == "121" and cmd[cmd.index("-d") + 1] == "1.5,0.25"
    assert cmd[-1] == rnahybrid.rna(GUIDE) and "-c" in cmd and cmd[cmd.index("-b") + 1] == "1"
    assert "-s" in rnahybrid.build_batch_command(GUIDE, (2, 7), 10, "t.fa")


# ---- RNAhybrid: fake tool, grouping, fallback --------------------------------------------

class FakeRNAhybrid:
    """Stands in for ``intarna.run_tool``. Reads the target FASTA of a batch call (or the
    sequence of a single call) and answers per window; windows containing 'N' crash the
    tool (non-zero exit); windows ending in 'GGGG' are silently dropped from batch output."""

    def __init__(self, drop=False):
        self.calls = []
        self.drop = drop

    @staticmethod
    def answer(name, seq, pos=3):
        if "A" * 12 in seq:
            return f"{name}:{len(seq)}:q:22:0.0:1.000000:0::::"
        # the duplex covers the whole window: always 'at site'
        return (f"{name}:{len(seq)}:q:22:{-10.0 - len(seq) / 10:.1f}:0.25:{pos}:" + "A" * len(seq)
                + ": " + "U" * 0 + ": : ")

    def __call__(self, cmd, timeout=None):
        self.calls.append(cmd)
        if "-t" in cmd:
            path = cmd[cmd.index("-t") + 1]
            recs = open(path).read().split(">")[1:]
            names_seqs = [r.strip().split("\n") for r in recs]
        else:
            names_seqs = [["command_line", cmd[-2]]]
        if any("N" in s for _, s in names_seqs):
            return SimpleNamespace(returncode=1, stdout="", stderr="boom")
        out = [self.answer(n, s) for n, s in names_seqs if not (self.drop and "-t" in cmd and s.endswith("GGGG"))]
        return SimpleNamespace(returncode=0, stdout="\n".join(out) + "\n", stderr="")


def sites_in(target):
    return [(target, s) for s in scan(GUIDE, target, (2, 8), ("8mer", "7mer-m8"))]


def test_score_sites_matches_score_site_on_a_fake_tool(monkeypatch, tmp_path):
    fake = FakeRNAhybrid()
    monkeypatch.setattr(intarna, "run_tool", fake)
    items = []
    for k in range(12):
        items += sites_in(rnd(60 + 7 * k, k) + M8 + "A" + rnd(50, 100 + k))
    assert len(items) >= 12
    batch = rnahybrid.score_sites(GUIDE, items, tmpdir=str(tmp_path))
    n_batch_calls = len(fake.calls)
    single = [rnahybrid.score_site(GUIDE, t, s) for t, s in items]
    assert batch == single and n_batch_calls == 1                 # ONE process for the chunk
    assert list(tmp_path.iterdir()) == []                        # temporary FASTA removed


def test_a_poisoned_window_costs_only_its_own_site(monkeypatch, tmp_path):
    fake = FakeRNAhybrid()
    monkeypatch.setattr(intarna, "run_tool", fake)
    items = []
    for k in range(40):
        items += sites_in(rnd(80, k) + M8 + "A" + rnd(60, 100 + k))
    bad = 17
    t_bad, s_bad = items[bad]
    p = s_bad.end + 4
    items[bad] = (t_bad[:p] + "N" + t_bad[p + 1:], s_bad)       # inside its window: the tool crashes
    got = rnahybrid.score_sites(GUIDE, items, tmpdir=str(tmp_path), min_split=4)
    assert isinstance(got[bad], RuntimeError) and "boom" in str(got[bad])
    ok = [i for i in range(len(items)) if i != bad]
    assert all(isinstance(got[i], dict) for i in ok)
    ref = [rnahybrid.score_site(GUIDE, *items[i]) for i in ok]
    assert [got[i] for i in ok] == ref
    assert len(fake.calls) < 2 * len(items)                       # bisection, not 1 call per site


def test_a_target_missing_from_the_output_is_rescored_not_dropped(monkeypatch, tmp_path):
    fake = FakeRNAhybrid(drop=True)
    monkeypatch.setattr(intarna, "run_tool", fake)
    t = rnd(60, 1) + M8 + "A" + "GGGG"                              # window ends in GGGG
    items = sites_in(t) + sites_in(rnd(60, 2) + M8 + "A" + rnd(30, 3))
    got = rnahybrid.score_sites(GUIDE, items, tmpdir=str(tmp_path))
    # batch output lacked the first window; the single call (which the fake answers) fills it
    assert got[0] == rnahybrid.score_site(GUIDE, *items[0]) and isinstance(got[0], dict)


def test_score_sites_refuses_mixed_seed_ranges(tmp_path):
    t = rnd(60, 1) + M8 + "A" + rnd(40, 2)
    s8 = scan(GUIDE, t, (2, 8), ("8mer",))[0]
    s6 = Site(s8.start, s8.end, s8.t1, "offset-6mer")
    with pytest.raises(ValueError):
        rnahybrid.score_sites(GUIDE, [(t, s8), (t, s6)], tmpdir=str(tmp_path))
    assert rnahybrid.score_sites(GUIDE, []) == []


def test_missing_binary_is_an_error_of_every_site(tmp_path):
    items = sites_in(rnd(60, 1) + M8 + "A" + rnd(40, 2))
    got = rnahybrid.score_sites(GUIDE, items, binary="RNAhybrid_does_not_exist",
                                tmpdir=str(tmp_path))
    assert all(isinstance(g, Exception) for g in got)


# ---- RNAhybrid: identity with the real tool on random windows ----------------------------

def random_windows(n, seed):
    """n (target, site) with sites of all classes in random regions of 40-400 nt; some sites
    sit near a region end (clamped windows)."""
    rng = random.Random(seed)
    items = []
    while len(items) < n:
        k = rng.randint(0, 5)
        body = rnd(rng.randint(5, 200), rng.randint(0, 10 ** 9))
        site = rng.choice([M8 + "A", M8 + "C", revcomp(GUIDE[1:7]) + "A", revcomp(GUIDE[1:7]) + "G",
                           revcomp(GUIDE[2:8]) + "CC"])
        t = (rnd(rng.randint(0, 150), rng.randint(0, 10 ** 9)) if k else "") + site + body
        items += [(t, s) for s in scan(GUIDE, t, (2, 8))]
    return items[:n]


@pytest.mark.skipif(not HAVE_RNAHYBRID, reason="RNAhybrid not on PATH")
@pytest.mark.parametrize("xi_theta", [None, (1.63, 0.297)])
def test_rnahybrid_batch_is_identical_to_one_process_per_site(tmp_path, xi_theta):
    n = 2000 if xi_theta is None else 400
    items = random_windows(n, 7)
    groups = {}
    for i, (t, s) in enumerate(items):
        groups.setdefault(rnahybrid.seed_range(s.site_class), []).append(i)
    assert len(groups) == 3                                         # (2,8), (2,7), (3,8)
    for idx in groups.values():
        got = rnahybrid.score_sites(GUIDE, [items[i] for i in idx], xi_theta,
                                    tmpdir=str(tmp_path))
        for i, g in zip(idx, got):
            ref = rnahybrid.score_site(GUIDE, *items[i], xi_theta)
            assert g == ref, (i, g, ref)
    # require_site=False keeps the duplexes elsewhere in the window, also identical
    idx = groups[(2, 8)][:200]
    got = rnahybrid.score_sites(GUIDE, [items[i] for i in idx], xi_theta, require_site=False,
                                tmpdir=str(tmp_path))
    assert got == [rnahybrid.score_site(GUIDE, *items[i], xi_theta, require_site=False)
                   for i in idx]


# ---- IntaRNA: parser, geometry, fake tool -------------------------------------------------

def test_parse_csv_batch():
    head = "id1;id2;start1;end1;start2;end2;E;ED1;ED2;Pu1;Pu2;E_hybrid;hybridDP\n"
    row = "{t};{q};159;180;1;22;-28.9;5.34;2.07;0.0001;0.03;-36.31;((&))\n"
    got = intarna.parse_csv_batch(head + row.format(t="t0", q="q0") + row.format(t="t1", q="q0")
                                  + row.format(t="t0", q="q1"))
    assert set(got) == {("t0", "q0"), ("t1", "q0"), ("t0", "q1")}
    assert got[("t0", "q0")] == intarna.parse_csv(head + row.format(t="x", q="y"))
    assert intarna.parse_csv_batch(head) == {} and intarna.parse_csv_batch("") == {}


def test_window_coords_agree_with_window():
    for n, t1, cls in [(500, 250, "8mer"), (500, 20, "7mer-m8"), (500, 498, "7mer-A1"),
                       (40, 20, "6mer"), (300, 0 + 30, "offset-6mer")]:
        target = rnd(n, 5)
        site = Site(t1 - 7, t1 + 1, t1, cls)
        seq, lo, treg, seedt = intarna.window(target, site, 22, 100)
        lo2, hi2, treg2, seedt2 = intarna.window_coords(n, site, 22, 100)
        assert (lo, treg, seedt) == (lo2, treg2, seedt2) and seq == target[lo2:hi2]


class FakeIntaRNA:
    """Answers every (target, query) pair of a call with a row derived from the sequences;
    a target containing 'N' crashes the tool."""

    HEAD = "id1;id2;start1;end1;start2;end2;E;ED1;ED2;Pu1;Pu2;E_hybrid;hybridDP\n"

    def __init__(self):
        self.calls = []

    @staticmethod
    def read(arg):
        if arg.endswith(".fa"):
            recs = open(arg).read().split(">")[1:]
            return [r.strip().split("\n") for r in recs]
        return [["command_line", arg]]

    def __call__(self, cmd, timeout=None):
        self.calls.append(cmd)
        ts, qs = self.read(cmd[cmd.index("-t") + 1]), self.read(cmd[cmd.index("-q") + 1])
        if any("N" in s for _, s in ts):
            return SimpleNamespace(returncode=1, stdout="", stderr="boom")
        rows = []
        for tn, ts_ in ts:
            for qn, qs_ in qs:
                if "AAAAA" in ts_ + qs_:      # no interaction under the constraints: no row
                    continue
                e = -(len(ts_) + 2 * len(qs_)) / 10
                rows.append(f"{tn};{qn};5;20;1;16;{e};1.0;2.0;0.1;0.2;{e - 1};((&))")
        return SimpleNamespace(returncode=0, stdout=self.HEAD + "\n".join(rows) + "\n", stderr="")


def test_score_pairs_cross_product_and_single_pair_equivalence(monkeypatch, tmp_path):
    fake = FakeIntaRNA()
    monkeypatch.setattr(intarna, "run_tool", fake)
    queries = [rnd(21 + k, k) for k in range(3)]
    windows = [rnd(120 + 5 * k, 50 + k) for k in range(4)]
    windows[2] = "AAAAAC" + windows[2]                              # no interaction for any query
    got = intarna.score_pairs(queries, windows, "8mer", (101, 133), (125, 131), tmpdir=str(tmp_path))
    assert len(fake.calls) == 1 and len(got) == 12
    for (wi, qi), r in got.items():
        ref = intarna.parse_csv(fake(intarna.build_command(queries[qi], windows[wi], "8mer",
                                                           (101, 133), (125, 131))).stdout)
        assert r == ref
    assert got[(2, 0)] is None
    assert list(tmp_path.iterdir()) == []


def test_score_pairs_halves_a_failing_call(monkeypatch, tmp_path):
    fake = FakeIntaRNA()
    monkeypatch.setattr(intarna, "run_tool", fake)
    windows = [rnd(120, k) for k in range(20)]
    windows[13] = windows[13][:50] + "N" + windows[13][51:]
    got = intarna.score_pairs([rnd(22, 1)], windows, "8mer", (101, 133), (125, 131),
                              tmpdir=str(tmp_path))
    assert isinstance(got[(13, 0)], RuntimeError) and "boom" in str(got[(13, 0)])
    assert sum(isinstance(v, Exception) for v in got.values()) == 1 and len(got) == 20
    assert len(fake.calls) <= 2 + 2 * 5 + 1                         # halving, not 20 single calls
    # the same, many queries and one window
    fake.calls.clear()
    qs = [rnd(22, k) for k in range(9)]
    qs[4] = qs[4][:5] + "N" + qs[4][6:]
    got = intarna.score_pairs(qs, [rnd(100, 3)], "7mer-A1", (101, 133), (125, 130),
                              tmpdir=str(tmp_path))
    assert len(got) == 9 and not any(isinstance(v, Exception) for v in got.values())
    assert intarna.score_pairs([], [rnd(50, 1)], "8mer", (1, 2), (1, 2)) == {}


# ---- IntaRNA: identity with the real tool -------------------------------------------------

@pytest.mark.skipif(not HAVE_INTARNA, reason="IntaRNA not on PATH")
def test_intarna_batch_is_identical_to_one_process_per_site(tmp_path):
    guides = [GUIDE, GUIDE[:8] + rnd(14, 1), GUIDE[:8] + rnd(14, 2)]       # same length, same seed
    target = rnd(260, 9) + M8 + "A" + rnd(260, 10)
    site = scan(GUIDE, target, (2, 8), ("8mer",))[0]
    # one window, three small RNAs (shared target accessibility)
    seq, lo, treg, seedt = intarna.window(target, site, len(GUIDE), 100)
    got = intarna.score_pairs(guides, [seq], "8mer", treg, seedt, tmpdir=str(tmp_path))
    for qi, g in enumerate(guides):
        ref = intarna.score_site(g, target, site, 100)
        r = dict(got[(0, qi)])
        r["intarna_t_start"] += lo
        r["intarna_t_end"] += lo
        assert r == ref
    # several windows of one layout, one small RNA, one of them with no room to interact
    wins, refs = [], []
    for k in range(5):
        t = rnd(260 + 9 * k, 20 + k) + M8 + "A" + rnd(250 + 3 * k, 40 + k)
        s = scan(GUIDE, t, (2, 8), ("8mer",))[0]
        w, lo2, treg2, seedt2 = intarna.window(t, s, len(GUIDE), 100)
        assert (treg2, seedt2) == (treg, seedt)
        wins.append(w)
        refs.append((intarna.score_site(GUIDE, t, s, 100), lo2))
    got = intarna.score_pairs([GUIDE], wins, "8mer", treg, seedt, tmpdir=str(tmp_path))
    for wi, (ref, lo2) in enumerate(refs):
        r = got[(wi, 0)]
        if r:
            r["intarna_t_start"] += lo2
            r["intarna_t_end"] += lo2
        assert r == ref


# ---- score.py: grouping of calls, chunks, planning, whole-run identity --------------------

def make_state(guides, unique):
    return {"guides": guides, "unique": unique}


def test_plan_intarna_calls_shares_windows_and_separates_layouts():
    t1 = rnd(300, 1) + M8 + "A" + rnd(300, 2)
    t2 = M8 + "A" + rnd(300, 3)                       # site at the region start: clamped window
    g2 = GUIDE[:8] + rnd(14, 4)                       # same seed and length: shares windows
    guides = {"a": GUIDE, "b": g2}
    unique = [t1, t2]
    jobs = []
    for gid in guides:
        for si, seq in enumerate(unique):
            for s in scan(guides[gid], seq, (2, 8), ("8mer",)):
                jobs.append((gid, si, s.start, s.end, s.t1, s.site_class))
    assert len(jobs) == 4
    cfg = score.ScoreConfig(intarna_flank=100)
    calls = score.plan_intarna_calls(jobs, range(len(jobs)), make_state(guides, unique), cfg, 50)
    by_call = sorted(sorted(jobs[i][0] + str(jobs[i][1]) for i in c) for c in calls)
    # same window of both small RNAs = one call; the clamped ones too (identical layout there)
    assert by_call == [["a0", "b0"], ["a1", "b1"]]
    # different guide lengths do not share a window
    guides2 = {"a": GUIDE, "b": GUIDE[:8] + rnd(15, 4)}
    jobs2 = [(g, 0, s.start, s.end, s.t1, s.site_class) for g in guides2
             for s in scan(guides2[g], unique[0], (2, 8), ("8mer",))]
    calls = score.plan_intarna_calls(jobs2, range(len(jobs2)), make_state(guides2, unique), cfg, 50)
    assert sorted(len(c) for c in calls) == [1, 1]
    # distinct layouts of one small RNA never share a call (no padding)
    calls = score.plan_intarna_calls([j for j in jobs if j[0] == "a"], [0, 1],
                                     make_state(guides, unique), cfg, 50)
    assert sorted(len(c) for c in calls) == [1, 1]


def test_plan_intarna_calls_respects_the_size_limit():
    guides = {"a": GUIDE}
    unique = [rnd(300, k) + M8 + "A" + rnd(300, 100 + k) for k in range(23)]
    jobs = [("a", si, s.start, s.end, s.t1, s.site_class)
            for si, seq in enumerate(unique) for s in scan(GUIDE, seq, (2, 8), ("8mer",))]
    cfg = score.ScoreConfig(intarna_flank=100)
    calls = score.plan_intarna_calls(jobs, range(len(jobs)), make_state(guides, unique), cfg, 10)
    assert sorted(len(c) for c in calls) == [3, 10, 10] and \
        sorted(i for c in calls for i in c) == list(range(len(jobs)))


def test_chunks_group_by_small_rna_and_seed_range():
    jobs = [("a", 0, 0, 8, 7, "8mer"), ("b", 0, 0, 8, 7, "8mer"), ("a", 1, 0, 7, 6, "7mer-A1"),
            ("a", 2, 0, 8, 7, "7mer-m8"), ("a", 3, 0, 8, 7, "8mer"), ("a", 4, 0, 6, 5, "offset-6mer")]
    chunks = score._chunks(jobs, 2)
    assert sorted(sorted(c) for c in chunks) == [[0, 3], [1], [2], [4], [5]]
    for c in chunks:
        assert len({(jobs[i][0], intarna.seed_range(jobs[i][5])) for i in c}) == 1 and len(c) <= 2


def test_plan_jobs_is_the_plain_scan_serial_or_parallel():
    guides = {f"g{i}": GUIDE[:8] + rnd(14, i) for i in range(5)}
    guides["real"] = GUIDE
    unique = [rnd(400, k) + M8 + "A" + rnd(300, 50 + k) for k in range(30)]
    classes = score.scored_classes(("8mer", "7mer-m8", "7mer-A1", "6mer", "offset-6mer"),
                                   score.ScoreConfig())
    ref = [(gid, si, s.start, s.end, s.t1, s.site_class)
           for gid, g in guides.items() for si, seq in enumerate(unique)
           for s in scan(g, seq, (2, 8), classes)]
    assert len(ref) > 30
    assert score.plan_jobs(guides, unique, (2, 8), classes) == ref
    assert score.plan_jobs(guides, unique, (2, 8), classes, score.ScoreConfig(workers=2)) == ref
    # the integers of the keys are shared objects, not one copy per job
    jobs = score.plan_jobs(guides, unique, (2, 8), classes)
    assert len({id(j[2]) for j in jobs}) < len(jobs)


def test_parallel_map_order_and_progress():
    log = []
    cfg = score.ScoreConfig(workers=1)
    out = score.parallel_map(lambda x: x * 2, list(range(100)), cfg, {}, log=log.append,
                             sizes=[50] * 100)               # 5000 'sites' in 100 items
    assert out == [x * 2 for x in range(100)]
    assert 15 <= len(log) <= 21 and log[-1].startswith("47")


def tiny_regions(n):
    shared = rnd(120, 11) + M8 + "A" + rnd(120, 12)
    return [Region(id=f"r{i}", transcript=f"t{i}", gene=f"g{i}", region="three_prime_UTR",
                   contig="c", strand="+", blocks=[(1, len(s))], sequence=s)
            for i, s in enumerate([shared] + [rnd(200, 70 + k) + M8 + "A" + rnd(150, 90 + k)
                                              for k in range(n - 1)])]


@pytest.mark.skipif(not (HAVE_RNAHYBRID and HAVE_INTARNA), reason="tools not on PATH")
@pytest.mark.parametrize("prefilter,predictors,workers",
                         [(0.6, ("intarna", "rnahybrid"), 1),     # two-stage
                          (1.0, ("intarna", "rnahybrid"), 1),     # one pass, both inside a chunk
                          (1.0, ("rnahybrid",), 2),               # worker processes
                          (1.0, ("intarna",), 1)])
def test_run_scoring_batch_equals_per_site(tmp_path, prefilter, predictors, workers):
    from smithhunter.targets.shuffle import make_decoys
    regions = tiny_regions(6)
    cands = [("c1", GUIDE)]
    decoys = make_decoys(cands, 3, mode="keep_seed")             # decoys share windows with c1
    res = {}
    for name, rb, ib in [("per_site", 0, 0), ("batch", 2000, 50)]:
        cfg = score.ScoreConfig(predictors=predictors, calibrate_samples=0, workers=workers,
                                intarna_prefilter=prefilter, rnahybrid_batch=rb, intarna_batch=ib,
                                intarna_flank=60, resume=False, timeout=60)
        wd = tmp_path / name
        out = score.run_scoring(cands, decoys, regions, (2, 8), None, cfg, str(wd),
                                log=lambda m: None)
        with open(wd / "checkpoint_scores.pkl", "rb") as fh:
            ck = pickle.load(fh)["scores"]
        res[name] = (out[1], out[2], out[3], ck)
    assert res["per_site"][3] is not None and len(res["per_site"][3]) >= 8
    assert res["batch"][3] == res["per_site"][3]
    assert list(res["batch"][3].items()) == list(res["per_site"][3].items())    # same order too
    assert res["batch"][0] == res["per_site"][0] and res["batch"][1] == res["per_site"][1]
    assert res["batch"][2] == res["per_site"][2]


def test_batch_worker_reports_a_failing_tool_as_site_errors(monkeypatch, tmp_path):
    """A missing RNAhybrid binary: every site of the chunk gets the same error text as in
    the per-site path, and nothing else."""
    guides = {"a": GUIDE}
    unique = [rnd(80, 1) + M8 + "A" + rnd(80, 2)]
    jobs = [("a", 0, s.start, s.end, s.t1, s.site_class) for s in scan(GUIDE, unique[0], (2, 8), ("8mer",))]
    cfg = score.ScoreConfig(predictors=("rnahybrid",), rnahybrid_bin="RNAhybrid_does_not_exist",
                            calibrate_samples=0)
    score._init({"guides": guides, "unique": unique, "xi": {}, "cfg": cfg, "tmpdir": str(tmp_path)})
    batch = score._work_batch(jobs)
    single = [score._run_job(j) for j in jobs]
    assert batch == single and all("_error" in r for r in batch)
