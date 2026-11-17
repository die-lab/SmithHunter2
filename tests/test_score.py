"""Module B step 2: IntaRNA / RNAhybrid wrappers, q-values and support.

The parsers are tested on stored outputs of the real tools (tests/data/targets), the
pipeline with fake predictors, so nothing here needs the tools. The last test runs them
and is skipped when they are not installed.
"""

import json
import shutil
from pathlib import Path

import pytest

from smithhunter.targets import intarna, rnahybrid
from smithhunter.targets.fdr import qvalues
from smithhunter.targets.pipeline import run_sites
from smithhunter.targets.regions import Region, revcomp
from statistics import mean, pstdev

from smithhunter.targets.score import (Job, add_qvalues, collect_jobs, run_scoring, standardise,
                                       scored_classes)
from smithhunter.targets.sites import Site, scan

DATA = Path(__file__).parent / "data" / "targets"
CASE = json.loads((DATA / "site_8mer.json").read_text())
GUIDE, TARGET = CASE["guide"], CASE["target"]
SITE = Site(*CASE["site"])


def test_seed_pairing_per_class():
    guide = "TGAGGTAGTAGGTTGTATAGTT"
    h28 = revcomp(guide[1:8])
    for target, cls, expect in (
            ("CCCC" + h28 + "ACCCC", "8mer", (2, 8, 4, 10)),
            ("CCCC" + h28 + "CCCCC", "7mer-m8", (2, 8, 4, 10)),
            ("CCCCG" + h28[1:] + "ACCCC", "7mer-A1", (2, 7, 5, 10)),
            ("CCCCG" + h28[1:] + "CCCCC", "6mer", (2, 7, 5, 10)),
            ("CCCC" + h28[:6] + "GCCCC", "offset-6mer", (3, 8, 4, 9))):
        site, = scan(guide, target)
        assert site.site_class == cls
        g_from, g_to, t_from, t_to = intarna.seed_pairing(site)
        assert (g_from, g_to, t_from, t_to) == expect
        # the paired stretch is the reverse complement of the guide stretch
        assert target[t_from:t_to + 1] == revcomp(guide[g_from - 1:g_to])


def test_seed_pairing_other_seed_range():
    guide = "TGAGGTAGTAGGTTGTATAGTT"
    target = "CCCC" + revcomp(guide[3:10]) + "CCCC"
    site, = scan(guide, target, (4, 10))
    assert intarna.seed_pairing(site, (4, 10)) == (4, 10, 4, 10)


def test_intarna_command_pins_seed_and_clips_to_the_window():
    cmd, offset = intarna.command(GUIDE, TARGET, SITE)
    opt = {cmd[i]: cmd[i + 1] for i in range(1, len(cmd) - 1) if cmd[i].startswith("-")}
    lo, hi = intarna.window(TARGET, SITE, len(GUIDE))
    # the window: guide length + 15 nt upstream of the seed match, 15 nt past t1
    assert lo == SITE.t1 - 7 - len(GUIDE) - 15 and hi == min(len(TARGET), SITE.t1 + 1 + 15)
    assert offset == lo and opt["-t"] == TARGET[lo:hi] and "--tRegion" not in cmd
    assert opt["--seedQRange"] == "2-8" and opt["--seedBP"] == "7" and "--seedNoGU" in cmd
    assert opt["--outMaxE"] == "999" and opt["--seedMaxE"] == "999"   # every site scored
    a, b = (int(x) for x in opt["--seedTRange"].split("-"))
    assert opt["-t"][a - 1:b] == revcomp(GUIDE[1:8])
    whole, no_offset = intarna.command(GUIDE, TARGET, SITE, temperature=25, clip=False)
    assert no_offset == 0 and whole[4] == TARGET and whole[-2:] == ["--temperature", "25"]
    assert whole[whole.index("--tRegion") + 1] == f"{lo + 1}-{hi}"
    a, b = (int(x) for x in whole[whole.index("--seedTRange") + 1].split("-"))
    assert TARGET[a - 1:b] == revcomp(GUIDE[1:8])


def test_intarna_parser_on_stored_output():
    offset = SITE.start - CASE["intarna_context"]   # the fixture was run on a clipped target
    got = intarna.parse_csv((DATA / "intarna_8mer.csv").read_text(), offset)
    assert got["intarna_E"] == -6.51 and got["intarna_E_hybrid"] == -10.9
    assert got["intarna_ED_target"] == 3.1 and got["intarna_ED_query"] == 1.29
    assert got["intarna_target"] == "124-150" and got["intarna_guide"] == "1-20"
    # total energy = hybridisation + accessibility penalties
    assert got["intarna_E"] == pytest.approx(-10.9 + 3.1 + 1.29)
    # the paired target stretch ends at t1 (1-based) and covers the seed match
    assert int(got["intarna_target"].split("-")[1]) == SITE.t1 + 1


def test_slice_accessibility():
    table = [" #i$\tl=1\t2\t3\n"] + [
        "\t".join([str(i)] + [f"{i}.{length}" if length <= i else "NA" for length in (1, 2, 3)])
        + "\t\n" for i in range(1, 7)]
    got = intarna.slice_accessibility(table, 2, 5).splitlines()
    assert got[0] == table[0].rstrip("\n")
    # positions 3-5 become 1-3; lengths reaching before the cut are NA
    assert got[1:] == ["1\t3.1\tNA\tNA\t", "2\t4.1\t4.2\tNA\t", "3\t5.1\t5.2\t5.3\t"]
    # the whole table comes back unchanged
    assert intarna.slice_accessibility(table, 0, 6) == "".join(table)
    cmd, _ = intarna.command(GUIDE, TARGET, SITE, accessibility="x.acc")
    assert cmd[cmd.index("--tAcc") + 1] == "E" and cmd[cmd.index("--tAccFile") + 1] == "x.acc"
    acc = intarna.accessibility_command(GUIDE, TARGET, "x.acc", temperature=25)
    assert "tAcc:x.acc" in acc and acc[-2:] == ["--temperature", "25"]


def test_intarna_parser_without_interaction():
    header = (DATA / "intarna_8mer.csv").read_text().splitlines()[0]
    assert intarna.parse_csv(header + "\n") == {}
    assert intarna.parse_csv("") == {}


def test_rnahybrid_window_and_command():
    win, offset = rnahybrid.window(TARGET, SITE, len(GUIDE))
    assert offset == CASE["window_offset"]
    assert revcomp(GUIDE[1:8]) in win and win == TARGET[offset:offset + len(win)]
    cmd = rnahybrid.command(GUIDE, win, SITE, params=(1.7, 0.1))
    assert cmd[cmd.index("-f") + 1] == "2,8" and cmd[cmd.index("-d") + 1] == "1.7,0.1"
    assert cmd[-2:] == [win, GUIDE]
    six = Site(SITE.start + 1, SITE.end - 1, SITE.t1, "6mer")
    assert "2,7" in rnahybrid.command(GUIDE, win, six)
    assert "3utr_human" in rnahybrid.command(GUIDE, win, six)  # no calibration: a preset


def test_rnahybrid_parser_on_stored_output():
    text = (DATA / "rnahybrid_8mer.txt").read_text()
    got = rnahybrid.parse_compact(text, SITE, CASE["window_offset"])
    assert got["rnahybrid_mfe"] == -16.3 and got["rnahybrid_p"] == 0.003338
    assert got["rnahybrid_duplex"].count("|") == 3
    assert got["rnahybrid_pairs"] == 17 and got["rnahybrid_on_site"] == 1
    assert got["rnahybrid_target"] == "124-150"  # same stretch as IntaRNA finds
    assert rnahybrid.parse_compact(text, SITE, CASE["window_offset"],
                                   with_pvalue=False)["rnahybrid_p"] == ""
    assert rnahybrid.parse_compact("", SITE, 0) == {}
    # a duplex elsewhere in the window is flagged
    moved = Site(SITE.start - 3, SITE.end - 3, SITE.t1 - 3, "8mer")
    assert rnahybrid.parse_compact(text, moved, CASE["window_offset"])["rnahybrid_on_site"] == 0


def test_rnahybrid_duplex_pairs():
    fields = (DATA / "rnahybrid_8mer.txt").read_text().strip().split(":")
    pairs = rnahybrid.duplex_pairs(fields[7:11], int(fields[6]), int(fields[3]))
    win, _ = rnahybrid.window(TARGET, SITE, len(GUIDE))
    wc = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G"), ("G", "T"), ("T", "G")}
    assert all((win[t - 1], GUIDE[q - 1]) in wc for t, q in pairs)
    assert {q for _, q in pairs} >= set(range(2, 9))  # the helix constraint


def test_rnacalibrate_parser_drops_failed_fits():
    got = rnahybrid.parse_calibrate((DATA / "rnacalibrate.txt").read_text())
    assert got == {"g": (1.613649, 0.140387), "bad": (1.718316, 0.092426)}


def test_calibration_windows():
    seqs = ["ACGT" * 50, "TTGCA" * 8, "ACG"]
    wins = rnahybrid.calibration_windows(seqs, 60, 200, random_seed=3)
    assert len(wins) == 200 and wins == rnahybrid.calibration_windows(seqs, 60, 200, 3)
    assert all(len(w) == 60 or w == seqs[1] for w in wins)
    assert rnahybrid.calibration_windows(["ACG"], 60, 10) == []


def job(gid, cls="8mer", start=0, decoy=False, cand="c1"):
    return Job(gid, cand, decoy, "ACGT", 0, start, start + 8, start + 7, cls)


def test_standardise_leaves_one_decoy_out():
    decoys = [job(f"c1__d{i:02d}", start=i, decoy=True) for i in range(25)]
    real = job("c1", start=99)
    values = [float(-5 - (i % 5)) for i in range(25)]          # -5 ... -9
    scores = {j: {"intarna_E": v} for j, v in zip(decoys, values)}
    scores[real] = {"intarna_E": -20.0}
    standardise(decoys + [real], scores, "intarna", "intarna_E")
    m, sd = mean(values), pstdev(values)
    assert scores[real]["intarna_z"] == pytest.approx((-20 - m) / sd, abs=1e-3)
    rest = values[1:]
    assert scores[decoys[0]]["intarna_z"] == pytest.approx(
        (values[0] - mean(rest)) / pstdev(rest), abs=1e-3)
    # fewer than MIN_DECOY_SCORES decoy scores: no z at all
    few = {j: {"intarna_E": v} for j, v in zip(decoys[:10], values)}
    standardise(decoys[:10], few, "intarna", "intarna_E")
    assert all("intarna_z" not in v for v in few.values())


def test_leave_one_out_keeps_the_null_calibrated():
    """Real sites of no interest drawn from the decoy distribution: with the decoys
    standardising themselves they would pass more often than the decoys do."""
    import random
    rng = random.Random(5)
    false_calls = calls = 0
    for rep in range(40):
        jobs, scores = [], {}
        for c in range(10):
            for d in range(21):
                j = job(f"c{c}__d{d}", start=d, decoy=True, cand=f"c{c}")
                jobs.append(j)
                scores[j] = {"intarna_E": rng.gauss(-8, 2)}
            j = job(f"c{c}", start=50, cand=f"c{c}")
            jobs.append(j)
            scores[j] = {"intarna_E": rng.gauss(-8, 2)}
        add_qvalues(jobs, scores, {"intarna": "intarna_E"}, n_decoys=21, fdr=0.1)
        hits = sum(1 for j in jobs if not j.decoy and scores[j]["intarna_q"] <= 0.1)
        false_calls += hits
        calls += hits
    assert false_calls <= 4   # all null: q <= 0.1 should almost never happen


def test_standardised_qvalues_per_class_candidate_and_support():
    """An AT-rich candidate (weak energies) and a GC-rich one (strong energies), 20
    decoys with one site each: on raw energies the GC-rich decoys would hide the AT-rich
    site that clearly beats its own decoys; standardised, it stands out and the GC-rich
    site that is average for its decoys does not."""
    spread = [-1.0, 0.0, 1.0, 0.0, 0.0] * 4 + [0.0]
    at_decoys = [job(f"c1__d{i:02d}", start=10 + i, decoy=True) for i in range(21)]
    gc_decoys = [job(f"c2__d{i:02d}", start=40 + i, decoy=True, cand="c2") for i in range(21)]
    at_real, gc_real = job("c1", start=0), job("c2", start=2, cand="c2")
    seven = job("c1", "7mer-m8", start=90)
    energies = {at_real: -9.0, gc_real: -15.0, seven: -7.0,
                **{j: -5.0 + d for j, d in zip(at_decoys, spread)},
                **{j: -15.0 + d for j, d in zip(gc_decoys, spread)}}
    jobs = list(energies)
    scores = {j: {"intarna_E": e} for j, e in energies.items()}
    curves = add_qvalues(jobs, scores, {"intarna": "intarna_E"}, n_decoys=21, fdr=0.1)
    assert scores[at_real]["intarna_z"] < -5 and abs(scores[gc_real]["intarna_z"]) < 0.01
    # no decoy reaches the AT-rich site: q = (0 + 1) / 21 / 1
    assert scores[at_real]["intarna_q"] == pytest.approx(1 / 21, abs=1e-6)
    assert scores[at_real]["support"] == 1
    assert scores[gc_real]["intarna_q"] > 0.5 and scores[gc_real]["support"] == 0
    # raw energies pooled: 21 GC-rich decoys at or below -9, q = (21 + 1) / 21 / 2
    assert qvalues([-9.0, -15.0], [energies[j] for j in at_decoys + gc_decoys], 21,
                   pseudo=1)[0] == pytest.approx(22 / 42, abs=1e-6)
    assert scores[at_real]["intarna_q_candidate"] == pytest.approx(1 / 21, abs=1e-6)
    # classes are separate: the 7mer has no decoy site in its class
    assert scores[seven]["intarna_q"] == pytest.approx(1 / 21, abs=1e-6)
    assert "support" not in scores[at_decoys[0]]
    top = [c for c in curves if c["site_class"] == "8mer"]
    assert top[0]["real_sites"] == 1 and top[0]["decoy_sites"] == 0
    assert all(c["score"] == "intarna_z" for c in top)


def test_weak_classes_scored_only_on_request():
    classes = ("8mer", "7mer-m8", "7mer-A1", "6mer", "offset-6mer")
    assert scored_classes(classes, False) == classes[:3]
    assert scored_classes(classes, True) == classes


def fake_pipeline(score_weak=False):
    guide = "TGAGGTAGTAGGTTGTATAGTT"
    h28 = revcomp(guide[1:8])
    decoy_seq = "TCCGCCGGTTGATAGTTGAGTA"
    decoy_site = revcomp(decoy_seq[1:8]) + "A"   # an 8mer of the decoys
    utr = ("CCGCC" * 6 + h28 + "A" + "CCGCC" * 6 + "G" + h28[1:] + "CCGCC"   # 8mer, 6mer
           + "TTTTT" + decoy_site + "TTTTT")
    regions = [Region("t1|u", "t1", "g1", "three_prime_UTR", utr),
               Region("t2|u", "t2", "g1", "three_prime_UTR", utr),     # isoform, same UTR
               Region("t3|u", "t3", "g2", "three_prime_UTR", "GGCGG" * 8 + h28 + "CGG")]
    candidates = [("c1", guide)]
    decoys = [(f"c1__decoy{i:03d}", "c1", decoy_seq) for i in range(1, 26)]
    decoys.append(("c1__decoy006", "c1", "TTTTTTTTTTTTTTTTTTTTTT"))
    _, site_rows, _, _ = run_sites(candidates, decoys, regions)
    site_rows = [{k: str(v) for k, v in r.items()} for r in site_rows]
    calls = []

    def make(unique):
        def predict_all(jobs, threads):
            calls.extend(jobs)
            return [{"intarna_E": -5.0 - 0.1 * int(j.guide_id[-2:]) if j.decoy
                     else -15.0 - j.target}
                    for j in jobs]
        return predict_all
    out = run_scoring(candidates, decoys, regions, site_rows, {"intarna": make},
                      {"intarna": "intarna_E"}, score_weak=score_weak)
    return out, calls


def test_run_scoring_shares_scores_between_isoforms():
    (rows, targets, decoy_rows, curves, columns), calls = fake_pipeline()
    assert columns[-3:] == ["intarna_q", "intarna_q_candidate", "support"]
    # unique sequences are scored once: two real sites (8mer in UTR 1, 7mer-m8 in UTR 3)
    assert sorted((j.target, j.site_class) for j in calls if not j.decoy) == \
        [(0, "8mer"), (1, "7mer-m8")]
    by = {(r["transcript"], r["site_class"]): r for r in rows}
    assert by[("t1", "8mer")]["intarna_E"] == by[("t2", "8mer")]["intarna_E"] == -15.0
    assert by[("t3", "7mer-m8")]["intarna_E"] == -16.0
    assert by[("t1", "8mer")]["intarna_q"] < 0.05 and by[("t1", "8mer")]["support"] == 1
    # the 6mer is listed but not scored
    assert by[("t1", "6mer")]["intarna_E"] == "" and by[("t1", "6mer")]["support"] == ""
    assert all(r["decoy"].startswith("c1__decoy") for r in decoy_rows)
    t1 = next(t for t in targets if t["transcript"] == "t1")
    assert t1["sites"] == 2 and t1["best_class"] == "8mer" and t1["best_intarna_E"] == -15.0
    assert t1["max_support"] == 1 and t1["sites_supported"] == 1


def test_run_scoring_weak_sites_on_request():
    (rows, *_), calls = fake_pipeline(score_weak=True)
    assert any(j.site_class == "6mer" for j in calls)
    assert next(r for r in rows if r["site_class"] == "6mer")["intarna_E"] != ""


def test_collect_jobs_ignores_decoys_of_other_candidates():
    jobs = collect_jobs([("c1", GUIDE)], [("x__decoy001", "x", GUIDE)], [TARGET])
    assert [j.guide_id for j in jobs] == ["c1"]


@pytest.mark.skipif(not (shutil.which("IntaRNA") and shutil.which("RNAhybrid")
                         and shutil.which("RNAcalibrate")),
                    reason="IntaRNA / RNAhybrid not installed")
def test_real_tools_agree_with_stored_outputs():
    from smithhunter.targets.score import calibrate, intarna_predictor, rnahybrid_predictor
    j = Job("g", "g", False, GUIDE, 0, SITE.start, SITE.end, SITE.t1, SITE.site_class)
    got, = intarna_predictor([TARGET], (2, 8), None)([j], 1)
    assert got["intarna_target"].endswith("-150") and got["intarna_E"] < -5
    # the sliced accessibility table gives exactly what IntaRNA computes on the whole region
    cmd, offset = intarna.command(GUIDE, TARGET, SITE, clip=False)
    import subprocess
    whole = intarna.parse_csv(subprocess.run(cmd, capture_output=True, text=True).stdout, offset)
    assert got == whole
    # a region longer than IntaRNA's accessibility window (150 nt), site in the middle
    long = TARGET[::-1] * 2 + TARGET + TARGET[::-1] * 2
    site = Site(*(x + 2 * len(TARGET) for x in CASE["site"][:3]), SITE.site_class)
    assert site in scan(GUIDE, long)
    jl = Job("g", "g", False, GUIDE, 0, site.start, site.end, site.t1, site.site_class)
    sliced, = intarna_predictor([long], (2, 8), None)([jl], 1)
    cmd, offset = intarna.command(GUIDE, long, site, clip=False)
    assert sliced == intarna.parse_csv(
        subprocess.run(cmd, capture_output=True, text=True).stdout, offset)
    params = calibrate({"g": GUIDE}, [TARGET * 20], 300, 1, 1)
    assert set(params) == {"g"}
    hyb, = rnahybrid_predictor([TARGET], (2, 8), params)([j], 1)
    assert hyb["rnahybrid_mfe"] == pytest.approx(-16.3, abs=0.05)
    assert hyb["rnahybrid_on_site"] == 1 and 0 <= hyb["rnahybrid_p"] <= 1
    # a weak site still gets an IntaRNA score (cut-offs lifted)
    weak = "C" * 40 + "G" + revcomp(GUIDE[1:7]) + "A" + "C" * 40
    ws = next(x for x in scan(GUIDE, weak) if x.site_class == "7mer-A1")
    jw = Job("g", "g", False, GUIDE, 0, ws.start, ws.end, ws.t1, ws.site_class)
    assert intarna_predictor([weak], (2, 8), None)([jw], 1)[0]["intarna_E"] > -5
    from smithhunter.targets.score import run
    assert run([], [], [], str(DATA / "empty_sites.tsv"), ["intarna", "rnahybrid"],
               log=lambda m: None)[0] == []
