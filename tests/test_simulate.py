import gzip

import pytest

from smithhunter.evaluate import (evaluate, expected_column, interval, mode_report,
                                  origin_by_feature, overlaps, read_yaml_scalars, run_settings)
from smithhunter.seqio import read_fasta, read_tsv, write_tsv
from smithhunter.simulate import ADAPTER_R1, FLANK, Params, Simulator, revcomp, simulate


def small_params(tmp_path, **kw):
    return Params(outdir=str(tmp_path / "sim"), reads=3000, nuclear_length=600_000,
                  bacterium_length=20_000, srnas=5, **kw)


def test_simulation_writes_consistent_truth(tmp_path):
    sim = simulate(small_params(tmp_path))
    out = tmp_path / "sim"
    mito = dict(read_fasta(str(out / "data/mito.fasta")))["chrM"]
    features = {f["feature_id"]: f for f in read_tsv(str(out / "truth/features.tsv"))}

    # canonical sequences are the genome (or its reverse complement) at the stated place
    for f in features.values():
        if f["category"] in ("srna", "srna_close_pair", "srna_numt_1mm"):
            s, e = int(f["start"]) - 1, int(f["end"])
            seq = mito[s:e] if f["strand"] == "+" else revcomp(mito[s:e])
            assert seq == f["sequence"]
    origin = features["srna_origin"]
    assert origin["wraps_origin"] == "1"
    assert origin["sequence"] == mito[int(origin["start"]) - 1:] + mito[:int(origin["end"])]

    # the identical NUMT is in the nuclear genome, the 1-mismatch one is not
    nuclear = dict(read_fasta(str(out / "data/nuclear.fasta")))["chrN"]
    assert features["srna_numt_identical"]["sequence"] in nuclear
    assert features["srna_numt_1mm"]["sequence"] not in nuclear

    # FASTQ: reads carry the adapter after the insert, counts match the truth
    with gzip.open(out / "data/sim1_1.fastq.gz", "rt") as fh:
        lines = fh.read().splitlines()
    reads = lines[1::4]
    assert len(reads) == sum(s.counts[0] for s in sim.sources)
    assert sum(ADAPTER_R1[:10] in r for r in reads) > 0.9 * len(reads)
    config = (out / "config/config.yaml").read_text()
    assert "data/mito.gff3" in config and "ambiguous: report" in config


def test_srna_isoforms_follow_the_model(tmp_path):
    sim = Simulator(small_params(tmp_path))
    mito = sim.build_mito()
    src = sim.srna(mito, "x", "srna", 1000, 22, "+", [100.0], "pass")
    inserts = [sim.insert(src, 0) for _ in range(5000)]
    canonical = sum(i == src.sequence for i in inserts) / len(inserts)
    five_exact = sum(i.startswith(src.sequence[:10]) for i in inserts) / len(inserts)
    assert 0.35 < canonical < 0.6  # 0.90 (5') x 0.55 (3') x 0.9 (no tail)
    assert five_exact > 0.85
    assert src.template[FLANK:FLANK + 22] == src.sequence


def write_results(tmp_path, loci):
    res = tmp_path / "results"
    for sub in ("discovery", "reference", "origin"):
        (res / sub).mkdir(parents=True, exist_ok=True)
    write_tsv(str(res / "reference/contigs.tsv"),
              [{"ref": "mito|chrM", "genome": "mito", "contig": "chrM", "length": 1000,
                "circular": 1}], ["ref", "genome", "contig", "length", "circular"])
    write_tsv(str(res / "reference/genomes.tsv"), [{"name": "mito", "role": "focus"}],
              ["name", "role"])
    cols = ["locus_id", "genome", "contig", "start", "end", "strand", "span", "wraps_origin",
            "rep_sequence", "rep_start", "rep_end", "total_count", "five_score", "three_score",
            "annotation_class", "pass"]
    write_tsv(str(res / "discovery/loci.tsv"),
              [dict(zip(cols, l)) for l in loci], cols)
    return str(res)


def test_evaluate_classifies_outcomes(tmp_path):
    truth = tmp_path / "truth"
    truth.mkdir()
    cols = ["feature_id", "category", "genome", "contig", "start", "end", "strand",
            "wraps_origin", "length", "sequence", "expected", "total_reads"]
    feats = [
        ("a", "srna", 101, 122, "+", 0, "AAAA", "pass"),
        ("b", "srna", 301, 322, "+", 0, "", "pass"),     # merged with c
        ("c", "srna", 330, 351, "+", 0, "", "pass"),
        ("d", "srna", 501, 522, "+", 0, "", "pass"),     # missed
        ("e", "srna", 990, 12, "+", 1, "", "pass"),      # across the origin
    ]
    write_tsv(str(truth / "features.tsv"),
              [dict(zip(cols, (f, c, "mito", "chrM", s, e, st, w, 22, q, x, 10)))
               for f, c, s, e, st, w, q, x in feats], cols)
    results = write_results(tmp_path, [
        ("m1", "mito", "chrM", 100, 125, "+", 26, 0, "AAAA", 101, 122, 50, 1, 1, "NA", 1),
        ("m2", "mito", "chrM", 290, 360, "+", 71, 0, "", 301, 322, 50, 1, 1, "NA", 1),
        ("m3", "mito", "chrM", 101, 122, "-", 22, 0, "", 101, 122, 50, 1, 1, "NA", 1),
        ("m4", "mito", "chrM", 991, 10, "+", 20, 1, "", 991, 10, 50, 1, 1, "NA", 0),
    ])
    rows, false_pos, _ = evaluate(str(truth), results)
    observed = {r["feature_id"]: r["observed"] for r in rows}
    assert observed == {"a": "pass", "b": "merged", "c": "merged", "d": "missed", "e": "fail"}
    a = next(r for r in rows if r["feature_id"] == "a")
    assert a["rep_is_canonical"] == 1 and a["rep_five_offset"] == 0 and a["extra_span"] == 4
    assert [fp["locus_id"] for fp in false_pos] == ["m3"]  # antisense, no truth on that strand


def test_circular_overlap():
    assert overlaps((995, 1010), (0, 5), 1000, True)
    assert not overlaps((995, 1000), (0, 5), 1000, False)
    assert interval({"start": "991", "end": "10", "wraps_origin": "1"}, 1000, True) == (990, 1010)


# --- genome-choice modes ------------------------------------------------------------------

def hamming_hits(seq, genome, max_diff):
    """Smallest number of differences of seq against genome (either strand), up to max_diff."""
    best = None
    for q in (seq, revcomp(seq)):
        for i in range(len(genome) - len(q) + 1):
            d = 0
            for a, b in zip(q, genome[i:i + len(q)]):
                d += a != b
                if d > max_diff:
                    break
            else:
                best = d if best is None else min(best, d)
    return best


def test_numt_scenarios_and_expected_outcome_per_mode(tmp_path):
    simulate(small_params(tmp_path))
    out = tmp_path / "sim"
    features = {f["feature_id"]: f for f in read_tsv(str(out / "truth/features.tsv"))}
    mito = dict(read_fasta(str(out / "data/mito.fasta")))["chrM"]
    nuclear = dict(read_fasta(str(out / "data/nuclear.fasta")))["chrN"]

    def expected(fid):
        f = features[fid]
        return f["expected_best"], f["expected_severe"], f["expected_evil"]

    def copy_in_nuclear(fid):
        # the copied segment is the smallRNA with 150 nt on each side: search a window of it
        f = features[fid]
        s, e = int(f["start"]) - 1, int(f["end"])
        left = mito[s - 30:s]
        at = nuclear.index(left) + len(left)
        return nuclear[at:at + e - s]

    def differences(fid):
        return sum(a != b for a, b in zip(features[fid]["sequence"], copy_in_nuclear(fid)))

    assert differences("srna_numt_identical") == 0
    assert expected("srna_numt_identical") == ("missed", "missed", "missed")
    assert differences("srna_numt_1mm") == 1
    assert expected("srna_numt_1mm") == ("pass", "pass", "missed")
    assert differences("srna_numt_2mm") == 2
    assert expected("srna_numt_2mm") == ("pass", "pass", "pass")
    # both differences lie where every isoform covers them (5' shift <= 2, 3' shift >= -3)
    seq, copy = features["srna_numt_2mm"]["sequence"], copy_in_nuclear("srna_numt_2mm")
    assert [i for i in range(22) if seq[i] != copy[i]] == [7, 14]
    # no other place of the nuclear genome fits these smallRNAs within 2 mismatches better
    assert hamming_hits(features["srna_numt_2mm"]["sequence"], nuclear, 2) == 2
    assert hamming_hits(features["srna_numt_1mm"]["sequence"], nuclear, 2) == 1

    # the expressed nuclear copy: reads carry the nuclear base, the truth is the mito place
    f = features["srna_numt_nuclear"]
    assert f["genome"] == "mito" and f["group"] == "nuclear_numt"
    assert differences("srna_numt_nuclear") == 1
    assert expected("srna_numt_nuclear") == ("fail", "fail", "fail")

    # ordinary smallRNAs expect the same in every mode; `expected` is the mode of the config
    assert expected("srna_001") == ("pass", "pass", "pass")
    assert all(f["expected"] == f["expected_best"] for f in features.values())
    config = (out / "config/config.yaml").read_text()
    assert "assignment:\n  mode: best\n  evil_max_diff: 1\n" in config
    assumed = {r["setting"]: r["value"] for r in read_tsv(str(out / "truth/assignment.tsv"))}
    assert assumed == {"mode": "best", "evil_max_diff": "1", "ambiguous": "report",
                       "mismatches": "1"}


def test_expressed_nuclear_copy_reads_match_the_nuclear_genome(tmp_path):
    sim = Simulator(small_params(tmp_path))
    sim.run()
    nuclear = dict(read_fasta(str(tmp_path / "sim/data/nuclear.fasta")))["chrN"]
    src = next(s for s in sim.sources if s.id == "srna_numt_nuclear")
    (pos, base), = src.alt.items()
    read = src.sequence[:pos - FLANK] + base + src.sequence[pos - FLANK + 1:]
    assert read in nuclear and src.sequence not in nuclear


def test_simulate_writes_the_chosen_mode(tmp_path):
    simulate(small_params(tmp_path, mode="evil"))
    out = tmp_path / "sim"
    assert "assignment:\n  mode: evil\n  evil_max_diff: 1\n" in (out / "config/config.yaml").read_text()
    features = {f["feature_id"]: f for f in read_tsv(str(out / "truth/features.tsv"))}
    assert features["srna_numt_1mm"]["expected"] == "missed"
    assert features["srna_numt_2mm"]["expected"] == "pass"


def test_evil_without_tolerance_expects_what_severe_expects(tmp_path):
    simulate(small_params(tmp_path, mode="evil", evil_max_diff=0))
    features = {f["feature_id"]: f for f in read_tsv(str(tmp_path / "sim/truth/features.tsv"))}
    assert features["srna_numt_1mm"]["expected_evil"] == "pass"
    assert features["srna_numt_1mm"]["expected"] == "pass"


def test_split_policy_changes_only_the_best_mode_expectation(tmp_path):
    simulate(small_params(tmp_path, ambiguous="split", mode="severe"))
    f = {f["feature_id"]: f for f in read_tsv(str(tmp_path / "sim/truth/features.tsv"))}["srna_numt_identical"]
    assert (f["expected_best"], f["expected_severe"], f["expected_evil"], f["expected"]) == \
        ("pass", "missed", "missed", "missed")


def test_simulate_rejects_unknown_mode_and_tolerance(tmp_path):
    with pytest.raises(ValueError, match="mode must be"):
        simulate(small_params(tmp_path, mode="nasty"))
    with pytest.raises(ValueError, match="evil_max_diff"):
        simulate(small_params(tmp_path, mode="evil", evil_max_diff=2))


def test_mode_scenarios_leave_everything_else_untouched(tmp_path, monkeypatch):
    with_modes = Simulator(Params(outdir=str(tmp_path / "a"), reads=3000, nuclear_length=600_000,
                                  bacterium_length=20_000, srnas=16))
    with_modes.run()
    monkeypatch.setattr(Simulator, "add_mode_scenarios", lambda self, mito, nuclear, taken: None)
    without = Simulator(Params(outdir=str(tmp_path / "b"), reads=3000, nuclear_length=600_000,
                               bacterium_length=20_000, srnas=16))
    without.run()
    for name in ("data/mito.fasta", "data/mito.gff3", "data/nuclear.gff3", "truth/targets.tsv"):
        assert (tmp_path / "a" / name).read_text() == (tmp_path / "b" / name).read_text()
    assert len(read_tsv(str(tmp_path / "a/truth/targets.tsv"))) == 15
    old = {s.id: (s.start, s.sequence, s.counts) for s in without.sources}
    new = {s.id: (s.start, s.sequence, s.counts) for s in with_modes.sources}
    assert set(new) - set(old) == {"srna_numt_2mm", "srna_numt_nuclear"}
    assert all(new[k] == v for k, v in old.items())


def test_run_settings_reads_provenance_then_config(tmp_path):
    res = tmp_path / "results"
    (res / "report").mkdir(parents=True)
    assert run_settings(str(res))["mode"] == "best"  # nothing found: defaults
    (tmp_path / "config").mkdir()
    (tmp_path / "config/config.yaml").write_text(
        "samples: config/samples.tsv\nambiguous: split\n\n# comment\nassignment:\n"
        "  mode: severe   # tie discarded\n  evil_max_diff: 1\nmapping:\n  mismatches: 1\n")
    got = run_settings(str(res))
    assert (got["mode"], got["evil_max_diff"], got["ambiguous"]) == ("severe", 1, "split")
    # provenance.yaml holds the configuration the run really used (e.g. --config overrides)
    (res / "report/provenance.yaml").write_text(
        "# header\nsoftware:\n  smithhunter: \"0.1\"\nsamples:\n  - \"s1\"\nconfig:\n"
        "  ambiguous: \"report\"\n  assignment:\n    mode: \"evil\"\n    evil_max_diff: 1\n"
        "  mapping:\n    mismatches: 1\n")
    got = run_settings(str(res))
    assert (got["mode"], got["evil_max_diff"], got["ambiguous"]) == ("evil", 1, "report")
    assert got["source"].endswith("provenance.yaml")
    values = read_yaml_scalars(str(res / "report/provenance.yaml"))
    assert values["config.mapping.mismatches"] == "1" and values["software.smithhunter"] == "0.1"


def test_expected_column_per_mode():
    assert expected_column({"mode": "best", "evil_max_diff": 0}) == "expected_best"
    assert expected_column({"mode": "severe", "evil_max_diff": 1}) == "expected_severe"
    assert expected_column({"mode": "evil", "evil_max_diff": 1}) == "expected_evil"
    assert expected_column({"mode": "evil", "evil_max_diff": 0}) == "expected_severe"


def mode_truth(tmp_path):
    truth = tmp_path / "truth"
    truth.mkdir()
    cols = ["feature_id", "category", "genome", "contig", "start", "end", "strand",
            "wraps_origin", "length", "sequence", "expected", "expected_best",
            "expected_severe", "expected_evil", "total_reads"]
    a, b = "ACGTTGCATCGATCGGATCCAT", "TTGACCAGTAGGCATTACGGCA"
    feats = [("kept", "srna", 101, 122, a, "pass", "pass", "pass", "pass"),
             ("numt", "srna_numt_1mm", 301, 322, b, "pass", "pass", "pass", "missed")]
    write_tsv(str(truth / "features.tsv"),
              [dict(zip(cols, (f, c, "mito", "chrM", s, e, "+", 0, 22, q, x, xb, xs, xe, 30)))
               for f, c, s, e, q, x, xb, xs, xe in feats], cols)
    return truth, a, b


def test_evaluate_uses_the_expected_outcome_of_the_run_mode(tmp_path):
    truth, a, b = mode_truth(tmp_path)
    results = write_results(tmp_path, [
        ("m1", "mito", "chrM", 101, 122, "+", 22, 0, a, 101, 122, 50, 1, 1, "NA", 1)])
    by_mode = {}
    for mode in ("best", "severe", "evil"):
        rows, _, _ = evaluate(str(truth), results, {"mode": mode, "evil_max_diff": 1})
        by_mode[mode] = {r["feature_id"]: (r["expected"], r["observed"], r["as_expected"])
                         for r in rows}
    # the locus of `numt` is gone: right under evil, wrong under best and severe
    assert by_mode["evil"]["numt"] == ("missed", "missed", 1)
    assert by_mode["best"]["numt"] == by_mode["severe"]["numt"] == ("pass", "missed", 0)
    assert all(v["kept"] == ("pass", "pass", 1) for v in by_mode.values())
    # evil without tolerance is judged as severe
    rows, _, _ = evaluate(str(truth), results, {"mode": "evil", "evil_max_diff": 0})
    assert {r["feature_id"]: r["expected"] for r in rows}["numt"] == "pass"


def test_evil_rejected_reads_are_followed_to_their_smallrna(tmp_path):
    truth, a, b = mode_truth(tmp_path)
    results = write_results(tmp_path, [
        ("m1", "mito", "chrM", 101, 122, "+", 22, 0, a, 101, 122, 50, 1, 1, "NA", 1)])
    reads = [("u1", a, "assigned", (20, 5)),
             ("u2", a[:-2], "assigned", (3, 0)),               # shorter isoform
             ("u3", b, "evil_rejected", (10, 8)),
             ("u4", b[1:] + "A", "evil_rejected", (2, 0)),      # 5' shift and a tail
             ("u5", b[:9] + ("A" if b[9] != "A" else "C") + b[10:], "ambiguous", (1, 0)),
             ("u6", "GGGGGGGGCCCCCCCCGGGGGG", "evil_rejected", (4, 0)),  # no simulated smallRNA
             ("u7", "CCCCCCCCGGGGGGGGCCCCCC", "unmapped", (7, 7))]
    with gzip.open(f"{results}/origin/read_origin.tsv.gz", "wt") as fh:
        fh.write("read_id\tclass\tgenomes\n")
        fh.writelines(f"{rid}\t{cls}\tmito\n" for rid, _, cls, _ in reads)
    (tmp_path / "results/collapsed").mkdir()
    with gzip.open(f"{results}/collapsed/counts.tsv.gz", "wt") as fh:
        fh.write("read_id\tsequence\tlength\ts1\ts2\n")
        fh.writelines(f"{rid}\t{seq}\t{len(seq)}\t{c[0]}\t{c[1]}\n" for rid, seq, _, c in reads)
    write_tsv(f"{results}/origin/origin_summary.tsv",
              [{"class": "assigned", "genomes": "mito", "s1": 23, "s2": 5},
               {"class": "evil_rejected", "genomes": "mito+nuclear", "s1": 16, "s2": 8}],
              ["class", "genomes", "s1", "s2"])
    settings = {"mode": "evil", "evil_max_diff": 1, "ambiguous": "report", "source": "test"}
    rows, _, features = evaluate(str(truth), results, settings)
    origin = {o["feature_id"]: o for o in origin_by_feature(features, results)}
    assert (origin["kept"]["assigned"], origin["kept"]["evil_rejected"]) == (28, 0)
    assert (origin["numt"]["evil_rejected"], origin["numt"]["ambiguous"],
            origin["numt"]["matched_reads"]) == (20, 1, 21)
    text = mode_report(rows, list(origin.values()), settings, str(truth), results)
    assert "Genome-choice mode of the run: evil (evil_max_diff 1)" in text
    assert "evil_rejected reads: 24 in the run, 20 from simulated smallRNAs, 4 from other" in text
    assert "smallRNAs without a passing locus after evil rejection: numt" in text
