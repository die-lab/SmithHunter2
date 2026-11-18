import os
import random
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from smithhunter import precursor as pc
from smithhunter.cli import main
from smithhunter.seqio import read_tsv

DATA = Path(__file__).parent / "data"
LET7 = ("TGAGGTAGTAGGTTGTATAGTTTTAGGGTCACACCCACCACTGGGAGATAACTATACAATCTACTGTCTTTC")


# ------------------------------------------------------------------ coordinates

def test_to_interval_is_one_based_inclusive():
    assert pc.to_interval(61, 94, 16000) == (60, 94)
    assert pc.to_interval(1, 1, 100) == (0, 1)
    # across the origin: loci.tsv writes end < start
    assert pc.to_interval(15992, 13, 16000) == (15991, 16013)
    with pytest.raises(ValueError):
        pc.to_interval(0, 10, 100)


def test_windows_plus_strand():
    contig = "".join(random.Random(1).choice("ACGT") for _ in range(500))
    three, five = pc.make_windows(100, 122, "+", 500, False, (15, 50))
    assert (three.name, three.start, three.end) == ("3p", 50, 137)
    assert (five.name, five.start, five.end) == ("5p", 85, 172)
    assert (three.cand_from, three.cand_to) == (50, 72)
    assert (five.cand_from, five.cand_to) == (15, 37)
    for w in (three, five):
        assert not w.wraps and not w.clipped and w.length == 22 + 65
        seq = pc.window_sequence(contig, w, "+")
        assert seq == contig[w.start:w.end]
        assert seq[w.cand_from:w.cand_to] == contig[100:122]


def test_windows_minus_strand_are_reverse_complemented():
    contig = "".join(random.Random(2).choice("ACGT") for _ in range(500))
    candidate = pc.reverse_complement(contig[100:122])
    three, five = pc.make_windows(100, 122, "-", 500, False, (15, 50))
    # 5' of the RNA is at higher forward coordinates: the long flank of 3p is on the right
    assert (three.start, three.end) == (85, 172)
    assert (five.start, five.end) == (50, 137)
    assert (three.cand_from, three.cand_to) == (50, 72)
    assert (five.cand_from, five.cand_to) == (15, 37)
    for w in (three, five):
        seq = pc.window_sequence(contig, w, "-")
        assert seq == pc.reverse_complement(contig[w.start:w.end])
        assert seq[w.cand_from:w.cand_to] == candidate


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("start,end", [(5, 27), (980, 1002), (995, 1017), (0, 22), (978, 1000)])
def test_windows_wrap_the_origin_of_circular_contigs(strand, start, end):
    contig = "".join(random.Random(3).choice("ACGT") for _ in range(1000))
    doubled = contig * 3
    candidate = doubled[start:end]
    if strand == "-":
        candidate = pc.reverse_complement(candidate)
    for w in pc.make_windows(start, end, strand, 1000, True, (15, 50)):
        assert 0 <= w.start < 1000 and w.length == 87 and not w.clipped
        assert w.wraps == (w.end > 1000)
        seq = pc.window_sequence(contig, w, strand)
        assert len(seq) == 87
        assert seq[w.cand_from:w.cand_to] == candidate
        expected = doubled[w.start:w.end]
        assert seq == (expected if strand == "+" else pc.reverse_complement(expected))
    assert any(w.wraps for w in pc.make_windows(start, end, strand, 1000, True, (15, 50)))


def test_wrapping_window_coordinates():
    three, five = pc.make_windows(5, 27, "+", 1000, True, (15, 50))
    assert (three.start, three.end, three.wraps) == (955, 1042, True)
    assert (five.start, five.end, five.wraps) == (990, 1077, True)


def test_windows_are_clipped_on_linear_contigs():
    contig = "".join(random.Random(4).choice("ACGT") for _ in range(200))
    three, five = pc.make_windows(5, 27, "+", 200, False, (15, 50))
    assert (three.start, three.end, three.clipped) == (0, 42, True)
    assert (five.start, five.end, five.clipped) == (0, 77, True)
    assert (three.cand_from, three.cand_to) == (5, 27)
    three, five = pc.make_windows(170, 192, "-", 200, False, (15, 50))
    assert (three.start, three.end, three.clipped) == (155, 200, True)   # long flank cut to 8
    assert (five.start, five.end, five.clipped) == (120, 200, True)
    assert (three.cand_from, five.cand_from) == (8, 8)
    for w in (three, five):
        seq = pc.window_sequence(contig, w, "-")
        assert seq[w.cand_from:w.cand_to] == pc.reverse_complement(contig[170:192])


def test_window_never_longer_than_a_small_circular_contig():
    contig = "".join(random.Random(5).choice("ACGT") for _ in range(60))
    for strand in "+-":
        for w in pc.make_windows(50, 72, strand, 60, True, (15, 50)):
            assert w.length == 60 and w.clipped
            seq = pc.window_sequence(contig, w, strand)
            cand = (contig * 2)[50:72]
            assert seq[w.cand_from:w.cand_to] == (cand if strand == "+" else pc.reverse_complement(cand))


def test_bad_candidates_are_rejected():
    with pytest.raises(ValueError):
        pc.make_windows(190, 210, "+", 200, False)   # beyond a linear contig
    with pytest.raises(ValueError):
        pc.make_windows(10, 30, "*", 200, False)


# ------------------------------------------------------------------ RNAfold

def test_rnafold_command():
    assert pc.rnafold_command(25) == ["RNAfold", "--noPS", "-T", "25"]
    assert pc.rnafold_command(37.5, "/opt/RNAfold") == ["/opt/RNAfold", "--noPS", "-T", "37.5"]


def test_parse_stored_rnafold_output():
    records = pc.parse_rnafold((DATA / "rnafold_T25.out").read_text())
    assert [r[0] for r in records] == ["w1", "w2", "let7"]
    assert records[0][2] == "(((((((((....))))))))).....((((((....)))))).."
    assert records[0][3] == -36.53
    assert records[1][2] == "." * 20 and records[1][3] == 0.0       # "(  0.00)"
    assert records[2][1] == LET7.replace("T", "U")
    assert records[2][3] == -36.32
    assert all(len(seq) == len(structure) for _, seq, structure, _ in records)


def test_parse_rnafold_without_headers_and_errors():
    assert pc.parse_rnafold("ACGUNNACGU\n(((....))) ( -0.41)\n") == [("", "ACGUNNACGU", "(((....)))", -0.41)]
    with pytest.raises(ValueError):
        pc.parse_rnafold(">a\nACGU\n")
    with pytest.raises(ValueError):
        pc.parse_rnafold(">a\nACGU\n.. (0.00)\n")


def test_missing_rnafold_gives_a_clear_error():
    with pytest.raises(RuntimeError, match="not found on PATH"):
        pc.run_rnafold(["ACGT"], 25, executable="RNAfold-not-installed")


# ------------------------------------------------------------------ measures

def test_amfe_and_mfei():
    seq = "GGCC" + "AATT" * 4                 # 20 nt, 4 G/C = 20%
    assert pc.gc_percent(seq) == pytest.approx(20.0)
    assert pc.amfe(-10.0, 20) == pytest.approx(50.0)       # -(-10) / 20 * 100
    assert pc.mfei(-10.0, seq) == pytest.approx(2.5)       # 50 / 20
    assert pc.mfei(0.0, seq) == 0.0
    assert pc.mfei(-3.0, "ATATATAT") is None               # no G or C
    assert pc.gc_percent("ggcc") == 100.0


def test_shuffle_pvalue():
    shuffled = [-5.0, -10.0, -12.0, -20.0]
    assert pc.shuffle_pvalue(-10.0, shuffled) == pytest.approx(4 / 5)   # ties count: k = 3
    assert pc.shuffle_pvalue(-30.0, shuffled) == pytest.approx(1 / 5)   # never 0
    assert pc.shuffle_pvalue(0.0, shuffled) == 1.0
    assert pc.shuffle_pvalue(-10.0, []) is None
    assert pc.shuffle_pvalue(-50.0, [-1.0] * 999) == pytest.approx(0.001)


def test_shuffles_are_reproducible_and_keep_dinucleotides():
    a = pc.shuffled_windows(LET7, 5, 1, "c1", "3p")
    assert a == pc.shuffled_windows(LET7, 5, 1, "c1", "3p")
    assert a != pc.shuffled_windows(LET7, 5, 1, "c1", "5p")
    assert a != pc.shuffled_windows(LET7, 5, 2, "c1", "3p")
    pairs = sorted(x + y for x, y in zip(LET7, LET7[1:]))
    assert all(sorted(x + y for x, y in zip(s, s[1:])) == pairs for s in a)


def test_shuffles_do_not_depend_on_the_hash_seed():
    # set iteration order changes between Python processes (PYTHONHASHSEED)
    code = ("from smithhunter import precursor as pc; "
            f"print(pc.shuffled_windows({LET7!r}, 3, 1, 'c1', '3p'))")
    src = str(Path(pc.__file__).parents[1])
    outs = {subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True,
                           env={**os.environ, "PYTHONPATH": src, "PYTHONHASHSEED": seed}).stdout
            for seed in ("1", "2", "3", "4")}
    assert len(outs) == 1


# ------------------------------------------------------------------ hairpin check

def test_perfect_hairpin_candidate_in_each_arm():
    structure = "((((((((((....))))))))))"
    five = pc.hairpin_check(structure, 0, 8)
    assert five["hairpin_loops"] == 1 and five["stem_pairs"] == 10
    assert five["stem_fraction"] == 1.0 and five["single_hairpin"] == 1
    assert (five["in_arm"], five["arm"], five["paired_fraction"]) == (1, "5p", 1.0)
    three = pc.hairpin_check(structure, 15, 24)
    assert (three["in_arm"], three["arm"], three["paired_fraction"]) == (1, "3p", 1.0)


def test_hairpin_with_bulge_and_unpaired_ends():
    structure = "..((((.((((....))))..))))..."
    hp = pc.hairpin_check(structure, 0, 8)      # two unpaired bases, a bulge
    assert hp["hairpin_loops"] == 1 and hp["stem_pairs"] == 8 and hp["single_hairpin"] == 1
    assert hp["arm"] == "5p" and hp["paired_fraction"] == pytest.approx(5 / 8)


def test_candidate_spanning_the_loop_is_in_no_arm():
    hp = pc.hairpin_check("((((((((((....))))))))))", 6, 18)
    assert hp["single_hairpin"] == 1
    assert (hp["in_arm"], hp["arm"]) == (0, "both")
    assert hp["paired_fraction"] == pytest.approx(8 / 12)


def test_multiloop_is_not_a_single_hairpin():
    structure = "((((..((((...))))..((((...))))..))))"
    hp = pc.hairpin_check(structure, 0, 4)       # candidate in the closing stem
    assert hp["hairpin_loops"] == 2
    assert hp["stem_pairs"] == 4 and hp["stem_fraction"] == pytest.approx(4 / 12)
    assert hp["single_hairpin"] == 0
    assert (hp["in_arm"], hp["arm"], hp["paired_fraction"]) == (0, "NA", 1.0)
    assert pc.hairpin_check(structure, 6, 10)["arm"] == "5p"    # arm of the first branch


def test_two_separate_hairpins():
    structure = "((((((....))))))..((((....)))).."
    hp = pc.hairpin_check(structure, 10, 16)
    assert hp["hairpin_loops"] == 2 and hp["stem_pairs"] == 6
    assert hp["stem_fraction"] == pytest.approx(0.6) and hp["single_hairpin"] == 0
    assert (hp["in_arm"], hp["arm"]) == (1, "3p")
    other = pc.hairpin_check(structure, 18, 22)   # paired, but not in the dominant stem-loop
    assert (other["in_arm"], other["arm"], other["paired_fraction"]) == (0, "NA", 1.0)


def test_small_side_hairpin_leaves_a_dominant_stem_loop():
    structure = "((((((((((((((((....))))))))))))))))..(((...)))"
    hp = pc.hairpin_check(structure, 0, 10)
    assert hp["hairpin_loops"] == 2 and hp["stem_fraction"] == pytest.approx(16 / 19)
    assert hp["single_hairpin"] == 1


def test_unfolded_window_and_bad_structures():
    hp = pc.hairpin_check("." * 30, 5, 25)
    assert hp == {"hairpin_loops": 0, "stem_pairs": 0, "stem_fraction": 0.0, "single_hairpin": 0,
                  "in_arm": 0, "arm": "NA", "paired_fraction": 0.0}
    with pytest.raises(ValueError):
        pc.hairpin_check("((..)", 0, 2)
    with pytest.raises(ValueError):
        pc.hairpin_check("(..))", 0, 2)


# ------------------------------------------------------------------ pipeline

def fake_fold(sequences, temperature, executable):
    """Stand-in for RNAfold: MFE = -(number of G), a hairpin when the sequence allows."""
    out = []
    for s in sequences:
        n = min(10, (len(s) - 4) // 2)
        out.append(("(" * n + "." * (len(s) - 2 * n) + ")" * n, -float(s.count("G"))))
    return out


def write_run(tmp_path, contig):
    (tmp_path / "g.fa").write_text(f">chr1 test\n{contig[:100]}\n{contig[100:]}\n")
    (tmp_path / "genomes.tsv").write_text(
        f"name\tfasta\trole\tcircular\nmito\t{tmp_path / 'g.fa'}\tfocus\t1\n")
    (tmp_path / "contigs.tsv").write_text(
        f"ref\tgenome\tcontig\tlength\tcircular\nmito|chr1\tmito\tchr1\t{len(contig)}\t1\n")
    head = "locus_id\tgenome\tcontig\tstart\tend\tstrand\trep_start\trep_end\tpass\n"
    (tmp_path / "loci.tsv").write_text(
        head + "L1\tmito\tchr1\t101\t125\t+\t101\t122\t1\n"
        "L2\tmito\tchr1\t295\t20\t-\t296\t17\t1\n"
        "L3\tmito\tchr1\t200\t230\t+\t200\t221\t0\n")


def test_select_loci(tmp_path):
    write_run(tmp_path, "A" * 300)
    loci = str(tmp_path / "loci.tsv")
    passing, unknown = pc.select_loci(loci)
    assert [c["candidate"] for c in passing] == ["L1", "L2"] and unknown == []
    assert (passing[0]["start"], passing[0]["end"]) == (101, 122)
    chosen, unknown = pc.select_loci(loci, ["L3", "X9"], anchor="locus")
    assert [(c["candidate"], c["start"], c["end"]) for c in chosen] == [("L3", 200, 230)]
    assert unknown == ["X9"]


def test_fold_candidates_table_and_thread_independence(tmp_path):
    contig = "".join(random.Random(7).choice("ACGT") for _ in range(300))
    write_run(tmp_path, contig)
    candidates, _ = pc.select_loci(str(tmp_path / "loci.tsv"))
    contigs = pc.load_contigs(str(tmp_path / "genomes.tsv"), str(tmp_path / "contigs.tsv"),
                              {("mito", "chr1")})
    assert contigs[("mito", "chr1")] == (contig, True)
    rows = pc.fold_candidates(candidates, contigs, 25, shuffles=20, threads=1, fold=fake_fold)
    again = pc.fold_candidates(candidates, contigs, 25, shuffles=20, threads=4, fold=fake_fold)
    assert rows == again
    assert [(r["candidate"], r["window"]) for r in rows] == [
        ("L1", "3p"), ("L1", "5p"), ("L2", "3p"), ("L2", "5p")]
    first = rows[0]
    assert (first["window_start"], first["window_end"], first["length"]) == (51, 137, 87)
    assert first["sequence"] == contig[50:137]
    assert first["sequence"][first["candidate_from"] - 1:first["candidate_to"]] == contig[100:122]
    wrap = rows[2]                    # minus strand, candidate 296..17 across the origin
    assert wrap["wraps_origin"] == 1
    assert (wrap["window_start"], wrap["window_end"]) == (281, 67)
    cand = pc.reverse_complement(contig[295:] + contig[:17])
    assert wrap["sequence"][wrap["candidate_from"] - 1:wrap["candidate_to"]] == cand
    for r in rows:
        assert float(r["amfe"]) == pytest.approx(-float(r["mfe"]) / r["length"] * 100, abs=0.01)
        assert r["shuffles"] == 20 and 1 / 21 <= float(r["pvalue"]) <= 1
        # the fake MFE depends on composition only, which shuffling keeps: every shuffle ties
        assert float(r["pvalue"]) == 1.0
    for i in (0, 2):
        assert rows[i]["best"] + rows[i + 1]["best"] == 1
    # equal p-values: the lower MFE (more G) wins
    for a, b in ((rows[0], rows[1]), (rows[2], rows[3])):
        if float(a["mfe"]) != float(b["mfe"]):
            assert (a if float(a["mfe"]) < float(b["mfe"]) else b)["best"] == 1
    # the result of one candidate does not depend on the others
    alone = pc.fold_candidates(candidates[1:], contigs, 25, shuffles=20, fold=fake_fold)
    assert alone == rows[2:]


def test_no_shuffles_gives_no_pvalue(tmp_path):
    contig = "".join(random.Random(8).choice("ACGT") for _ in range(300))
    write_run(tmp_path, contig)
    candidates, _ = pc.select_loci(str(tmp_path / "loci.tsv"))
    contigs = {("mito", "chr1"): (contig, True)}
    rows = pc.fold_candidates(candidates, contigs, 25, shuffles=0, fold=fake_fold)
    assert all(r["pvalue"] == "NA" and r["shuffles"] == 0 for r in rows)


def test_cli_requires_the_temperature(tmp_path):
    write_run(tmp_path, "ACGT" * 75)
    args = ["precursors", "--loci", str(tmp_path / "loci.tsv"), "--genomes",
            str(tmp_path / "genomes.tsv"), "--contigs", str(tmp_path / "contigs.tsv"),
            "--out", str(tmp_path / "p.tsv")]
    with pytest.raises(SystemExit):
        main(args)
    with pytest.raises(SystemExit, match="not found on PATH"):
        main(args + ["--temperature", "25", "--rnafold", "RNAfold-not-installed"])


# ------------------------------------------------------------------ real RNAfold

@pytest.mark.skipif(shutil.which("RNAfold") is None, reason="RNAfold (ViennaRNA) not installed")
def test_real_rnafold_end_to_end(tmp_path):
    stored = pc.parse_rnafold((DATA / "rnafold_T25.out").read_text())
    assert pc.run_rnafold([LET7], 25) == [(stored[2][2], stored[2][3])]
    assert pc.run_rnafold([LET7], 37)[0][1] > stored[2][3]      # less stable when warmer

    # a let-7 hairpin planted in a random circular contig, candidate = its 5' arm
    rng = random.Random(11)
    flank = "".join(rng.choice("ACGT") for _ in range(150))
    contig = flank[:100] + LET7 + flank[100:]
    write_run(tmp_path, contig)
    (tmp_path / "loci.tsv").write_text(
        "locus_id\tgenome\tcontig\tstart\tend\tstrand\trep_start\trep_end\tpass\n"
        "let7\tmito\tchr1\t101\t122\t+\t101\t122\t1\n")
    out = tmp_path / "precursors.tsv"
    args = ["precursors", "--loci", str(tmp_path / "loci.tsv"), "--genomes",
            str(tmp_path / "genomes.tsv"), "--contigs", str(tmp_path / "contigs.tsv"),
            "--temperature", "25", "--shuffles", "99", "--out", str(out)]
    main(args + ["--threads", "1"])
    one = out.read_text()
    main(args + ["--threads", "3"])
    assert out.read_text() == one
    rows = read_tsv(str(out))
    assert list(rows[0]) == pc.COLUMNS and len(rows) == 2
    best = [r for r in rows if r["best"] == "1"]
    assert len(best) == 1 and best[0]["window"] == "5p"
    assert float(best[0]["pvalue"]) <= 0.05
    assert best[0]["single_hairpin"] == "1" and best[0]["arm"] == "5p"
    assert float(best[0]["paired_fraction"]) > 0.8
    assert float(best[0]["mfei"]) > 0.85
