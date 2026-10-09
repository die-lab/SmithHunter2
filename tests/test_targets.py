import random
from collections import Counter

import pytest

from smithhunter.targets.fdr import qvalues, site_count_fdr
from smithhunter.targets.pipeline import read_candidates, run_sites
from smithhunter.targets.regions import (Region, read_regions, regions_from_genome,
                                         regions_from_transcripts, revcomp, write_regions)
from smithhunter.targets.shuffle import dinucleotide_shuffle, make_decoys
from smithhunter.targets.sites import scan, supplementary_pairing

GUIDE = "TAGCTTATCAGACTGATGTTGA"  # miR-21: seed 2-8 AGCTTAT
M8 = revcomp(GUIDE[1:8])           # ATAAGCT, pairs with guide 2-8
M27 = revcomp(GUIDE[1:7])          # TAAGCT, guide 2-7
M38 = revcomp(GUIDE[2:8])          # ATAAGC, guide 3-8
PAD = "CCCCCCCCCC"


def classes(target, **kw):
    return [(s.site_class, s.start, s.end) for s in scan(GUIDE, target, **kw)]


def test_canonical_site_classes():
    assert classes(PAD + M8 + "A" + PAD) == [("8mer", 10, 18)]
    assert classes(PAD + M8 + "C" + PAD) == [("7mer-m8", 10, 17)]
    # guide 2-7 matched, position 8 not: G instead of A before TAAGCT
    assert classes(PAD + "G" + M27 + "A" + PAD) == [("7mer-A1", 11, 18)]
    assert classes(PAD + "G" + M27 + "C" + PAD) == [("6mer", 11, 17)]
    # guide 3-8 matched, guide 2 not (T instead of the A pairing with guide 2)
    assert classes(PAD + M38 + "C" + "C" + PAD) == [("offset-6mer", 10, 16)]
    assert classes(PAD * 3) == []


def test_each_site_counted_once_and_classes_filter():
    target = PAD + M8 + "A" + PAD + "G" + M27 + "C" + PAD
    assert [c for c, _, _ in classes(target)] == ["8mer", "6mer"]
    assert [c for c, _, _ in classes(target, classes=("6mer",))] == ["6mer"]


def test_other_seed_range():
    pattern = revcomp(GUIDE[3:10])
    assert classes(PAD + pattern + PAD, seed=(4, 10)) == [("seed4-10", 10, 17)]


def test_supplementary_pairing():
    # guide 13-16 (CTGA) pairs with the target 12-15 nt upstream of t1
    t1 = 30
    target = list("C" * 40)
    for k, base in zip(range(13, 17), GUIDE[12:16]):
        target[t1 - (k - 1)] = {"A": "T", "C": "G", "G": "C", "T": "A"}[base]
    assert supplementary_pairing(GUIDE, "".join(target), t1) == 4
    assert supplementary_pairing(GUIDE, "C" * 40, t1) == 1  # only guide 15 (G) pairs with C


def write(path, text):
    path.write_text(text)
    return str(path)


def test_regions_from_genome_minus_strand_spliced(tmp_path):
    rng = random.Random(3)
    genome = "".join(rng.choice("ACGT") for _ in range(200))
    fasta = write(tmp_path / "g.fa", f">chr1\n{genome}\n")
    # minus-strand mRNA with exons 21-60 and 101-150 (1-based), CDS 41-60 + 101-120
    gff = write(tmp_path / "g.gff3", "\n".join([
        "chr1\tt\tgene\t21\t150\t.\t-\t.\tID=g1;Name=GENE1",
        "chr1\tt\tmRNA\t21\t150\t.\t-\t.\tID=t1;Parent=g1",
        "chr1\tt\texon\t21\t60\t.\t-\t.\tParent=t1",
        "chr1\tt\texon\t101\t150\t.\t-\t.\tParent=t1",
        "chr1\tt\tCDS\t41\t60\t.\t-\t0\tParent=t1",
        "chr1\tt\tCDS\t101\t120\t.\t-\t0\tParent=t1",
    ]) + "\n")
    regions = {r.region: r for r in regions_from_genome(
        fasta, gff, ["five_prime_UTR", "CDS", "three_prime_UTR"])}
    transcript = revcomp(genome[20:60] + genome[100:150])
    assert regions["five_prime_UTR"].sequence == revcomp(genome[120:150])
    assert regions["CDS"].sequence == revcomp(genome[40:60] + genome[100:120])
    assert regions["three_prime_UTR"].sequence == revcomp(genome[20:40])
    assert "".join(regions[k].sequence for k in ("five_prime_UTR", "CDS", "three_prime_UTR")) \
        == transcript
    utr3 = regions["three_prime_UTR"]
    assert utr3.gene == "GENE1" and utr3.blocks == [(20, 40)]
    assert utr3.genomic(0) == 39 and utr3.genomic(19) == 20  # minus strand runs backwards
    cds = regions["CDS"]
    assert cds.genomic(0) == 119 and cds.genomic(20) == 59

    path = str(tmp_path / "regions.tsv")
    write_regions(path, list(regions.values()))
    back = {r.region: r for r in read_regions(path)}
    assert back["CDS"].blocks == cds.blocks and back["CDS"].sequence == cds.sequence


def test_regions_without_exon_features(tmp_path):
    genome = "A" * 10 + "C" * 30 + "G" * 10
    fasta = write(tmp_path / "g.fa", f">c\n{genome}\n")
    gff = write(tmp_path / "g.gff3",
                "c\tt\tmRNA\t1\t50\t.\t+\t.\tID=m\n"
                "c\tt\tfive_prime_UTR\t1\t10\t.\t+\t.\tParent=m\n"
                "c\tt\tCDS\t11\t40\t.\t+\t0\tParent=m\n"
                "c\tt\tthree_prime_UTR\t41\t50\t.\t+\t.\tParent=m\n")
    [utr] = regions_from_genome(fasta, gff, ["three_prime_UTR"])
    assert utr.sequence == "G" * 10
    # no CDS at all: whole transcript only when asked
    gff2 = write(tmp_path / "n.gff3", "c\tt\tmRNA\t1\t50\t.\t+\t.\tID=m\n"
                                      "c\tt\texon\t1\t50\t.\t+\t.\tParent=m\n")
    assert regions_from_genome(fasta, gff2, ["three_prime_UTR"]) == []
    [whole] = regions_from_genome(fasta, gff2, ["three_prime_UTR"], whole_if_missing=True)
    assert whole.region == "transcript" and whole.sequence == genome


def test_regions_from_transcripts(tmp_path):
    fasta = write(tmp_path / "t.fa", ">tA\nACGUACGUAC\n>tB\nGGGGG\n")
    assert [r.sequence for r in regions_from_transcripts(fasta)] == ["ACGTACGTAC", "GGGGG"]
    table = write(tmp_path / "r.tsv", "transcript\tregion\tstart\tend\n"
                                      "tA\tthree_prime_UTR\t3\t6\n")
    [r] = regions_from_transcripts(fasta, table, ["three_prime_UTR"])
    assert r.sequence == "GTAC" and r.id == "tA|three_prime_UTR"


def dinucs(seq):
    return Counter(seq[i:i + 2] for i in range(len(seq) - 1))


def test_dinucleotide_shuffle_preserves_counts():
    rng = random.Random(7)
    for _ in range(50):
        seq = "".join(rng.choice("ACGT") for _ in range(rng.randint(3, 30)))
        shuffled = dinucleotide_shuffle(seq, rng)
        assert dinucs(shuffled) == dinucs(seq)
        assert shuffled[0] == seq[0] and shuffled[-1] == seq[-1]
    assert len({dinucleotide_shuffle(GUIDE, rng) for _ in range(30)}) > 5


def test_decoys_never_share_a_real_seed():
    cands = [("a", GUIDE), ("b", "TGAGGTAGTAGGTTGTATAGTT")]
    decoys = make_decoys(cands, 25, random_seed=2)
    assert len(decoys) == 50
    seeds = {s[1:8] for _, s in cands}
    assert all(d[2][1:8] not in seeds for d in decoys)
    assert {d[1] for d in decoys} == {"a", "b"}


def test_qvalues():
    # 4 real scores, 2 decoys per candidate; lower is better
    q = qvalues([-30, -25, -20, -10], [-22, -12, -11, -5], n_decoys=2)
    # at -30, -25: no decoys pass -> 0; at -20: 1 decoy / 2 / 3 real = 0.1667
    # at -10: 3 decoys / 2 / 4 = 0.375
    assert q == [0.0, 0.0, pytest.approx(0.166667), 0.375]
    # higher is better: 5 beats the decoy (q 0), 1 does not (1 decoy / 2 real)
    assert qvalues([5, 1], [4], n_decoys=1, higher_is_better=True) == [0.0, 0.5]
    assert qvalues([], [1], 1) == []


def test_site_count_fdr():
    rows = {r["site_class"]: r for r in site_count_fdr(
        {"8mer": 4, "6mer": 10}, [{"8mer": 1, "6mer": 10}, {"8mer": 1, "6mer": 12}],
        ["8mer", "6mer"])}
    assert rows["8mer"]["fdr"] == 0.25 and rows["8mer"]["enrichment"] == 4.0
    assert rows["6mer"]["fdr"] == 1.0
    assert rows["any"]["real_sites"] == 14


def test_run_sites_end_to_end(tmp_path):
    utr = PAD + M8 + "A" + PAD
    regions = [Region("t1|three_prime_UTR", "t1", "G1", "three_prime_UTR", utr, "c", "+",
                      [(100, 100 + len(utr))]),
               # an isoform with the same UTR: one site in the counts, two rows in the table
               Region("t2|three_prime_UTR", "t2", "G1", "three_prime_UTR", utr, "c", "+",
                      [(100, 100 + len(utr))])]
    decoys = make_decoys([("smith1", GUIDE)], 5)
    classes, sites, targets, fdr = run_sites([("smith1", GUIDE)], decoys, regions)
    assert [(s["transcript"], s["site_class"], s["genomic"]) for s in sites] == \
        [("t1", "8mer", "c:111-118(+)"), ("t2", "8mer", "c:111-118(+)")]
    real = {r["site_class"]: r["real_sites"] for r in fdr if r["candidate"] == "smith1"}
    assert real["8mer"] == 1 and real["any"] == 1
    assert {t["transcript"] for t in targets} == {"t1", "t2"}
    assert targets[0]["best_class"] == "8mer" and targets[0]["n_8mer"] == 1


def test_read_candidates(tmp_path):
    path = write(tmp_path / "c.fa",
                 ">mito_00001 mito|chrM:10-31(+) locus=10-40 count=50\nACGUACGUACGUACGUACGU\n"
                 ">custom\nTTTTTTTTTTTTTTTTTTTT\n")
    assert read_candidates(path) == [("mito_00001", "mito", "ACGTACGTACGTACGTACGT"),
                                     ("custom", "", "T" * 20)]
