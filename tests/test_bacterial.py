"""Bacterial mode: start-codon windows (regions) and IntaRNA scoring with decoy q-values."""

import os
import random
import shutil

import pytest

from smithhunter.cli import main
from smithhunter.targets import bacterial, intarna
from smithhunter.targets.regions import (Region, bacterial_windows, parse_cds, read_regions,
                                         revcomp, write_regions)
from smithhunter.targets.shuffle import make_decoys

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def rnd(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def gff_line(contig, ftype, start, end, strand, attrs, phase="."):
    return "\t".join([contig, "test", ftype, str(start), str(end), ".", strand, phase, attrs])


def write_genome(tmp_path, contigs, gff_lines, name="g"):
    fa, gff = tmp_path / f"{name}.fa", tmp_path / f"{name}.gff3"
    fa.write_text("".join(f">{n}\n{s}\n" for n, s in contigs.items()))
    gff.write_text("##gff-version 3\n" + "\n".join(gff_lines) + "\n")
    return str(fa), str(gff)


def cds(contig, start, end, strand, name, phase="0", extra=""):
    """A CDS line (1-based, inclusive) in NCBI style: gene + CDS with locus_tag."""
    return [gff_line(contig, "gene", start, end, strand, f"ID=gene-{name};Name={name};locus_tag=L_{name}"),
            gff_line(contig, "CDS", start, end, strand,
                     f"ID=cds-{name};Parent=gene-{name};locus_tag=L_{name}{extra}", phase)]


def check_mapping(r: Region, genome: str, circular: bool):
    """Every window position maps to the genomic base it was cut from (complement on '-')."""
    assert sum(e - s for s, e in r.blocks) == len(r.sequence)
    for i, base in enumerate(r.sequence):
        g = r.genomic(i)
        assert g is not None
        expected = genome[g] if r.strand == "+" else COMP[genome[g]]
        assert base == expected, (r.id, i)


# ---- windows ------------------------------------------------------------------------------

def genome_with_start_codons(length, plus_starts, minus_starts, seed=1):
    """Random contig with ATG at the 0-based start positions given (minus: the position of the
    A, i.e. CAT at start-2..start on the plus strand)."""
    g = list(rnd(length, seed))
    for p in plus_starts:
        g[p:p + 3] = "ATG"
    for p in minus_starts:
        g[p - 2:p + 1] = "CAT"
    return "".join(g)


def test_window_on_both_strands(tmp_path):
    genome = genome_with_start_codons(3000, [1000], [2299])
    fa, gff = write_genome(tmp_path, {"c1": genome},
                           cds("c1", 1001, 1300, "+", "plus") + cds("c1", 2001, 2300, "-", "minus"))
    regions = {r.gene: r for r in bacterial_windows(fa, gff, (-200, 100))}
    p, m = regions["plus"], regions["minus"]
    assert p.region == "start_window" and p.id == "L_plus|start_window" and p.transcript == "L_plus"
    assert p.sequence == genome[800:1100] and p.anchor == 200 and p.blocks == [(800, 1100)]
    assert p.sequence[200:203] == "ATG" and p.cds_nt == 100
    assert m.sequence == revcomp(genome[2200:2500]) and m.anchor == 200
    assert m.sequence[200:203] == "ATG" and m.blocks == [(2200, 2500)]
    for r in (p, m):
        assert len(r.sequence) == 300
        check_mapping(r, genome, False)
    assert p.genomic(200) == 1000 and m.genomic(200) == 2299      # the A of the start codon


def test_window_size_is_configurable(tmp_path):
    genome = genome_with_start_codons(2000, [1000], [])
    fa, gff = write_genome(tmp_path, {"c1": genome}, cds("c1", 1001, 1400, "+", "a"))
    (r,) = bacterial_windows(fa, gff, (-50, 20))
    assert r.sequence == genome[950:1020] and r.anchor == 50
    with pytest.raises(ValueError):
        bacterial_windows(fa, gff, (0, 100))
    with pytest.raises(ValueError):
        bacterial_windows(fa, gff, (100, -200))


def test_linear_contig_ends_clip_the_window(tmp_path):
    genome = genome_with_start_codons(500, [49], [449], seed=2)
    # + gene starting 49 nt from the contig start; - gene whose upstream (higher coordinates)
    # leaves only 50 nt, and whose CDS ends 120 nt into the contig
    fa, gff = write_genome(tmp_path, {"c1": genome},
                           cds("c1", 50, 300, "+", "left") + cds("c1", 121, 450, "-", "right"))
    stats = {}
    regions = {r.gene: r for r in bacterial_windows(fa, gff, (-200, 100), stats=stats)}
    left, right = regions["left"], regions["right"]
    assert left.sequence == genome[0:149] and left.anchor == 49          # 49 nt upstream only
    assert right.sequence == revcomp(genome[350:500]) and right.anchor == 50   # 50 nt upstream
    assert len(right.sequence) == 150 and right.sequence[50:53] == "ATG"
    for r in regions.values():
        check_mapping(r, genome, False)
    assert stats["clipped"] == 2 and stats["wrapped"] == 0 and stats["windows"] == 2
    # a window that is too short after cutting is dropped
    stats = {}
    assert bacterial_windows(fa, gff, (-200, 100), min_length=200, stats=stats) == []
    assert stats["too_short"] == 2


@pytest.mark.parametrize("length", [1000, 150])
def test_circular_contig_wraps_both_strands(tmp_path, length):
    # the + gene is near the origin (its upstream wraps), the - gene near the end (same)
    plus_start, minus_start = 50, length - 50
    genome = genome_with_start_codons(length, [plus_start], [minus_start], seed=3)
    fa, gff = write_genome(
        tmp_path, {"c1": genome},
        cds("c1", plus_start + 1, plus_start + 90, "+", "p") +
        cds("c1", minus_start - 89, minus_start + 1, "-", "m"))
    stats = {}
    regions = {r.gene: r for r in bacterial_windows(fa, gff, (-200, 100), circular=True,
                                                    stats=stats)}
    assert stats["wrapped"] == 2 and stats["clipped"] == 0
    for r in regions.values():
        assert len(r.sequence) == 300 and r.anchor == 200
        assert r.sequence[200:203] == "ATG"
        check_mapping(r, genome, True)
    if length == 1000:
        p, m = regions["p"], regions["m"]
        assert p.blocks == [(850, 1000), (0, 150)]
        assert p.sequence == genome[850:] + genome[:150]
        # minus strand: transcript order starts at the wrapped (highest) part
        assert m.blocks == [(0, 151), (851, 1000)]
        assert m.sequence == revcomp(genome[851:] + genome[:151])
    # the same genes on a linear contig are cut
    stats = {}
    for r in bacterial_windows(fa, gff, (-200, 100), circular=False, stats=stats):
        assert len(r.sequence) < 300
    assert stats["clipped"] == 2 and stats["wrapped"] == 0


def test_circular_contigs_can_be_chosen_by_name(tmp_path):
    g1 = genome_with_start_codons(600, [20], [], seed=4)
    g2 = genome_with_start_codons(600, [20], [], seed=5)
    fa, gff = write_genome(tmp_path, {"chr": g1, "plasmid": g2},
                           cds("chr", 21, 200, "+", "a") + cds("plasmid", 21, 200, "+", "b"))
    regions = {r.gene: r for r in bacterial_windows(fa, gff, (-100, 50), circular={"plasmid"})}
    assert len(regions["a"].sequence) == 20 + 50 and regions["a"].anchor == 20
    assert len(regions["b"].sequence) == 150 and regions["b"].anchor == 100
    assert len(regions["b"].blocks) == 2


def test_upstream_cds_counts_the_preceding_gene_of_the_operon(tmp_path):
    genome = genome_with_start_codons(2000, [700, 1400], [1000], seed=6)
    # A: 301-600; B starts at 701, 100 nt after A ends; C (minus strand) overlaps B's upstream
    # and must not count
    fa, gff = write_genome(tmp_path, {"c1": genome},
                           cds("c1", 301, 600, "+", "A") + cds("c1", 701, 1000, "+", "B") +
                           cds("c1", 501, 1001, "-", "C"))
    regions = {r.gene: r for r in bacterial_windows(fa, gff, (-200, 100))}
    assert regions["B"].upstream_cds == 100          # B window upstream = 500..700, A ends at 600
    assert regions["A"].upstream_cds == 0
    assert regions["C"].upstream_cds == 0


def test_gff_flavours_and_skipped_cds(tmp_path):
    genome = genome_with_start_codons(3000, [100, 400, 800, 1200, 1600, 2000], [], seed=7)
    lines = (
        cds("c1", 101, 300, "+", "ncbi") +
        # mRNA / CDS style (Parent = mRNA, mRNA Parent = gene)
        [gff_line("c1", "gene", 401, 600, "+", "ID=g2;Name=mrna_gene"),
         gff_line("c1", "mRNA", 401, 600, "+", "ID=t2;Parent=g2"),
         gff_line("c1", "CDS", 401, 600, "+", "ID=cds2;Parent=t2", "0")] +
        # Prokka style: no gene line, no Parent; gene symbol only on the CDS
        [gff_line("c1", "CDS", 801, 1000, "+", "ID=PROKKA_001;gene=symA;locus_tag=PROKKA_001", "0")] +
        # 5' partial (start_range) and phase 1 and pseudogene: skipped
        [gff_line("c1", "CDS", 1201, 1400, "+", "ID=part1;locus_tag=P1;start_range=.,1201;partial=true", "0"),
         gff_line("c1", "CDS", 1601, 1800, "+", "ID=part2;locus_tag=P2", "1"),
         gff_line("c1", "CDS", 2001, 2200, "+", "ID=ps;locus_tag=P3;pseudo=true", "0"),
         gff_line("c1", "CDS", 2401, 2600, "?", "ID=nostrand;locus_tag=P4", "0")])
    fa, gff = write_genome(tmp_path, {"c1": genome}, lines)
    stats = {}
    regions = bacterial_windows(fa, gff, (-50, 50), stats=stats)
    assert {(r.transcript, r.gene) for r in regions} == {
        ("L_ncbi", "ncbi"), ("cds2", "mrna_gene"), ("PROKKA_001", "symA")}
    assert stats["partial"] == 2 and stats["pseudo"] == 1 and stats["windows"] == 3


def test_cds_in_two_segments_is_one_window_at_the_first_start_codon(tmp_path):
    genome = genome_with_start_codons(1500, [300], [1100], seed=8)
    lines = [gff_line("c1", "CDS", 301, 500, "+", "ID=fs;locus_tag=FS", "0"),
             gff_line("c1", "CDS", 502, 700, "+", "ID=fs;locus_tag=FS", "0"),     # frameshift
             gff_line("c1", "CDS", 901, 1000, "-", "ID=fm;locus_tag=FM", "0"),
             gff_line("c1", "CDS", 1002, 1101, "-", "ID=fm;locus_tag=FM", "0")]
    fa, gff = write_genome(tmp_path, {"c1": genome}, lines)
    (fs,) = [c for c in parse_cds(gff) if c.id == "FS"]
    (fm,) = [c for c in parse_cds(gff) if c.id == "FM"]
    assert fs.start == 300 and fm.start == 1100 and fs.length == 200 + 199
    regions = {r.transcript: r for r in bacterial_windows(fa, gff, (-20, 30))}
    assert len(regions) == 2 and regions["FS"].sequence[20:23] == "ATG"
    assert regions["FM"].sequence[20:23] == "ATG"


def test_duplicate_locus_tags_get_distinct_ids(tmp_path):
    genome = genome_with_start_codons(2000, [100, 600], [], seed=9)
    lines = [gff_line("c1", "CDS", 101, 300, "+", "ID=a;locus_tag=X", "0"),
             gff_line("c1", "CDS", 601, 800, "+", "ID=b;locus_tag=X", "0")]
    fa, gff = write_genome(tmp_path, {"c1": genome}, lines)
    ids = [r.id for r in bacterial_windows(fa, gff, (-50, 50))]
    assert len(set(ids)) == 2


def test_regions_file_roundtrip_keeps_the_anchor_and_leaves_ago_tables_alone(tmp_path):
    genome = genome_with_start_codons(1000, [50], [], seed=10)
    fa, gff = write_genome(tmp_path, {"c1": genome}, cds("c1", 51, 300, "+", "a"))
    regions = bacterial_windows(fa, gff, (-100, 50), circular=True)
    path = tmp_path / "regions.tsv"
    write_regions(str(path), regions)
    header = path.read_text().splitlines()[0].split("\t")
    assert header[:9] == ["region_id", "transcript", "gene", "region", "length", "contig",
                          "strand", "blocks", "sequence"] and "anchor" in header
    back = read_regions(str(path))[0]
    assert back.sequence == regions[0].sequence and back.anchor == regions[0].anchor == 100
    assert back.blocks == regions[0].blocks and back.upstream_cds == 0 and back.cds_nt == 50
    ago = tmp_path / "ago.tsv"
    write_regions(str(ago), [Region("r", "t", "g", "three_prime_UTR", "ACGT", "c", "+", [(0, 4)])])
    assert ago.read_text().splitlines()[0].split("\t")[-1] == "sequence"
    assert read_regions(str(ago))[0].anchor is None


# ---- coordinates --------------------------------------------------------------------------

def test_relative_coordinates_have_no_position_zero():
    # window of 5 upstream nt then the start codon: index 5 is the A (+1), index 4 is -1
    assert [bacterial.to_rel(p, 5) for p in range(1, 9)] == [-5, -4, -3, -2, -1, 1, 2, 3]
    cfg = bacterial.BacterialConfig()
    assert bacterial.derive({"intarna_t_start": 1, "intarna_t_end": 4}, 200, cfg) == {
        "rel_start": -200, "rel_end": -197, "hit_region": "5UTR", "hit_rbs": 0}
    d = bacterial.derive({"intarna_t_start": 180, "intarna_t_end": 190}, 200, cfg)
    assert (d["rel_start"], d["rel_end"], d["hit_region"], d["hit_rbs"]) == (-21, -11, "5UTR", 1)
    d = bacterial.derive({"intarna_t_start": 195, "intarna_t_end": 210}, 200, cfg)
    assert (d["rel_start"], d["rel_end"], d["hit_region"]) == (-6, 10, "spans_start")
    d = bacterial.derive({"intarna_t_start": 201, "intarna_t_end": 230}, 200, cfg)
    assert (d["rel_start"], d["rel_end"], d["hit_region"], d["hit_rbs"]) == (1, 30, "CDS", 0)
    assert bacterial.derive({}, 200, cfg) == {}


def test_command_line():
    cfg = bacterial.BacterialConfig(temperature=30.0, seed_bp=6)
    cmd = bacterial.build_command("ACGTTGCA", "A" * 300, cfg)
    assert cmd[0] == "IntaRNA" and cmd[cmd.index("-q") + 1] == "ACGUUGCA"
    assert cmd[cmd.index("--seedBP") + 1] == "6" and "--noSeed" not in cmd
    assert cmd[cmd.index("--temperature") + 1] == "30.0" and cmd[cmd.index("--threads") + 1] == "1"
    assert cmd[cmd.index("--tAccW") + 1] == "150" and cmd[cmd.index("--tAccL") + 1] == "100"
    assert "--seedQRange" not in cmd and "--tRegion" not in cmd     # whole window, any position
    short = bacterial.build_command("ACGTTGCA", "A" * 60, cfg)       # clipped window
    assert short[short.index("--tAccW") + 1] == "60" and short[short.index("--tAccL") + 1] == "60"
    cmd = bacterial.build_command("ACGTTGCA", "A" * 300, bacterial.BacterialConfig(no_seed=True))
    assert "--noSeed" in cmd and "--seedBP" not in cmd


# ---- scoring and q-values with IntaRNA replaced ------------------------------------------

SRNA = "TGCCAAGGCTTACGATCCGTAA"      # 22 nt
PLANTED = 7                          # index of the window carrying the planted interaction


def make_windows(n=40, anchor=200, length=300):
    regions = []
    for i in range(n):
        seq = rnd(length, 100 + i)
        regions.append(Region(f"g{i}|start_window", f"g{i}", f"gene{i}", "start_window", seq,
                              "c1", "+", [(i * 1000, i * 1000 + length)], anchor, 0, 100))
    return regions


def fake_window_score(guide, window, cfg):
    """Background best E around -9 for everyone; the real sRNA binds the planted window at the
    RBS (window positions 180-191) with -28."""
    jitter = (sum(map(ord, guide)) * 31 + sum(map(ord, window[:40])) * 17) % 40 / 10
    start, end, e = 20, 31, -9.0 - jitter
    if guide == SRNA and window.startswith(PLANTED_PREFIX):
        start, end, e = 180, 191, -28.0
    return {"intarna_E": e, "intarna_ED_target": 1.5, "intarna_ED_query": 0.5,
            "intarna_Pu_target": 0.1, "intarna_Pu_query": 0.2, "intarna_E_hybrid": e - 2,
            "intarna_t_start": start, "intarna_t_end": end, "intarna_q_start": 1,
            "intarna_q_end": 12, "intarna_duplex": "((((((&))))))"}


PLANTED_PREFIX = make_windows()[PLANTED].sequence[:40]


@pytest.fixture
def fake_intarna(monkeypatch):
    monkeypatch.setattr(bacterial, "score_window", fake_window_score)
    monkeypatch.setattr(bacterial, "check_tools", lambda cfg: None)


def run_fake(tmp_path, regions=None, n_decoys=20, cfg=None, **kw):
    decoys = make_decoys([("s1", SRNA)], n_decoys, (2, 8), 5)
    return bacterial.run_window_scoring([("s1", SRNA)], decoys, regions or make_windows(),
                                        cfg or bacterial.BacterialConfig(), str(tmp_path),
                                        log=lambda *_: None, **kw)


def test_planted_window_is_supported_and_the_others_are_not(fake_intarna, tmp_path):
    rows, targets, curves = run_fake(tmp_path)
    assert len(rows) == 40
    top = rows[0]
    assert top["transcript"] == f"g{PLANTED}" and top["intarna_E"] == -28.0
    assert top["intarna_q"] <= 0.05 and top["intarna_q_cand"] <= 0.05 and top["pass_fdr"] == 1
    assert sum(r["pass_fdr"] for r in rows) == 1
    assert all(r["intarna_q"] > 0.05 for r in rows[1:])
    # coordinates relative to the start codon (A = +1 at window index 200) and the part hit
    assert (top["rel_start"], top["rel_end"]) == (-21, -10)
    assert top["hit_region"] == "5UTR" and top["hit_rbs"] == 1
    assert top["genomic"] == f"c1:{PLANTED * 1000 + 180}-{PLANTED * 1000 + 191}(+)"
    assert "support" not in rows[0]
    assert targets[0]["gene"] == f"gene{PLANTED}" and targets[0]["windows_q_pass"] == 1
    assert {c["predictor"] for c in curves} == {"intarna"}
    assert all(0 <= c["fdr"] <= 1 for c in curves)


def test_the_unit_is_one_best_interaction_per_window(fake_intarna, tmp_path, monkeypatch):
    calls = []

    def counting(guide, window, cfg):
        calls.append((guide, window))
        return fake_window_score(guide, window, cfg)

    monkeypatch.setattr(bacterial, "score_window", counting)
    regions = make_windows(10)
    rows, _, curves = run_fake(tmp_path, regions, n_decoys=5)
    assert len(calls) == 10 * 6                                  # (1 + 5 decoys) x 10 windows
    # each real window enters the FDR once: the largest "real" count on the curve is 10
    assert max(c["real_sites"] for c in curves) <= 10


def test_decoys_that_score_better_give_high_q(monkeypatch, tmp_path):
    def reversed_scores(guide, window, cfg):
        out = fake_window_score(guide, window, cfg)
        out["intarna_E"] = -10.0 if guide == SRNA else -30.0
        return out

    monkeypatch.setattr(bacterial, "score_window", reversed_scores)
    monkeypatch.setattr(bacterial, "check_tools", lambda cfg: None)
    rows, _, _ = run_fake(tmp_path)
    assert all(r["intarna_q"] > 0.05 and r["pass_fdr"] == 0 for r in rows)


def test_identical_windows_are_scored_once_and_listed_for_each_gene(fake_intarna, tmp_path,
                                                                    monkeypatch):
    calls = []
    monkeypatch.setattr(bacterial, "score_window",
                        lambda g, w, c: calls.append(1) or fake_window_score(g, w, c))
    regions = make_windows(5)
    copy = Region("dup|start_window", "dup", "dupgene", "start_window", regions[0].sequence,
                  "c2", "-", [(0, 300)], 200, 0, 100)
    rows, _, _ = run_fake(tmp_path, regions + [copy], n_decoys=3)
    assert len(rows) == 6 and len(calls) == 5 * 4
    assert {r["transcript"] for r in rows if r["window_length"] == 300} >= {"g0", "dup"}


def test_candidates_without_decoys_get_no_q_and_missing_decoys_fail(fake_intarna, tmp_path):
    decoys = make_decoys([("s1", SRNA)], 5, (2, 8), 5)
    msgs = []
    rows, _, _ = bacterial.run_window_scoring(
        [("s1", SRNA), ("s2", rnd(22, 3))], decoys, make_windows(5), bacterial.BacterialConfig(),
        str(tmp_path), log=msgs.append)
    assert all(r["intarna_q"] == "" for r in rows if r["candidate"] == "s2")
    assert any("no decoys" in m for m in msgs)
    with pytest.raises(ValueError, match="needs decoys"):
        bacterial.run_window_scoring([("s1", SRNA)], [], make_windows(3),
                                     bacterial.BacterialConfig(), str(tmp_path))


def test_regions_without_anchor_are_refused(fake_intarna, tmp_path):
    bad = [Region("r", "t", "g", "three_prime_UTR", rnd(100, 1), "c", "+", [(0, 100)])]
    with pytest.raises(ValueError, match="start-codon window"):
        run_fake(tmp_path, bad)


def test_checkpoint_resume_and_changed_settings(fake_intarna, tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(bacterial, "score_window",
                        lambda g, w, c: calls.append(1) or fake_window_score(g, w, c))
    first, _, _ = run_fake(tmp_path, n_decoys=4)
    n = len(calls)
    assert n == 40 * 5
    second, _, _ = run_fake(tmp_path, n_decoys=4)
    assert len(calls) == n and [r["intarna_q"] for r in second] == [r["intarna_q"] for r in first]
    run_fake(tmp_path, n_decoys=4, cfg=bacterial.BacterialConfig(temperature=25.0))
    assert len(calls) == 2 * n                       # another temperature: other signature
    run_fake(tmp_path, n_decoys=4, cfg=bacterial.BacterialConfig(resume=False))
    assert len(calls) == 3 * n


def test_interrupted_run_resumes_from_the_partial_checkpoint(fake_intarna, tmp_path, monkeypatch):
    calls = []

    def dies_late(guide, window, cfg):
        calls.append(1)
        if len(calls) > 120:
            raise KeyboardInterrupt        # not an Exception: the run aborts, like a kill
        return fake_window_score(guide, window, cfg)

    monkeypatch.setattr(bacterial, "score_window", dies_late)
    cfg = bacterial.BacterialConfig(checkpoint_seconds=0)
    with pytest.raises(KeyboardInterrupt):
        run_fake(tmp_path, n_decoys=4, cfg=cfg)
    monkeypatch.setattr(bacterial, "score_window",
                        lambda g, w, c: calls.append(1) or fake_window_score(g, w, c))
    before = len(calls)
    rows, _, _ = run_fake(tmp_path, n_decoys=4, cfg=cfg)
    assert 0 < len(calls) - before < 40 * 5          # only the missing part was scored
    assert rows[0]["transcript"] == f"g{PLANTED}"


def test_a_failing_window_is_left_unscored_not_fatal(monkeypatch, tmp_path):
    state = {"n": 0}

    def flaky(guide, window, cfg):
        state["n"] += 1
        if state["n"] == 3:
            raise RuntimeError("IntaRNA failed (1): boom")
        return fake_window_score(guide, window, cfg)

    monkeypatch.setattr(bacterial, "score_window", flaky)
    monkeypatch.setattr(bacterial, "check_tools", lambda cfg: None)
    msgs = []
    decoys = make_decoys([("s1", SRNA)], 20, (2, 8), 5)
    rows, _, _ = bacterial.run_window_scoring(
        [("s1", SRNA)], decoys, make_windows(), bacterial.BacterialConfig(), str(tmp_path),
        log=msgs.append)
    assert any("failed and are left unscored" in m for m in msgs)
    assert len(rows) == 40


def test_too_many_failures_stop_the_run(monkeypatch, tmp_path):
    def broken(guide, window, cfg):
        raise RuntimeError("tool crashed")

    monkeypatch.setattr(bacterial, "score_window", broken)
    monkeypatch.setattr(bacterial, "check_tools", lambda cfg: None)
    with pytest.raises(RuntimeError, match="windows failed"):
        run_fake(tmp_path)


def test_dry_run_and_missing_binary(fake_intarna, tmp_path, capsys):
    decoys = make_decoys([("s1", SRNA)], 5, (2, 8), 5)
    assert bacterial.run_window_scoring([("s1", SRNA)], decoys, make_windows(3),
                                        bacterial.BacterialConfig(), str(tmp_path),
                                        dry_run=True) is None
    assert "IntaRNA runs" in capsys.readouterr().out


def test_missing_binary_is_reported_before_scoring():
    with pytest.raises(SystemExit, match="IntaRNA_not_installed"):
        bacterial.check_tools(bacterial.BacterialConfig(intarna_bin="IntaRNA_not_installed"))


def test_write_windows_roundtrip(fake_intarna, tmp_path):
    rows, targets, curves = run_fake(tmp_path)
    out = [str(tmp_path / f) for f in ("w.tsv", "t.tsv", "c.tsv")]
    bacterial.write_windows(rows, targets, curves, *out)
    header = (tmp_path / "w.tsv").read_text().splitlines()[0].split("\t")
    for col in ("intarna_E", "intarna_q", "intarna_q_cand", "rel_start", "hit_region", "hit_rbs",
                "upstream_cds"):
        assert col in header
    assert "support" not in header and "rnahybrid_mfe" not in header
    assert len((tmp_path / "w.tsv").read_text().splitlines()) == len(rows) + 1
    # validate-targets reads the same columns
    from smithhunter.seqio import read_tsv
    from smithhunter.targets.validate import validate
    _, summary = validate(read_tsv(out[1]), [{"candidate": "s1", "gene": f"gene{PLANTED}"}])
    assert summary["auc"]["intarna"] == 1.0


def test_cli_regions_and_score_windows_with_a_faked_intarna(tmp_path, monkeypatch, capsys):
    genome, lines, truth = simulate_genome(12, seed=5)
    fa, gff = write_genome(tmp_path, {"chr1": genome}, lines)
    cand = tmp_path / "candidates.fasta"
    cand.write_text(f">s1{chr(10)}{SRNA_REAL}{chr(10)}")
    out = lambda name: str(tmp_path / name)
    main(["regions", "--mode", "bacterial", "--genome", fa, "--annotation", gff,
          "--window=-150,60", "--circular-contigs", "chr1", "--out", out("regions.tsv")])
    assert "12 start-codon windows" in capsys.readouterr().out
    regions = read_regions(out("regions.tsv"))
    assert all(len(r.sequence) == 210 and r.anchor == 150 for r in regions)
    with pytest.raises(SystemExit):
        main(["regions", "--mode", "bacterial", "--genome", fa, "--annotation", gff,
              "--window=100,-50", "--out", out("x.tsv")])
    main(["decoys", "--candidates", str(cand), "--n", "5", "--out", out("decoys.fasta")])
    planted = "L_" + truth["gene"]
    seen = []

    def fake(guide, window, cfg):
        seen.append(cfg)
        planted_win = next(r.sequence for r in regions if r.transcript == planted)
        e = -30.0 if (guide == SRNA_REAL and window == planted_win) else -8.0
        return {"intarna_E": e, "intarna_t_start": 100, "intarna_t_end": 118,
                "intarna_q_start": 1, "intarna_q_end": 19, "intarna_duplex": "(&)"}

    monkeypatch.setattr(bacterial, "score_window", fake)
    monkeypatch.setattr(bacterial, "check_tools", lambda cfg: None)
    main(["score-windows", "--candidates", str(cand), "--decoys", out("decoys.fasta"),
          "--regions", out("regions.tsv"), "--temperature", "30", "--seed-bp", "6",
          "--acc-window", "120", "--workdir", out("work"), "--windows", out("w.tsv"),
          "--targets", out("t.tsv"), "--curves", out("c.tsv"), "--fdr", "0.2"])
    assert seen[0].temperature == 30.0 and seen[0].seed_bp == 6 and seen[0].acc_window == 120
    from smithhunter.seqio import read_tsv
    rows = read_tsv(out("w.tsv"))
    assert rows[0]["transcript"] == planted and rows[0]["pass_fdr"] == "1"
    assert rows[0]["hit_region"] == "5UTR" and rows[0]["rel_start"] == "-51"
    assert sum(int(r["pass_fdr"]) for r in rows) == 1
    assert len(read_tsv(out("t.tsv"))) == 12 and read_tsv(out("c.tsv"))


# ---- with the real IntaRNA ----------------------------------------------------------------

needs_intarna = pytest.mark.skipif(shutil.which("IntaRNA") is None, reason="IntaRNA is not on PATH")

SRNA_REAL = "GCAUGCCUUGAACGGAUAGCUC".replace("U", "T")     # 22 nt


@needs_intarna
def test_real_intarna_finds_a_planted_rbs_site():
    site = revcomp(SRNA_REAL[:12])                           # pairs with the 5' end of the sRNA
    window = rnd(175, 1) + site + rnd(3, 2) + "AGGAGG" + rnd(4, 3) + "ATG" + rnd(97, 4)
    assert window.index("ATG", 195) == 200 and len(window) == 300
    cfg = bacterial.BacterialConfig(temperature=37.0)
    res = bacterial.score_window(SRNA_REAL, window, cfg)
    assert res["intarna_E"] < -15
    d = bacterial.derive(res, 200, cfg)
    assert d["hit_region"] == "5UTR" and d["hit_rbs"] == 1
    other = bacterial.score_window(SRNA_REAL, rnd(300, 5), cfg)
    assert other is None or other["intarna_E"] > res["intarna_E"] + 8


# ---- synthetic genome, end to end ---------------------------------------------------------

def sense_codons(rng, n):
    stops = {"TAA", "TAG", "TGA"}
    out = []
    while len(out) < n:
        c = "".join(rng.choice("ACGT") for _ in range(3))
        if c not in stops:
            out.append(c)
    return "".join(out)


PLANT_LEN = 16


def simulate_genome(n_genes=30, seed=1, srna=SRNA_REAL, plant=("last", -34), plant_len=PLANT_LEN):
    """A circular bacterial-like genome on both strands. Gene 0 (+) starts 20 nt after the origin
    and the last gene (-) starts 20 nt before the end, so both of their windows cross the origin.
    The sRNA's 5' ``plant_len`` nt are planted (as the complementary target, in the orientation of the mRNA)
    in the window of one gene at relative positions ``plant[1]``..``plant[1] + plant_len - 1``.
    Returns (genome, gff lines, truth) where truth gives gene, strand, genomic 0-based positions
    of the planted target in genome order."""
    rng = random.Random(seed)
    parts, features, pos = [rnd(20, seed + 1000)], [], 20      # leading spacer
    genes = []
    for i in range(n_genes):
        strand = "+" if i == 0 else "-" if i == n_genes - 1 else rng.choice("+-")
        body = "ATG" + sense_codons(rng, rng.randint(40, 120)) + "TAA"
        if i == n_genes - 1:
            gap = ""
        else:
            gap = rnd(rng.randint(120, 250), rng.randint(0, 10**6))
        genes.append((strand, body, gap))
    # lay the genes out (mRNA orientation on the + strand for '+', reverse complement for '-')
    for i, (strand, body, gap) in enumerate(genes):
        seg = body if strand == "+" else revcomp(body)
        start = pos
        parts.append(seg)
        pos += len(seg)
        features.append((f"G{i:03d}", start, pos, strand))
        parts.append(gap)
        pos += len(gap)
    # last gene: its start codon (high coordinate, '-') must sit 20 nt before the end
    parts.append(rnd(20, seed + 2000))
    genome = list("".join(parts))
    L = len(genome)
    gene = features[-1] if plant[0] == "last" else features[plant[0]]
    name, s, e, strand = gene
    g0 = s if strand == "+" else e - 1                             # 0-based start codon A
    target = revcomp(srna[:plant_len])
    truth = []
    for k, base in enumerate(target):
        rel = plant[1] + k                                         # negative: upstream
        assert rel < 0
        if strand == "+":
            gp, b = (g0 + rel) % L, base
        else:
            gp, b = (g0 - rel) % L, COMP[base]
        genome[gp] = b
        truth.append(gp)
    lines = []
    for name, s, e, strand in features:
        lines += cds("chr1", s + 1, e, strand, name)
    return "".join(genome), lines, {"gene": name, "strand": strand, "positions": truth,
                                    "length": L, "g0": g0}


def test_simulated_genome_layout(tmp_path):
    genome, lines, truth = simulate_genome(12, seed=3)
    fa, gff = write_genome(tmp_path, {"chr1": genome}, lines)
    regions = bacterial_windows(fa, gff, (-200, 100), circular=True)
    assert len(regions) == 12
    for r in regions:
        assert len(r.sequence) == 300 and r.sequence[200:203] == "ATG"
        check_mapping(r, genome, True)
    wrapped = [r for r in regions if len(r.blocks) > 1]
    assert {r.transcript for r in wrapped} >= {"L_G000", "L_G011"}
    planted = [r for r in regions if r.transcript == "L_" + truth["gene"]][0]
    assert revcomp(SRNA_REAL[:PLANT_LEN]) in planted.sequence[:200]       # in the 5' part, mRNA sense


@needs_intarna
def test_end_to_end_on_a_simulated_genome(tmp_path, capsys):
    n_genes = int(os.environ.get("SMITHHUNTER_BACT_GENES", "30"))
    genome, lines, truth = simulate_genome(n_genes, seed=11)
    fa, gff = write_genome(tmp_path, {"chr1": genome}, lines)
    cand = tmp_path / "candidates.fasta"
    cand.write_text(f">srna_planted\n{SRNA_REAL}\n>srna_random\n{rnd(22, 77)}\n")
    out = lambda name: str(tmp_path / name)
    workers = os.environ.get("SMITHHUNTER_BACT_WORKERS", "2")
    main(["regions", "--mode", "bacterial", "--genome", fa, "--annotation", gff,
          "--window=-200,100", "--topology", "circular", "--out", out("regions.tsv")])
    main(["decoys", "--candidates", str(cand), "--n", "20", "--out", out("decoys.fasta")])
    main(["score-windows", "--candidates", str(cand), "--decoys", out("decoys.fasta"),
          "--regions", out("regions.tsv"), "--temperature", "37", "--workers", workers,
          "--workdir", out("work"), "--windows", out("windows.tsv"),
          "--targets", out("targets.tsv"), "--curves", out("curves.tsv")])
    from smithhunter.seqio import read_tsv
    rows = read_tsv(out("windows.tsv"))
    planted = [r for r in rows if r["candidate"] == "srna_planted"]
    assert len(planted) == n_genes
    top = planted[0]
    assert top["transcript"] == "L_" + truth["gene"], top
    assert float(top["intarna_q"]) <= 0.05 and top["pass_fdr"] == "1"
    assert top["hit_region"] == "5UTR" and top["hit_rbs"] == "1"
    assert -45 <= int(top["rel_start"]) and int(top["rel_end"]) <= -5
    # the interaction maps back to the planted positions in the genome; the planted target
    # straddles the origin, so the reported span runs clockwise with start > end
    L = truth["length"]
    a, b = (int(x) for x in top["genomic"].split(":")[1].split("(")[0].split("-"))
    assert truth["strand"] == "-" and top["genomic"].endswith("(-)") and a > b
    covered = {(a - 1 + k) % L + 1 for k in range((b - a) % L + 1)}
    assert len(covered & {p + 1 for p in truth["positions"]}) >= 7
    n_pass = sum(int(r["pass_fdr"]) for r in planted)
    n_rand = sum(int(r["pass_fdr"]) for r in rows if r["candidate"] == "srna_random")
    print(f"\nsimulated genome: {n_genes} genes, planted candidate: {n_pass} window(s) with "
          f"q <= 0.05 (planted q = {top['intarna_q']}, E = {top['intarna_E']}); "
          f"random candidate: {n_rand}")
    assert n_pass <= 1 + max(1, n_genes // 20) and n_rand <= max(1, n_genes // 20)
    curves = read_tsv(out("curves.tsv"))
    assert curves and all(c["site_class"] == "window" for c in curves)
