"""Structural RNA classification and the loci.exclude_structural filter."""
import pytest

from smithhunter.annotate import Annotation, Feature, describe, is_structural, structural
from smithhunter.loci import discover


def test_prokaryotic_rna_types_are_classes_not_exon():
    # NCBI prokaryotic GFF: gene -> RNase_P_RNA -> exon; tmRNA likewise
    for rna in ("RNase_P_RNA", "tmRNA", "SRP_RNA", "antisense_RNA"):
        feats = [Feature(0, 300, "+", "gene", "rnpB"), Feature(0, 300, "+", rna, "rna-x"),
                 Feature(0, 300, "+", "exon", "exon-x")]
        assert describe(feats, "+")["annotation_class"] == rna


@pytest.mark.parametrize("feature, expected", [
    (Feature(0, 70, "+", "tRNA", "trnM"), True),
    (Feature(0, 1500, "+", "rRNA", "rrs"), True),
    (Feature(0, 360, "+", "tmRNA", "ssrA"), True),
    (Feature(0, 180, "+", "ncRNA", "ssrS"), True),                  # 6S RNA by gene name
    (Feature(0, 180, "+", "ncRNA", "rna-1", "6S RNA"), True),      # by product
    (Feature(0, 100, "+", "exon", "exon-1", "16S ribosomal RNA", gbkey="rRNA"), True),
    (Feature(0, 70, "+", "exon", "exon-2", "tRNA-Leu", gbkey="tRNA"), True),
    (Feature(0, 300, "+", "ncRNA", "rna-3", "", ncrna_class="RNase_MRP_RNA"), True),
    (Feature(0, 300, "+", "ncRNA", "rna-4", "RNase MRP RNA"), True),
    (Feature(0, 300, "+", "gene", "7SLRNA:CR42652"), True),               # FlyBase names
    (Feature(0, 300, "+", "gene", "snRNA:7SK"), True),
    (Feature(0, 900, "+", "CDS", "rsmA", "16S rRNA (adenine) methyltransferase"), False),
    (Feature(0, 900, "+", "CDS", "ffh", "signal recognition particle protein"), False),
    (Feature(0, 200, "+", "ncRNA", "rna-2", "small regulatory RNA RsmZ"), False),
    (Feature(0, 900, "+", "gene", "ftsZ"), False),
    # eukaryotic mRNA exons and transcripts carry the protein product (review finding)
    (Feature(0, 300, "+", "exon", "Srp54", "signal recognition particle 54 kDa protein, "
             "transcript variant X1", gbkey="mRNA"), False),
    (Feature(0, 300, "+", "mRNA", "Rrp1", "ribosomal RNA processing protein 1 homolog",
             gbkey="mRNA"), False),
    (Feature(0, 300, "+", "exon", "x", "ribosomal RNA small subunit methyltransferase"), False),
    (Feature(0, 300, "+", "exon", "smpB", "tmRNA-binding protein SmpB"), False),
    (Feature(0, 70, "+", "lnc_RNA", "lncRNA:CR1", "tRNA-derived lncRNA"), False),
])
def test_is_structural(feature, expected):
    assert is_structural(feature) is expected


def test_structural_orientation():
    trna_plus = Feature(10, 80, "+", "tRNA", "trnS")
    cds_minus = Feature(0, 900, "-", "CDS", "x")
    assert structural([trna_plus, cds_minus], "+") == "sense"
    assert structural([trna_plus, cds_minus], "-") == "antisense"
    assert structural([cds_minus], "+") == ""
    assert structural([Feature(10, 80, ".", "rRNA", "rrf")], "-") == "sense"   # unstranded
    both = [trna_plus, Feature(10, 80, "-", "tRNA", "trnL")]
    assert structural(both, "-") == "sense"


def test_gff_product_is_read_and_unquoted(tmp_path):
    p = tmp_path / "a.gff3"
    p.write_text("##gff-version 3\n"
                 "c1\tRefSeq\tncRNA\t1\t180\t.\t+\t.\tID=rna-1;product=6S%20RNA\n"
                 "c1\tRefSeq\tncRNA\t300\t600\t.\t+\t.\tID=rna-2;gbkey=ncRNA;ncrna_class=SRP_RNA\n"
                 "c1\tRefSeq\texon\t700\t900\t.\t+\t.\tID=exon-3;gbkey=mRNA;"
                 "product=signal%20recognition%20particle%2054%20kDa%20protein\n")
    ann = Annotation.from_gff(str(p))
    f = ann.overlaps("c1", 0, 10)[0]
    assert f.product == "6S RNA" and is_structural(f)
    g = ann.overlaps("c1", 400, 410)[0]
    assert g.gbkey == "ncRNA" and g.ncrna_class == "SRP_RNA" and is_structural(g)
    assert not is_structural(ann.overlaps("c1", 800, 810)[0])


def write(path, text):
    path.write_text(text)
    return str(path)


def inputs(tmp_path):
    counts = write(tmp_path / "counts.tsv",
                   "read_id\tsequence\tlength\ts1\ts2\n"
                   "u1\tACGTACGTACGTACGTACGT\t20\t80\t90\n"
                   "u2\tACGTACGTACGTACGTACG\t19\t20\t10\n"
                   "u3\tGGGGGGGGGGGGGGGGGGGG\t20\t60\t70\n")
    libs = write(tmp_path / "libs.tsv", "sample\treads_in\treads_kept\ns1\t200\t101\ns2\t200\t100\n")
    gff = write(tmp_path / "b.gff3",
                "c1\tRefSeq\ttRNA\t1\t60\t.\t+\t.\tName=trnM\n"           # u1/u2 sense
                "c1\tRefSeq\trRNA\t490\t600\t.\t+\t.\tName=rrf\n")       # u3 antisense
    genomes = write(tmp_path / "genomes.tsv",
                    "name\tfasta\trole\tcircular\tannotation\tmin_length\tmax_length\ttargets\n"
                    f"bact\tb.fa\tfocus\t1\t{gff}\t18\t35\tnone\n")
    contigs = write(tmp_path / "contigs.tsv",
                    "ref\tgenome\tcontig\tlength\tcircular\nbact|c1\tbact\tc1\t1000\t1\n")
    alignments = write(tmp_path / "aln.tsv",
                       "read_id\tgenome\tcontig\tstart\tend\tstrand\tweight\tmismatches\n"
                       "u1\tbact\tc1\t10\t30\t+\t1\t0\n"
                       "u2\tbact\tc1\t10\t29\t+\t1\t0\n"
                       "u3\tbact\tc1\t500\t520\t-\t1\t0\n")
    return alignments, counts, libs, genomes, contigs


@pytest.mark.parametrize("mode, passing", [("none", [1, 1]), ("sense", [0, 1]), ("any", [0, 0])])
def test_exclude_structural_modes(tmp_path, mode, passing):
    _, rows, _ = discover(*inputs(tmp_path), min_rpm=5, min_samples=2, exclude_structural=mode)
    assert [r["structural"] for r in rows] == ["sense", "antisense"]
    assert [r["expression_pass"] and r["ends_pass"] for r in rows] == [1, 1]
    assert [r["pass"] for r in rows] == passing
    assert [r["structural_pass"] for r in rows] == passing


def test_no_annotation_means_not_structural(tmp_path):
    alignments, counts, libs, genomes, contigs = inputs(tmp_path)
    text = open(genomes).read().split("\n")
    fields = text[1].split("\t")
    fields[4] = ""                       # no annotation
    text[1] = "\t".join(fields)
    open(genomes, "w").write("\n".join(text))
    _, rows, _ = discover(alignments, counts, libs, genomes, contigs, min_rpm=5, min_samples=2,
                          exclude_structural="any")
    assert [r["structural"] for r in rows] == ["", ""]
    assert [r["pass"] for r in rows] == [1, 1]
    assert {r["annotation_class"] for r in rows} == {"NA"}


def test_min_rpm_genome(tmp_path):
    # genome totals: s1 = 160 reads, s2 = 170; locus 1 has 100 + 100, locus 2 has 60 + 70
    args = inputs(tmp_path)
    _, rows, _ = discover(*args, min_rpm=5, min_samples=2, min_rpm_genome=0)
    assert [r["expression_pass"] for r in rows] == [1, 1]
    _, rows, _ = discover(*args, min_rpm=5, min_samples=2, min_rpm_genome=500_000)
    assert [r["expression_pass"] for r in rows] == [1, 0]
    assert rows[0]["rpm_genome_s1"] == 625000.0


def test_exact_fraction(tmp_path):
    alignments, counts, libs, genomes, contigs = inputs(tmp_path)
    text = open(alignments).read().replace("u2\tbact\tc1\t10\t29\t+\t1\t0", "u2\tbact\tc1\t10\t29\t+\t1\t1")
    text = text.replace("u3\tbact\tc1\t500\t520\t-\t1\t0", "u3\tbact\tc1\t500\t520\t-\t1\t1")
    open(alignments, "w").write(text)
    args = (alignments, counts, libs, genomes, contigs)
    _, rows, _ = discover(*args, min_rpm=5, min_samples=2)
    assert [r["exact_fraction"] for r in rows] == [round(170 / 200, 4), 0.0]   # u1 170 of 200
    assert [r["pass"] for r in rows] == [1, 1]                                 # off by default
    _, rows, _ = discover(*args, min_rpm=5, min_samples=2, min_exact_fraction=0.5)
    assert [r["exact_pass"] for r in rows] == [1, 0]
    assert [r["pass"] for r in rows] == [1, 0]


def test_circular_wrap_keeps_mismatches(tmp_path):
    from smithhunter.loci import Member, build_loci, fold_origin, Locus
    ms = [Member("a", 990, 1010, 1.0, 5.0, (5.0,), 1), Member("b", 0, 15, 1.0, 5.0, (5.0,), 0)]
    loci = build_loci("g", "c", "+", ms, 0, 1000, True)
    assert sorted(m.mismatches for l in loci for m in l.members) == [0, 1]
    shifted = fold_origin(Locus("g", "c", "+", 1005, 1020, [Member("x", 1005, 1020, 1.0, 1.0, (1.0,), 1)]),
                          1000, True)
    assert shifted.members[0].mismatches == 1 and shifted.start == 5


def test_min_samples_scope_group(tmp_path):
    alignments, counts, libs, genomes, contigs = inputs(tmp_path)
    text = open(counts).read().replace("u3\tGGGGGGGGGGGGGGGGGGGG\t20\t60\t70",
                                       "u3\tGGGGGGGGGGGGGGGGGGGG\t20\t60\t0")   # only in s1
    open(counts, "w").write(text)
    args = (alignments, counts, libs, genomes, contigs)
    _, rows, _ = discover(*args, min_rpm=5)                       # n - 1 = 1 of 2
    assert [r["expression_pass"] for r in rows] == [1, 1]
    _, rows, _ = discover(*args, min_rpm=5, min_samples=2)
    assert [r["expression_pass"] for r in rows] == [1, 0]
    groups = {"s1": "male", "s2": "female"}
    _, rows, _ = discover(*args, min_rpm=5, min_samples=2, groups=groups, min_samples_scope="group")
    assert [r["expression_pass"] for r in rows] == [1, 1]          # min(2, group size 1) in males
    assert rows[1]["samples_passing_rpm"] == 1
    with pytest.raises(ValueError):
        discover(*args, groups={"s1": "male"}, min_samples_scope="group")
    with pytest.raises(ValueError):
        discover(*args, min_samples_scope="sex")


def test_read_groups(tmp_path):
    from smithhunter.cli import read_groups
    p = tmp_path / "samples.tsv"
    p.write_text("sample\tfq1\tfq2\tgroup\nW1M1\ta.fq\t\tmale\nW1F1\tb.fq\t\tfemale\n")
    assert read_groups(str(p)) == {"W1M1": "male", "W1F1": "female"}
    assert read_groups(None) is None


def test_exclude_structural_rejects_unknown_mode(tmp_path):
    with pytest.raises(ValueError):
        discover(*inputs(tmp_path), exclude_structural="tRNA")
