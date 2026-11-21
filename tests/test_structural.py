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
    (Feature(0, 100, "+", "exon", "exon-1", "16S ribosomal RNA"), True),
    (Feature(0, 900, "+", "CDS", "rsmA", "16S rRNA (adenine) methyltransferase"), False),
    (Feature(0, 900, "+", "CDS", "ffh", "signal recognition particle protein"), False),
    (Feature(0, 200, "+", "ncRNA", "rna-2", "small regulatory RNA RsmZ"), False),
    (Feature(0, 900, "+", "gene", "ftsZ"), False),
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
                 "c1\tRefSeq\tncRNA\t1\t180\t.\t+\t.\tID=rna-1;product=6S%20RNA\n")
    f = Annotation.from_gff(str(p)).overlaps("c1", 0, 10)[0]
    assert f.product == "6S RNA" and is_structural(f)


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


def test_exclude_structural_rejects_unknown_mode(tmp_path):
    with pytest.raises(ValueError):
        discover(*inputs(tmp_path), exclude_structural="tRNA")
