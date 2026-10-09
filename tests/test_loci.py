import pytest

from smithhunter.annotate import Annotation, Feature, describe
from smithhunter.ends import end_score
from smithhunter.loci import Member, build_loci, discover


def test_end_score():
    assert end_score({10: 100}) == (1.0, 1.0)
    # the top position alone holds half of the reads: it is enough, as in v0
    assert end_score({10: 60, 11: 40}) == (1.0, 0.6)
    score, dominant = end_score({10: 40, 11: 35, 12: 25})
    assert score == pytest.approx(0.5) and dominant == pytest.approx(0.4)
    # two kept positions three bases apart: 2 positions + 2 gaps * 0.1
    score, _ = end_score({10: 30, 13: 30, 20: 1})
    assert score == pytest.approx(1 / 2.2)
    assert end_score({}) == (0.0, 0.0)


def test_build_loci_merges_overlaps_only():
    members = [Member("a", 10, 30, 1), Member("b", 25, 45, 1), Member("c", 50, 70, 1)]
    loci = build_loci("mt", "chrM", "+", members, 0, 1000, False)
    assert [(l.start, l.end, len(l.members)) for l in loci] == [(10, 45, 2), (50, 70, 1)]
    loci = build_loci("mt", "chrM", "+", members, 5, 1000, False)
    assert [(l.start, l.end) for l in loci] == [(10, 70)]


def test_build_loci_joins_across_circular_origin():
    members = [Member("a", 90, 110, 1), Member("b", 2, 22, 1), Member("c", 40, 60, 1)]
    loci = build_loci("mt", "chrM", "+", members, 0, 100, True)
    spans = sorted((l.start, l.end, l.wraps) for l in loci)
    assert spans == [(40, 60, False), (90, 122, True)]


def test_describe_picks_most_specific_class():
    feats = [Feature(0, 100, "+", "gene", "COX1"), Feature(10, 80, "-", "tRNA", "trnS")]
    d = describe(feats, "+")
    assert d["annotation_class"] == "tRNA" and d["annotation_orientation"] == "antisense"
    assert describe([], "+")["annotation_class"] == "intergenic"


def test_annotation_overlap_query():
    ann = Annotation({"chrM": [Feature(0, 70, "+", "tRNA", "trnM"),
                               Feature(100, 200, "+", "CDS", "ND2")]})
    assert [f.name for f in ann.overlaps("chrM", 60, 110)] == ["trnM", "ND2"]
    assert ann.overlaps("chrM", 70, 100) == []
    assert ann.overlaps("other", 0, 10) == []


def write(path, text):
    path.write_text(text)
    return str(path)


def test_discover_end_to_end(tmp_path):
    counts = write(tmp_path / "counts.tsv",
                   "read_id\tsequence\tlength\ts1\ts2\n"
                   "u1\tACGTACGTACGTACGTACGT\t20\t80\t90\n"
                   "u2\tACGTACGTACGTACGTACG\t19\t20\t10\n"
                   "u3\tGGGGGGGGGGGGGGGGGGGG\t20\t1\t0\n")
    libs = write(tmp_path / "libs.tsv", "sample\treads_in\treads_kept\ns1\t200\t101\ns2\t200\t100\n")
    gff = write(tmp_path / "mt.gff3", "chrM\tmitos\ttRNA\t1\t60\t.\t+\t.\tName=trnM\n")
    genomes = write(tmp_path / "genomes.tsv",
                    "name\tfasta\trole\tcircular\tannotation\tmin_length\tmax_length\ttargets\n"
                    f"mt\tmt.fa\tfocus\t1\t{gff}\t18\t35\tago\n")
    contigs = write(tmp_path / "contigs.tsv",
                    "ref\tgenome\tcontig\tlength\tcircular\nmt|chrM\tmt\tchrM\t1000\t1\n")
    alignments = write(tmp_path / "aln.tsv",
                       "read_id\tgenome\tcontig\tstart\tend\tstrand\tweight\tmismatches\n"
                       "u1\tmt\tchrM\t10\t30\t+\t1\t0\n"
                       "u2\tmt\tchrM\t10\t29\t+\t1\t0\n"
                       "u3\tmt\tchrM\t500\t520\t-\t1\t0\n")
    samples, rows, reads = discover(alignments, counts, libs, genomes, contigs,
                                    min_rpm=5, min_samples=2)
    assert samples == ["s1", "s2"]
    first, second = rows
    assert (first["start"], first["end"], first["strand"]) == (11, 30, "+")
    assert first["rep_id"] == "u1" and first["count_s1"] == 100
    assert first["five_score"] == 1.0  # every read starts at position 11
    assert first["annotation_class"] == "tRNA" and first["annotation_orientation"] == "sense"
    assert first["pass"] == 1
    assert second["annotation_class"] == "intergenic"
    assert second["expression_pass"] == 0  # only one read
    assert len(reads) == 3
