import gzip

import pytest

from smithhunter.collapse import collapse
from smithhunter.reference import build_reference, parse_genomes
from smithhunter.seqio import read_fasta


def write_fastq(path, seqs):
    with gzip.open(path, "wt") as fh:
        for i, s in enumerate(seqs):
            fh.write(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n")


def test_collapse_counts_per_sample_and_filters_length(tmp_path):
    a, b = tmp_path / "a.fq.gz", tmp_path / "b.fq.gz"
    write_fastq(a, ["ACGTACGTACGTACGTACGT"] * 3 + ["ACGT"] + ["TTTTTTTTTTTTTTTTTTTT"])
    write_fastq(b, ["ACGTACGTACGTACGTACGT", "ACGTNCGTACGTACGTACGT"])
    names, ordered, stats = collapse({"a": [str(a)], "b": [str(b)]}, 18, 35)
    assert names == ["a", "b"]
    assert ordered[0] == ("ACGTACGTACGTACGTACGT", [3, 1])
    assert ordered[1] == ("TTTTTTTTTTTTTTTTTTTT", [1, 0])
    assert len(ordered) == 2
    assert stats[0] == {"sample": "a", "reads_in": 5, "reads_kept": 4}
    assert stats[1]["reads_kept"] == 1


def test_parse_genomes_validates_roles():
    with pytest.raises(ValueError, match="at least one genome must have role 'focus'"):
        parse_genomes({"nuc": {"fasta": "x.fa", "role": "exclude"}}, 18, 35)
    with pytest.raises(ValueError, match="role must be one of"):
        parse_genomes({"nuc": {"fasta": "x.fa", "role": "study"}}, 18, 35)
    g = parse_genomes({"mt": {"fasta": "m.fa", "topology": "circular", "length": [20, 40]}}, 18, 35)
    assert g[0].circular and (g[0].min_length, g[0].max_length) == (20, 40)


def test_build_reference_prefixes_and_extends_circular(tmp_path):
    mt, nuc = tmp_path / "mt.fa", tmp_path / "nuc.fa"
    mt.write_text(">chrM desc\nACGTTGCA\n")
    nuc.write_text(">chr1\nGGGG\nCCCC\n")
    genomes = parse_genomes({"mt": {"fasta": str(mt), "topology": "circular"},
                             "nuc": {"fasta": str(nuc), "role": "exclude"}}, 18, 35)
    out = tmp_path / "combined.fa"
    contigs = build_reference(genomes, str(out), overhang=3)
    seqs = dict(read_fasta(str(out)))
    assert seqs == {"mt|chrM": "ACGTTGCAACG", "nuc|chr1": "GGGGCCCC"}
    assert contigs[0] == {"ref": "mt|chrM", "genome": "mt", "contig": "chrM", "length": 8,
                          "circular": 1}
