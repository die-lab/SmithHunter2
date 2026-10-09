import pytest

from smithhunter.assign import Contig, Hit, classify, group_hits, parse_sam, reference_length

CONTIGS = {
    "mt|chrM": Contig("mt|chrM", "mt", "chrM", 100, True),
    "nuc|chr1": Contig("nuc|chr1", "nuc", "chr1", 1000, False),
    "bac|chr": Contig("bac|chr", "bac", "chr", 500, True),
}
ROLES = {"mt": "focus", "nuc": "exclude", "bac": "focus"}


def sam(read, flag, ref, pos, cigar="20M", nm=0):
    return f"{read}\t{flag}\t{ref}\t{pos}\t255\t{cigar}\t*\t0\t0\tACGT\tIIII\tXA:i:0\tNM:i:{nm}\n"


def test_reference_length():
    assert reference_length("20M") == 20
    assert reference_length("5S10M2D3M") == 15


def test_parse_sam_skips_unmapped_and_reads_strand():
    lines = ["@HD\tVN:1.0\n", sam("u1", 0, "mt|chrM", 11, nm=1), sam("u2", 4, "*", 0),
             sam("u3", 16, "nuc|chr1", 1)]
    hits = list(parse_sam(lines))
    assert hits == [("u1", Hit("mt|chrM", 10, 30, "+", 1)), ("u3", Hit("nuc|chr1", 0, 20, "-", 0))]


def test_group_hits_requires_contiguous_reads():
    lines = [sam("u1", 0, "mt|chrM", 1), sam("u2", 0, "mt|chrM", 5), sam("u1", 0, "mt|chrM", 9)]
    with pytest.raises(ValueError, match="not contiguous"):
        list(group_hits(lines))


def h(ref, start, strand="+"):
    return Hit(ref, start, start + 20, strand, 0)


def test_single_focus_genome_is_assigned_with_fractional_weights():
    cls, genomes, placed = classify([h("mt|chrM", 0), h("mt|chrM", 50)], CONTIGS, ROLES,
                                    "report", "fractional")
    assert cls == "assigned" and genomes == ["mt"]
    assert [w for _, _, w in placed] == [0.5, 0.5]


def test_unique_multimap_policy_drops_multi_hits():
    cls, _, placed = classify([h("mt|chrM", 0), h("mt|chrM", 50)], CONTIGS, ROLES,
                              "report", "unique")
    assert cls == "multimap_dropped" and placed == []


def test_exclude_only_reads_are_excluded():
    cls, genomes, placed = classify([h("nuc|chr1", 5)], CONTIGS, ROLES, "report", "fractional")
    assert (cls, genomes, placed) == ("excluded", ["nuc"], [])


def test_focus_and_exclude_tie_is_ambiguous():
    hits = [h("mt|chrM", 0), h("nuc|chr1", 5)]
    cls, genomes, placed = classify(hits, CONTIGS, ROLES, "report", "fractional")
    assert cls == "ambiguous" and genomes == ["mt", "nuc"] and placed == []
    cls, _, placed = classify(hits, CONTIGS, ROLES, "split", "fractional")
    # the exclude genome keeps its share, so the focus hit gets half of the read
    assert cls == "ambiguous" and [(g, w) for g, _, w in placed] == [("mt", 0.5)]


def test_overhang_copy_is_folded_and_deduplicated():
    # the same read aligned at 5 and in the appended copy at 105
    cls, _, placed = classify([h("mt|chrM", 5), h("mt|chrM", 105)], CONTIGS, ROLES,
                              "report", "fractional")
    assert cls == "assigned"
    assert [(hit.start, w) for _, hit, w in placed] == [(5, 1.0)]
