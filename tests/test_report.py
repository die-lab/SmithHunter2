import gzip
import json

import pytest

from smithhunter.report import (coverage, env_pins, funnel, provenance, write_mqc_funnel,
                                yaml_lines)


def write(path, text):
    if str(path).endswith(".gz"):
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return str(path)


def tsv(*rows):
    return "".join("\t".join(map(str, r)) + "\n" for r in rows)


@pytest.fixture
def run(tmp_path):
    fastp = []
    for s, raw in (("s1", 100), ("s2", 50)):
        fastp.append(write(tmp_path / f"{s}.json", json.dumps({
            "summary": {"fastp_version": "0.23.4"},
            "read1_before_filtering": {"total_reads": raw}})))
    files = {
        "fastp": fastp,
        "libraries": write(tmp_path / "libraries.tsv", tsv(
            ("sample", "reads_in", "reads_kept"), ("s1", 90, 80), ("s2", 45, 40))),
        "summary": write(tmp_path / "summary.tsv", tsv(
            ("class", "genomes", "s1", "s2"), ("assigned", "mt", 50, 30),
            ("ambiguous", "mt+nuc", 10, 2), ("excluded", "nuc", 5, 3), ("unmapped", "", 15, 5))),
        "loci": write(tmp_path / "loci.tsv", tsv(
            ("locus_id", "genome", "locus_type", "pass", "count_s1", "count_s2"),
            ("mt_00001", "mt", "peak", 1, 30, 20), ("mt_00002", "mt", "peak", 0, 15.5, 6),
            ("mt_00003", "mt", "background", 0, 4.5, 4))),
    }
    return files


def test_funnel_counts_every_stage(run):
    rows = funnel(["s1", "s2"], run["fastp"], run["libraries"], run["summary"], run["loci"])
    got = {r["stage"]: (r["s1"], r["s2"], r["total"]) for r in rows}
    assert got == {
        "raw": (100, 50, 150), "trimmed": (90, 45, 135), "length_range": (80, 40, 120),
        "mapped": (65, 35, 100), "focus": (50, 30, 80), "in_peaks": (45.5, 26, 71.5),
        "in_passing": (30, 20, 50)}


def test_mqc_funnel_has_one_row_per_sample(run, tmp_path):
    rows = funnel(["s1", "s2"], run["fastp"], run["libraries"], run["summary"], run["loci"])
    path = tmp_path / "funnel_mqc.tsv"
    write_mqc_funnel(str(path), ["s1", "s2"], rows)
    lines = [l for l in path.read_text().splitlines() if not l.startswith("#")]
    assert lines[0].split("\t")[:2] == ["Sample", "Raw reads"]
    assert lines[2].split("\t")[:3] == ["s2", "50", "45"]


def test_coverage_wraps_circular_origin_and_skips_exclude_genomes(tmp_path):
    counts = write(tmp_path / "counts.tsv.gz", tsv(
        ("id", "sequence", "total", "s1", "s2"), ("r1", "A" * 20, 3, 2, 1),
        ("r2", "C" * 20, 1, 0, 1), ("r3", "G" * 20, 5, 5, 0), ("r4", "T" * 10, 1, 1, 0)))
    genomes = write(tmp_path / "genomes.tsv", tsv(
        ("name", "role", "min_length", "max_length"), ("mt", "focus", 18, 35),
        ("nuc", "exclude", 18, 35)))
    contigs = write(tmp_path / "contigs.tsv", tsv(
        ("ref", "genome", "contig", "length", "circular"), ("mt|chrM", "mt", "chrM", 100, 1),
        ("nuc|chr1", "nuc", "chr1", 1000, 0)))
    alignments = write(tmp_path / "alignments.tsv.gz", tsv(
        ("read_id", "genome", "contig", "start", "end", "strand", "weight", "mismatches"),
        ("r1", "mt", "chrM", 10, 30, "+", 1, 0),
        ("r2", "mt", "chrM", 90, 110, "-", 0.5, 0),     # across the origin
        ("r3", "nuc", "chr1", 10, 30, "+", 1, 0),       # exclude genome
        ("r4", "mt", "chrM", 10, 20, "+", 1, 0)))       # too short for mt
    rows = {(r[3], r[2]): r for r in coverage(alignments, counts, genomes, contigs)}
    plus = sorted(p for s, p in rows if s == "+")
    assert plus == list(range(11, 31))
    assert rows[("+", 11)][4:] == [3.0, 3.0, 2.0, 1.0]
    assert rows[("+", 12)][5] == 0.0
    minus = sorted(p for s, p in rows if s == "-")
    assert minus == list(range(1, 11)) + list(range(91, 101))
    # 5' end of a minus-strand read is its last base: position 110 -> 10 after wrapping
    assert rows[("-", 10)][5] == 0.5 and rows[("-", 91)][5] == 0.0
    assert rows[("-", 95)][4:] == [0.5, 0.0, 0.0, 0.5]


def test_env_pins(tmp_path):
    env = write(tmp_path / "env.yaml", "channels:\n  - conda-forge\n  - bioconda\n"
                "dependencies:\n  - bowtie=1.3.1\n  - python>=3.10  # comment\n  - pytest\n")
    assert env_pins([env]) == {"bowtie": "=1.3.1", "python": ">=3.10", "pytest": ""}


def test_yaml_lines_round_trip(run):
    yaml = pytest.importorskip("yaml")
    config = {"genomes": {"mito": {"fasta": "data/mt.fa", "role": "focus", "length": [18, 40]},
                          "nuclear": {"role": "exclude"}},
              "ambiguous": "report", "loci": {"min_samples": None, "peak_pvalue": 1e-05,
                                              "min_rpm": 5.0, "peaks": True},
              "trimming": {"adapter_r1": "", "extra": "say \"hi\": yes", "version": "1.20"},
              "yes": "no", "target_sets": {}, "list_of_dicts": [{"a": 1, "b": [1, 2]}, {"c": []}]}
    data = provenance(config, ["s1", "s2"], [], run["fastp"], run["loci"], {"smithhunter": "0.1"})
    assert yaml.safe_load("\n".join(yaml_lines(data))) == json.loads(json.dumps(data))
    assert data["results"] == {"mt": {"loci": 3, "peak_loci": 2, "passing": 1}}
    assert data["software"]["fastp (ran)"] == "0.23.4"
