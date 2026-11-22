"""Regression tests on real data: small cuts of the Wolbachia benchmarks (tests/data/).

Each fixture holds the reads around the known small RNAs and the known false positives of one
line (infected or tetracycline-cured), made by benchmarks/wolbachia_multi/make_fixture.py and
checked there against the full run (same loci, counts and filters). The expectations are the
results of docs/module_a_benchmarks.md: a change that breaks them changes the benchmark.
"""
import csv
from pathlib import Path

import pytest
import yaml

from smithhunter import cli
from smithhunter.settings import discover_args

REPO = Path(__file__).resolve().parent.parent
DATA = Path(__file__).resolve().parent / "data"
PRECISE = yaml.safe_load((REPO / "config/profiles/symbiont_precise.yaml").read_text())

WSN46 = ("NC_002978.6", 1204493, 1204642, "+")      # wMel, Rfam RF02625 precursor
WSN59 = ("NC_002978.6", 582143, 582231, "-")        # RF02626
JW18_CROSSTALK = ("NC_002978.6", 229440, 229489, "-")
WSN7542 = ("NZ_CP099592.1", 489055, 489072, "-")    # spider mite wTtru
WSN744 = ("NZ_CP099592.1", 248655, 248679, "+")     # male-biased
MITE_TRNA = ("NZ_CP099592.1", 549825, 549860, "-")  # another bacterium's tRNA reads
MITE_DNAA = ("NZ_CP099592.1", 1395447, 1395468, "-")  # dnaA antisense, 0% exact reads


def merged(*layers: dict) -> dict:
    out: dict = {}
    for layer in layers:
        for section, values in (layer or {}).items():
            out[section] = {**out.get(section, {}), **(values or {})}
    return out


def run(name: str, tmp_path, monkeypatch, *layers: dict) -> list[dict]:
    d = DATA / name
    if not d.exists():
        pytest.skip(f"fixture {name} missing")
    cfg = merged(yaml.safe_load((d / "settings.yaml").read_text()), *layers)
    monkeypatch.chdir(d)          # genomes.tsv names the annotation relative to the fixture
    loci = tmp_path / f"{name}.tsv"
    cli.main(["discover", "--alignments", "alignments.tsv.gz", "--counts", "counts.tsv.gz",
              "--libraries", "libraries.tsv", "--genomes", "genomes.tsv", "--contigs", "contigs.tsv",
              *discover_args(cfg, "samples.tsv"), "--loci", str(loci),
              "--reads", str(tmp_path / "reads.tsv.gz"), "--fasta", str(tmp_path / "c.fasta")])
    with open(loci, newline="") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t") if r["contig"] != "filler"]


def passing(rows, where) -> list[dict]:
    contig, start, end, strand = where
    return [r for r in rows if r["pass"] == "1" and r["contig"] == contig and r["strand"] == strand
            and int(r["start"]) <= end and int(r["end"]) >= start]


def test_jw18_known_small_rnas_pass_with_the_precise_profile(tmp_path, monkeypatch):
    rows = run("jw18_infected", tmp_path, monkeypatch, PRECISE)
    assert passing(rows, WSN46) and passing(rows, WSN59)


def test_jw18_cured_line_has_no_candidate_with_the_precise_profile(tmp_path, monkeypatch):
    rows = run("jw18_tet", tmp_path, monkeypatch, PRECISE)
    assert not [r for r in rows if r["pass"] == "1"]


def test_jw18_cured_line_cross_talk_passes_at_one_rpm(tmp_path, monkeypatch):
    # known limitation (finding 10): the most abundant symbiont small RNAs reach the cured
    # libraries; only the library threshold keeps them out
    rows = run("jw18_tet", tmp_path, monkeypatch, PRECISE, {"loci": {"min_rpm": 1}})
    assert passing(rows, JW18_CROSSTALK)


def test_mite_cured_line_false_positives_at_the_defaults(tmp_path, monkeypatch):
    rows = run("mite_tet", tmp_path, monkeypatch)
    trna, dnaa = passing(rows, MITE_TRNA), passing(rows, MITE_DNAA)
    assert trna and dnaa
    assert trna[0]["structural"] == "sense"
    assert float(dnaa[0]["exact_fraction"]) < 0.5


def test_mite_cured_line_is_clean_with_the_precise_profile(tmp_path, monkeypatch):
    rows = run("mite_tet", tmp_path, monkeypatch, PRECISE)
    assert not [r for r in rows if r["pass"] == "1"]


def test_mite_male_biased_small_rna_needs_the_group_scope(tmp_path, monkeypatch):
    rows = run("mite_infected", tmp_path, monkeypatch, PRECISE)
    assert passing(rows, WSN7542) and not passing(rows, WSN744)
    rows = run("mite_infected", tmp_path, monkeypatch, PRECISE,
               {"loci": {"min_samples_scope": "group"}})
    assert passing(rows, WSN7542) and passing(rows, WSN744)


def test_mite_cured_line_stays_clean_with_the_group_scope(tmp_path, monkeypatch):
    rows = run("mite_tet", tmp_path, monkeypatch, PRECISE, {"loci": {"min_samples_scope": "group"}})
    assert not [r for r in rows if r["pass"] == "1"]
