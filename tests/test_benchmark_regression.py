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


# Chlamydia trachomatis in HeLa vs human THP-1 cells with M. tuberculosis (Furuse 2014);
# known sRNAs of Albrecht 2010 (benchmarks/wolbachia_multi/positive_ctr.tsv)
CTR = {name: ("NC_010285.1" if name == "pL2-sRNA1" else "NC_010280.1", s, e, st) for name, s, e, st in [
    ("ihtA", 53029, 53137, "+"), ("ctrR1", 53220, 53363, "+"), ("ctrR2", 126536, 126638, "-"),
    ("ctrR3", 342838, 342947, "-"), ("ctrR4", 497963, 498126, "-"), ("ctrR5", 586269, 586341, "+"),
    ("ctrR6", 725268, 725411, "-"), ("ctrR7", 876384, 876694, "-"), ("ctrR8", 938403, 938509, "+"),
    ("ctrR0332", 416302, 416641, "+"), ("pL2-sRNA1", 213, 292, "+")]}
CTR_PRECISE = {"ctrR1", "ctrR3", "ctrR4", "ctrR6", "ctrR7", "ctrR8", "ctrR0332", "pL2-sRNA1"}
SENSITIVE = yaml.safe_load((REPO / "config/profiles/symbiont_sensitive.yaml").read_text())


def found(rows) -> set[str]:
    return {name for name, where in CTR.items() if passing(rows, where)}


def test_chlamydia_known_srnas_with_the_profiles(tmp_path, monkeypatch):
    assert found(run("chlamydia_infected", tmp_path, monkeypatch, PRECISE)) == CTR_PRECISE
    assert found(run("chlamydia_infected", tmp_path, monkeypatch, SENSITIVE)) == set(CTR)


def test_chlamydia_min_rpm_genome_loses_known_srnas_at_high_load(tmp_path, monkeypatch):
    # why the profiles no longer set min_rpm_genome: C. trachomatis holds 23% of the library
    rows = run("chlamydia_infected", tmp_path, monkeypatch, PRECISE, {"loci": {"min_rpm_genome": 100}})
    assert found(rows) == {"ctrR1", "ctrR8", "ctrR0332", "pL2-sRNA1"}


def test_chlamydia_negative_control_rrna_reads_of_another_bacterium(tmp_path, monkeypatch):
    # M. tuberculosis rRNA reads land on the Chlamydia rRNA operons, on both strands
    rows = run("chlamydia_tet", tmp_path, monkeypatch)
    hits = [r for r in rows if r["pass"] == "1"]
    assert hits and all(r["annotation_class"] == "rRNA" for r in hits)
    sense = run("chlamydia_tet", tmp_path, monkeypatch, PRECISE, {"loci": {"exclude_structural": "sense"}})
    assert {r["structural"] for r in sense if r["pass"] == "1"} == {"antisense"}
    assert not [r for r in run("chlamydia_tet", tmp_path, monkeypatch, PRECISE) if r["pass"] == "1"]


# M. tuberculosis H37Rv in THP-1 vs C. trachomatis in HeLa (Furuse 2014); 19 RefSeq sRNAs
# (benchmarks/wolbachia_multi/positive_mtb.tsv)
def mtb_known() -> dict[str, tuple]:
    path = REPO / "benchmarks/wolbachia_multi/positive_mtb.tsv"
    with open(path, newline="") as fh:
        return {r["name"]: (r["contig"], int(r["start"]), int(r["end"]), r["strand"])
                for r in csv.DictReader(fh, delimiter="\t")}


MTB_PRECISE = {"B55", "mcr3", "mcr11", "mcr15", "G2", "MTS1338", "B11", "MTS2823"}


def test_mtb_known_srnas_with_the_profiles(tmp_path, monkeypatch):
    known = mtb_known()
    for profile in (PRECISE, SENSITIVE):
        rows = run("mtb_infected", tmp_path, monkeypatch, profile)
        assert {n for n, w in known.items() if passing(rows, w)} == MTB_PRECISE


def test_mtb_negative_control_false_positives_go_with_the_precise_profile(tmp_path, monkeypatch):
    # Chlamydia / HeLa reads on the Mtb rRNA operon (sense) and on conserved CDS (0% exact)
    hits = [r for r in run("mtb_tet", tmp_path, monkeypatch) if r["pass"] == "1"]
    assert {r["annotation_class"] for r in hits} >= {"rRNA", "CDS"}
    assert all(r["structural"] or float(r["exact_fraction"]) < 0.5 for r in hits)
    assert not [r for r in run("mtb_tet", tmp_path, monkeypatch, PRECISE) if r["pass"] == "1"]


def test_mtb_long_srnas_only_as_transcript_loci(tmp_path, monkeypatch):
    # mcr7 (350 nt) and F6 are split into weak peaks and background fragments: no candidate,
    # but an expressed transcript locus with loci.transcript_gap (report only)
    known = mtb_known()
    rows = run("mtb_infected", tmp_path, monkeypatch, PRECISE, {"loci": {"transcript_gap": 10}})
    for name in ("mcr7", "F6"):
        contig, start, end, strand = known[name]
        assert not passing(rows, known[name])
        tx = [r for r in rows if r["locus_type"] == "transcript" and r["strand"] == strand
              and int(r["start"]) <= end and int(r["end"]) >= start]
        assert any(r["expression_pass"] == "1" for r in tx) and not any(r["pass"] == "1" for r in tx)
