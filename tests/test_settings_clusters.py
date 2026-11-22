"""One list of discover settings (smithhunter.settings) and the cluster columns of loci.tsv."""
import argparse
import re
from pathlib import Path

import pytest
import yaml

from smithhunter import cli
from smithhunter.loci import Locus, peak_clusters
from smithhunter.settings import DISCOVER_SETTINGS, KNOWN, discover_args, effective

REPO = Path(__file__).resolve().parent.parent


def discover_options() -> dict[str, argparse.Action]:
    sub = next(a for a in cli.build_parser()._actions
               if isinstance(a, argparse._SubParsersAction))
    parser = sub.choices["discover"]
    return {s: a for a in parser._actions for s in a.option_strings if s.startswith("--")}


def test_every_setting_is_a_discover_flag_with_the_same_default():
    options = discover_options()
    for flag, section, key, default in DISCOVER_SETTINGS:
        assert flag in options, flag
        assert options[flag].default == default, flag


def test_every_discover_option_is_reachable_from_the_config():
    io = {"--alignments", "--counts", "--libraries", "--genomes", "--contigs", "--loci",
          "--reads", "--fasta", "--help",
          "--no-peaks", "--min-samples-scope", "--samples-table"}   # loci.peaks, ..._scope
    assert set(discover_options()) - io == {f for f, *_ in DISCOVER_SETTINGS}


def test_main_config_defaults_match_the_settings():
    cfg = yaml.safe_load((REPO / "config/config.yaml").read_text())
    for flag, section, key, default in DISCOVER_SETTINGS:
        if key in (cfg.get(section) or {}):
            got = cfg[section][key]
            if key == "min_count":   # the workflow default (5) is stricter than the CLI one (1)
                continue
            assert got == default, f"{section}.{key}: config {got!r}, settings {default!r}"


def test_discover_args_from_config():
    args = discover_args({"loci": {"min_rpm": 2, "min_samples": 3.0, "peaks": False,
                                   "exclude_structural": "sense"}})
    joined = " ".join(args)
    assert "--min-rpm 2" in joined and "--min-samples 3" in joined
    assert "--exclude-structural sense" in joined and "--no-peaks" in args
    assert "--min-samples-scope" not in args
    assert "--min-samples" not in discover_args({})       # null: n - 1 inside discover


def test_discover_args_group_scope_needs_the_table():
    cfg = {"loci": {"min_samples_scope": "group"}}
    with pytest.raises(ValueError):
        discover_args(cfg)
    args = discover_args(cfg, "config/samples.tsv")
    assert args[-4:] == ["--min-samples-scope", "group", "--samples-table", "config/samples.tsv"]
    with pytest.raises(ValueError):
        discover_args({"loci": {"min_samples_scope": "sex"}}, "x")


def test_profiles_use_known_keys_only():
    for path in (REPO / "config/profiles").glob("*.yaml"):
        cfg = yaml.safe_load(path.read_text())
        for section in ("loci", "ends"):
            for key in (cfg.get(section) or {}):
                assert (section, key) in KNOWN, f"{path.name}: {section}.{key}"
        assert effective(cfg)["loci.min_rpm"] == cfg["loci"]["min_rpm"]


def test_variants_script_has_no_copy_of_the_settings():
    text = (REPO / "benchmarks/wolbachia_walbb/05_variants.py").read_text()
    assert not re.search(r'\("--min-rpm"', text)
    smk = (REPO / "workflow/rules/discovery.smk").read_text()
    assert "--min-rpm" not in smk


def _peak(start, end, strand="+", contig="c", kind="peak"):
    return Locus("g", contig, strand, start, end, kind=kind)


def test_peak_clusters_chain_overlapping_peaks_only():
    loci = [_peak(0, 30), _peak(20, 50), _peak(49, 70), _peak(80, 100), _peak(10, 40, "-"),
            _peak(25, 45, kind="background"), _peak(0, 30, contig="d")]
    c = peak_clusters(loci)
    assert c[0] == c[1] == c[2]                  # 0-30, 20-50, 49-70 chain
    assert c[3] != c[0]                          # 80-100 is 10 bases away
    assert c[4] not in (c[0], c[3])              # other strand
    assert c[5] is None                          # background loci are never clustered
    assert c[6] not in (c[0], c[3], c[4])        # other contig
    assert peak_clusters(loci, gap=10)[3] == peak_clusters(loci, gap=10)[0]


def test_cluster_columns_in_discover(tmp_path):
    from smithhunter.loci import discover
    d = tmp_path
    (d / "genomes.tsv").write_text("name\tfasta\trole\tcircular\tannotation\tmin_length\tmax_length\n"
                                   "g\tg.fa\tfocus\t0\t\t18\t35\n")
    (d / "contigs.tsv").write_text("ref\tgenome\tcontig\tlength\tcircular\ng|c\tg\tc\t10000\t0\n")
    (d / "libraries.tsv").write_text("sample\treads_in\treads_kept\ns1\t1000\t1000\ns2\t1000\t1000\n")
    # two 5' isoforms 4 nt apart (each above the satellite fraction) and one distant peak
    counts = ["read_id\tsequence\tlength\ts1\ts2"]
    aln = ["read_id\tgenome\tcontig\tstart\tend\tstrand\tweight\tmismatches"]
    for rid, start, n in (("a", 1000, 60), ("b", 1004, 40), ("c", 3000, 50)):
        counts.append(f"{rid}\t{'A' * 22}\t22\t{n}\t{n}")
        aln.append(f"{rid}\tg\tc\t{start}\t{start + 22}\t+\t1\t0")
    import gzip
    with gzip.open(d / "counts.tsv.gz", "wt") as fh:
        fh.write("\n".join(counts) + "\n")
    with gzip.open(d / "aln.tsv.gz", "wt") as fh:
        fh.write("\n".join(aln) + "\n")
    _, rows, _ = discover(str(d / "aln.tsv.gz"), str(d / "counts.tsv.gz"), str(d / "libraries.tsv"),
                          str(d / "genomes.tsv"), str(d / "contigs.tsv"), min_rpm=1, min_samples=1)
    peaks = [r for r in rows if r["locus_type"] == "peak"]
    assert len(peaks) == 3
    a, b, c = sorted(peaks, key=lambda r: r["start"])
    assert a["cluster_id"] == b["cluster_id"] == a["locus_id"]   # the stronger isoform leads
    assert a["cluster_size"] == b["cluster_size"] == 2 and a["cluster_count"] == 200
    assert c["cluster_id"] == c["locus_id"] and c["cluster_size"] == 1
