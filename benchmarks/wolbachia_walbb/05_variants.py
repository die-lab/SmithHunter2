#!/usr/bin/env python3
"""Re-run only `discover` with other settings on a finished benchmark, without touching it.

    python3 05_variants.py --bench DIR --name sense --set loci.exclude_structural=sense \
        [--set loci.min_rpm=10 ...] [--lines infected tet] [--genome wolbachia]

For each line (DIR/<line>, a finished module A run with config/config.yaml and results/), the
settings are those of the line's config with the --set overrides; outputs go to
DIR/variants/<name>/<line>/ (discovery/ plus links to the inputs 04_compare.py reads), then
04_compare.py writes DIR/variants/<name>/compare_report.md. Cheap: no alignment, minutes.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent

# discover flag <- (config section, key, default): mirrors workflow/rules/discovery.smk
PARAMS = [
    ("--merge-gap", "loci", "merge_gap", 0), ("--min-rpm", "loci", "min_rpm", 5),
    ("--min-count", "loci", "min_count", 1), ("--min-samples", "loci", "min_samples", None),
    ("--min-rpm-genome", "loci", "min_rpm_genome", 0),
    ("--peak-window", "loci", "peak_window", 2), ("--peak-pvalue", "loci", "peak_pvalue", 0.001),
    ("--background-flank", "loci", "background_flank", 50),
    ("--satellite-distance", "loci", "satellite_distance", 10),
    ("--satellite-fraction", "loci", "satellite_fraction", 0.1),
    ("--exclude-structural", "loci", "exclude_structural", "none"),
    ("--n-thre", "ends", "n_thre", 0.5), ("--penalty", "ends", "penalty", 0.1),
    ("--min-five", "ends", "min_five_score", 0.5), ("--min-three", "ends", "min_three_score", 0.0),
    ("--min-isolation", "ends", "min_isolation", 0.5),
    ("--five-filter", "ends", "five_filter", "dominant"),
    ("--min-five-fraction", "ends", "min_five_fraction", 0.5),
    ("--five-jitter", "ends", "five_jitter", 1),
]


KNOWN = {(section, key) for _, section, key, _ in PARAMS} | {("loci", "peaks")}


def parse_set(items: list[str]) -> dict[tuple[str, str], object]:
    """section.key=value; only keys that `discover` uses (upstream steps are not re-run)."""
    out = {}
    for item in items:
        key, eq, value = item.partition("=")
        section, dot, name = key.partition(".")
        if not (eq and dot and name):
            sys.exit(f"--set needs section.key=value, got {item!r}")
        if (section, name) not in KNOWN:
            sys.exit(f"--set {key}: not a discover setting (known: "
                     f"{', '.join(sorted(f'{s}.{k}' for s, k in KNOWN))})")
        out[(section, name)] = yaml.safe_load(value)
    return out


def setting(cfg: dict, overrides: dict, section: str, key: str, default):
    return overrides.get((section, key), (cfg.get(section) or {}).get(key, default))


def link(src: Path, dst: Path) -> None:
    if not dst.exists():
        dst.symlink_to(src)


def run_line(line_dir: Path, out: Path, overrides: dict):
    """Run discover for one line; returns its effective loci.min_samples (None = n - 1)."""
    cfg = yaml.safe_load((line_dir / "config/config.yaml").read_text()) or {}
    res = line_dir / cfg.get("outdir", "results")
    if not (res / "origin/alignments.tsv.gz").exists():
        sys.exit(f"{res}: not a finished module A run")
    for sub in ("collapsed", "origin", "reference"):
        out.mkdir(parents=True, exist_ok=True)
        link(res / sub, out / sub)
    (out / "discovery").mkdir(exist_ok=True)
    cmd = [sys.executable, "-m", "smithhunter", "discover",
           "--alignments", str(res / "origin/alignments.tsv.gz"),
           "--counts", str(res / "collapsed/counts.tsv.gz"),
           "--libraries", str(res / "collapsed/libraries.tsv"),
           "--genomes", str(res / "reference/genomes.tsv"),
           "--contigs", str(res / "reference/contigs.tsv")]
    used = {}
    for flag, section, key, default in PARAMS:
        value = setting(cfg, overrides, section, key, default)
        if value is None:
            continue
        cmd += [flag, str(value)]
        used[f"{section}.{key}"] = value
    peaks = setting(cfg, overrides, "loci", "peaks", True)
    used["loci.peaks"] = peaks
    if not peaks:
        cmd.append("--no-peaks")
    cmd += ["--loci", str(out / "discovery/loci.tsv"),
            "--reads", str(out / "discovery/locus_reads.tsv.gz"),
            "--fasta", str(out / "discovery/candidates.fasta")]
    (out / "settings.yaml").write_text(yaml.safe_dump(used, sort_keys=True))
    env = dict(os.environ, PYTHONPATH=str(REPO / "src") + os.pathsep + os.environ.get("PYTHONPATH", ""))
    print(f"== {line_dir.name}: discover", flush=True)
    # genomes.tsv keeps annotation paths as written in the config (relative to the run dir)
    subprocess.run(cmd, check=True, env=env, cwd=line_dir)
    return used.get("loci.min_samples")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bench", required=True, type=Path)
    ap.add_argument("--name", required=True, help="variant name (output folder)")
    ap.add_argument("--set", action="append", default=[], help="section.key=value (repeatable)")
    ap.add_argument("--lines", nargs="+", default=["infected", "tet"])
    ap.add_argument("--genome", default="wolbachia")
    ap.add_argument("--genome-fasta", type=Path)
    ap.add_argument("--positive", type=Path, help="known small RNAs, see 04_compare.py --positive")
    args = ap.parse_args()
    bench = args.bench.expanduser().resolve()
    overrides = parse_set(args.set)
    vdir = bench / "variants" / args.name
    min_samples = {line: run_line(bench / line, vdir / line, overrides) for line in args.lines}
    if set(args.lines) >= {"infected", "tet"}:
        cmd = [sys.executable, str(HERE / "04_compare.py"), "--infected", str(vdir / "infected"),
               "--tet", str(vdir / "tet"), "--genome", args.genome]
        if args.genome_fasta:
            cmd += ["--genome-fasta", str(args.genome_fasta)]
        if args.positive:
            cmd += ["--positive", str(args.positive)]
        ms = min_samples.get("infected")
        if ms is not None:
            cmd += ["--min-samples", str(ms)]
        report = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
        (vdir / "compare_report.md").write_text(report)
        print(f"report: {vdir / 'compare_report.md'}")


if __name__ == "__main__":
    main()
