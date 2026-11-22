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
sys.path.insert(0, str(REPO / "src"))
# the same list of settings as workflow/rules/discovery.smk
from smithhunter.settings import KNOWN, discover_args, effective  # noqa: E402


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


def with_overrides(cfg: dict, overrides: dict) -> dict:
    out = {k: dict(v) if isinstance(v, dict) else v for k, v in cfg.items()}
    for (section, key), value in overrides.items():
        out.setdefault(section, {})
        out[section] = dict(out[section] or {}, **{key: value})
    return out


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
    merged = with_overrides(cfg, overrides)
    used = {k: v for k, v in effective(merged).items() if v is not None}
    cmd += discover_args(merged, str(line_dir / cfg.get("samples", "config/samples.tsv")))
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
