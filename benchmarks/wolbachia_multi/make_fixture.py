#!/usr/bin/env python3
"""Cut a small regression fixture out of a finished module A run (tests/data/<name>/).

    python3 make_fixture.py --line BENCH/jw18/infected --genome wolbachia \
        --region NC_002978.6:1204493-1204642 --region ... --out tests/data/jw18_infected

Kept: the alignments of the focus genome within --flank bases of each region, the reads they
use, the annotation features there, the library sizes, and the `group` column of samples.tsv.
`discover` divides by the reads of the whole genome (min_rpm_genome): the reads left out are
summed into one `filler` read on a separate contig `filler`, so rpm_genome stays exact. With
--check, discover runs on the fixture with the line's own settings and the loci overlapping the
regions must equal those of the line's loci.tsv (otherwise the flank is too short).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
from smithhunter.settings import KNOWN, discover_args  # noqa: E402

FILLER = "filler"
CHECK = ["start", "end", "strand", "locus_type", "total_count", "background_count",
         "samples_passing_rpm", "five_dominant_fraction", "isolation", "structural",
         "exact_fraction", "expression_pass", "ends_pass", "pass"]


def region(text: str):
    contig, _, span = text.rpartition(":")
    start, _, end = span.partition("-")
    return contig, int(start), int(end)


def read_tsv(path: Path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], cols: list[str]):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--line", required=True, type=Path)
    ap.add_argument("--genome", default="wolbachia")
    ap.add_argument("--region", action="append", required=True, type=region,
                    help="contig:start-end, 1-based inclusive (repeatable)")
    ap.add_argument("--flank", type=int, default=400)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    line = args.line.resolve()
    cfg = yaml.safe_load((line / "config/config.yaml").read_text()) or {}
    res = line / cfg.get("outdir", "results")
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)

    genome = next(g for g in read_tsv(res / "reference/genomes.tsv") if g["name"] == args.genome)
    lo, hi = int(genome["min_length"]), int(genome["max_length"])
    windows = defaultdict(list)
    for contig, s, e in args.region:
        windows[contig].append((s - 1 - args.flank, e + args.flank))

    def inside(contig, start, end):
        return any(start < we and end > ws for ws, we in windows.get(contig, ()))

    # 1. alignments of the genome: keep those in the windows, weigh every read for the totals
    weight_all: dict[str, float] = defaultdict(float)
    kept_aln, kept_weight = [], defaultdict(float)
    with gzip.open(res / "origin/alignments.tsv.gz", "rt") as fh:
        header = fh.readline()
        for line_ in fh:
            f = line_.rstrip("\n").split("\t")
            if f[1] != args.genome:
                continue
            weight_all[f[0]] += float(f[6])
            if inside(f[2], int(f[3]), int(f[4])):
                kept_aln.append(line_)
                kept_weight[f[0]] += float(f[6])

    # 2. counts: kept reads verbatim, everything else of the genome into the filler read
    with gzip.open(res / "collapsed/counts.tsv.gz", "rt") as fh:
        cheader = fh.readline()
        samples = cheader.rstrip("\n").split("\t")[3:]
        residual = [0.0] * len(samples)
        kept_counts = []
        for line_ in fh:
            rid, seq, rest = line_.split("\t", 2)
            w = weight_all.get(rid)
            if not w or not lo <= len(seq) <= hi:
                continue
            row = [float(x) for x in rest.rstrip("\n").split("\t")[1:]]
            if rid in kept_weight:
                kept_counts.append(line_)
            left = w - kept_weight.get(rid, 0.0)
            for i, c in enumerate(row):
                residual[i] += left * c
    filler_seq = "A" * lo
    with gzip.open(out / "counts.tsv.gz", "wt") as fh:
        fh.write(cheader)
        fh.writelines(kept_counts)
        fh.write(f"{FILLER}\t{filler_seq}\t{lo}\t" + "\t".join(f"{x:.6f}" for x in residual) + "\n")
    with gzip.open(out / "alignments.tsv.gz", "wt") as fh:
        fh.write(header)
        fh.writelines(kept_aln)
        fh.write(f"{FILLER}\t{args.genome}\t{FILLER}\t0\t{lo}\t+\t1\t0\n")

    # 3. reference tables, annotation, libraries, samples
    contigs = [c for c in read_tsv(res / "reference/contigs.tsv") if c["genome"] == args.genome]
    contigs.append({"ref": f"{args.genome}|{FILLER}", "genome": args.genome, "contig": FILLER,
                    "length": str(lo), "circular": "0"})
    write_tsv(out / "contigs.tsv", contigs, ["ref", "genome", "contig", "length", "circular"])
    ann = genome.get("annotation") or ""
    if ann:
        src = Path(ann) if Path(ann).is_absolute() else line / ann
        with open(src) as fh, open(out / "annotation.gff3", "w") as dst:
            for gl in fh:
                if gl.startswith("#"):
                    if gl.startswith("##gff-version"):
                        dst.write(gl)
                    continue
                f = gl.split("\t")
                if len(f) >= 5 and inside(f[0], int(f[3]) - 1, int(f[4])):
                    dst.write(gl)
    gcols = list(genome)
    write_tsv(out / "genomes.tsv", [dict(genome, fasta="-", annotation="annotation.gff3" if ann else "")], gcols)
    (out / "libraries.tsv").write_text((res / "collapsed/libraries.tsv").read_text())
    srows = read_tsv(line / cfg.get("samples", "config/samples.tsv"))
    write_tsv(out / "samples.tsv", srows, ["sample"] + (["group"] if "group" in srows[0] else []))
    # the line's discover settings: tests add a profile on top, as 05_variants.py does
    own = {s: {k: v for k, v in (cfg.get(s) or {}).items() if (s, k) in KNOWN}
           for s in ("loci", "ends")}
    (out / "settings.yaml").write_text(yaml.safe_dump(own, sort_keys=True))
    size = sum(p.stat().st_size for p in out.iterdir())
    print(f"{out}: {len(kept_aln)} alignments, {len(kept_counts)} reads, {size / 1e3:.0f} kB")

    if args.check:
        check(out, cfg, res, args)


def same(a: str, b: str) -> bool:
    """Equal, or equal numbers up to the rounding of loci.tsv (sums in another order)."""
    if a == b:
        return True
    try:
        return abs(float(a) - float(b)) <= max(0.002, 1e-6 * abs(float(a)))
    except ValueError:
        return False


def check(out: Path, cfg: dict, res: Path, args) -> None:
    loci = out / "check_loci.tsv"
    cmd = [sys.executable, "-m", "smithhunter", "discover",
           "--alignments", str(out / "alignments.tsv.gz"), "--counts", str(out / "counts.tsv.gz"),
           "--libraries", str(out / "libraries.tsv"), "--genomes", "genomes.tsv",
           "--contigs", str(out / "contigs.tsv"), *discover_args(cfg, str(out / "samples.tsv")),
           "--loci", str(loci), "--reads", str(out / "check_reads.tsv.gz"),
           "--fasta", str(out / "check.fasta")]
    env = {"PYTHONPATH": str(REPO / "src"), "PATH": __import__("os").environ["PATH"]}
    subprocess.run(cmd, check=True, cwd=out, env=env, stdout=subprocess.DEVNULL)

    def near(rows):
        keep = {}
        for r in rows:
            for contig, s, e in args.region:
                if r["contig"] == contig and int(r["start"]) <= e and int(r["end"]) >= s:
                    keep[(r["start"], r["end"], r["strand"], r["locus_type"])] = r
        return keep

    full = near(read_tsv(res / "discovery/loci.tsv"))
    small = near(read_tsv(loci))
    bad = []
    for k in sorted(set(full) | set(small)):
        a, b = full.get(k), small.get(k)
        if not a or not b:
            # background loci never pass; their chains can be cut at the edge of a window
            if k[3] != "background":
                bad.append(f"{k}: only in {'fixture' if b else 'full run'}")
            continue
        for col in CHECK + [c for c in a if c.startswith("rpm_genome_")]:
            if col in a and col in b and not same(a[col], b[col]):
                bad.append(f"{k} {col}: full {a[col]} fixture {b[col]}")
    for p in (loci, out / "check_reads.tsv.gz", out / "check.fasta"):
        p.unlink()
    if bad:
        sys.exit("fixture differs from the full run:\n  " + "\n  ".join(bad[:30]))
    print(f"check: {len(full)} loci overlapping the regions identical to the full run")


if __name__ == "__main__":
    main()
