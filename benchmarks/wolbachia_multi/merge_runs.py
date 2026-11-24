#!/usr/bin/env python3
"""Merge the infected and tet runs of a benchmark into one run and run `discover` on it.

The two module A runs of a study share their reference (focus, host and exclude genomes), so a
read sequence has the same alignments and weights in both. Reads are matched by sequence
(read ids are per run); only reads aligned to the focus genome are kept, which leaves the
focus loci unchanged. The samples table marks each library `infected` or `tet`, so
`loci.control_group: tet` and `loci.min_samples_scope: group` can be tried on real data.

    merge_runs.py --bench DIR --genome NAME --out DIR [--set loci.key=value ...]

Settings: the infected line's config, then --set. Output: OUT/{collapsed,origin,reference}
(merged inputs), OUT/config/samples.tsv, OUT/discovery/loci.tsv.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import subprocess
import sys
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
sys.path.insert(0, str(REPO / "src"))
from smithhunter.settings import KNOWN, discover_args  # noqa: E402


def results(line_dir: Path) -> tuple[dict, Path]:
    cfg = yaml.safe_load((line_dir / "config/config.yaml").read_text()) or {}
    return cfg, line_dir / cfg.get("outdir", "results")


def overrides(pairs: list[str]) -> dict:
    out: dict = {}
    for p in pairs:
        key, _, value = p.partition("=")
        section, _, name = key.partition(".")
        if (section, name) not in KNOWN:
            sys.exit(f"--set {key}: not a discover setting")
        out.setdefault(section, {})[name] = yaml.safe_load(value)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bench", required=True, type=Path)
    ap.add_argument("--genome", required=True, help="focus genome name")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--set", action="append", default=[], dest="sets")
    args = ap.parse_args()
    out = args.out.resolve()
    for sub in ("collapsed", "origin", "reference", "config", "discovery"):
        (out / sub).mkdir(parents=True, exist_ok=True)

    lines = {}
    for line in ("infected", "tet"):
        cfg, res = results(args.bench / line)
        lines[line] = (args.bench / line, cfg, res)

    # 1. focus-genome alignments of each run, by read id
    seq_id: dict[str, str] = {}          # sequence -> merged id
    aln_out: list[str] = []
    per_line_ids: dict[str, dict[str, str]] = {}
    for line, (_, _, res) in lines.items():
        focus_ids: set[str] = set()
        with gzip.open(res / "origin/alignments.tsv.gz", "rt") as fh:
            header = fh.readline()
            rows = [l for l in fh if l.split("\t", 3)[1] == args.genome]
        for l in rows:
            focus_ids.add(l.split("\t", 1)[0])
        ids: dict[str, str] = {}
        with gzip.open(res / "collapsed/counts.tsv.gz", "rt") as fh:
            fh.readline()
            for l in fh:
                rid, seq, _ = l.split("\t", 2)
                if rid in focus_ids:
                    ids[rid] = seq
        per_line_ids[line] = ids
        known = set(seq_id)               # sequences aligned from an earlier run
        for l in rows:
            rid, rest = l.split("\t", 1)
            seq = ids[rid]
            if seq in known:
                continue                  # same reference: same alignments already written
            mid = seq_id.setdefault(seq, f"m{len(seq_id) + 1:09d}")
            aln_out.append(f"{mid}\t{rest}")
        lines[line] += (header,)
        print(f"{line}: {len(focus_ids)} reads on {args.genome}, {len(rows)} alignments", flush=True)
    with gzip.open(out / "origin/alignments.tsv.gz", "wt") as fh:
        fh.write(lines["infected"][3])
        fh.writelines(aln_out)

    # 2. counts: one column per library of both runs
    samples: list[tuple[str, str]] = []
    counts: dict[str, list[float]] = {}
    for line, (_, _, res, _) in lines.items():
        with gzip.open(res / "collapsed/counts.tsv.gz", "rt") as fh:
            names = fh.readline().rstrip("\n").split("\t")[3:]
            offset = len(samples)
            samples += [(n, line) for n in names]
            ids = per_line_ids[line]
            for l in fh:
                f = l.rstrip("\n").split("\t")
                if f[0] not in ids:
                    continue
                row = counts.setdefault(ids[f[0]], [])
                row += [0.0] * (offset - len(row)) + [float(x) for x in f[3:]]
    n = len(samples)
    with gzip.open(out / "collapsed/counts.tsv.gz", "wt") as fh:
        fh.write("read_id\tsequence\tlength\t" + "\t".join(s for s, _ in samples) + "\n")
        for seq, mid in sorted(seq_id.items(), key=lambda kv: kv[1]):
            row = counts.get(seq, [])
            row += [0.0] * (n - len(row))
            fh.write(f"{mid}\t{seq}\t{len(seq)}\t" + "\t".join(f"{c:g}" for c in row) + "\n")

    # 3. libraries, reference (annotation paths made absolute), samples table
    with open(out / "collapsed/libraries.tsv", "w") as fh:
        first = True
        for line, (_, _, res, _) in lines.items():
            text = (res / "collapsed/libraries.tsv").read_text().splitlines()
            fh.write("\n".join(text if first else text[1:]) + "\n")
            first = False
    inf_dir, inf_cfg, inf_res, _ = lines["infected"]
    with open(inf_res / "reference/genomes.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        if r.get("annotation") and not os.path.isabs(r["annotation"]):
            r["annotation"] = str((inf_dir / r["annotation"]).resolve())
    with open(out / "reference/genomes.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    (out / "reference/contigs.tsv").write_text((inf_res / "reference/contigs.tsv").read_text())
    with open(out / "config/samples.tsv", "w") as fh:
        fh.write("sample\tfq1\tfq2\tgroup\n")
        fh.writelines(f"{s}\t-\t\t{g}\n" for s, g in samples)

    # 4. discover
    cfg = dict(inf_cfg)
    for section, values in overrides(args.sets).items():
        cfg[section] = dict(cfg.get(section) or {}, **values)
    (out / "config/settings.yaml").write_text(yaml.safe_dump((cfg.get("loci") or {}) | (cfg.get("ends") or {})))
    cmd = [sys.executable, "-m", "smithhunter", "discover",
           "--alignments", str(out / "origin/alignments.tsv.gz"),
           "--counts", str(out / "collapsed/counts.tsv.gz"),
           "--libraries", str(out / "collapsed/libraries.tsv"),
           "--genomes", str(out / "reference/genomes.tsv"),
           "--contigs", str(out / "reference/contigs.tsv")]
    cmd += discover_args(cfg, str(out / "config/samples.tsv"))
    cmd += ["--loci", str(out / "discovery/loci.tsv"),
            "--reads", str(out / "discovery/locus_reads.tsv.gz"),
            "--fasta", str(out / "discovery/candidates.fasta")]
    env = dict(os.environ, PYTHONPATH=str(REPO / "src") + os.pathsep + os.environ.get("PYTHONPATH", ""))
    print("discover:", " ".join(cmd[3:]), flush=True)
    subprocess.run(cmd, check=True, env=env, cwd=out)


if __name__ == "__main__":
    main()
