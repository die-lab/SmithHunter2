#!/usr/bin/env python3
"""One table over benchmarks and variants: candidates, false positives, structural share, positives.

    python3 06_summary.py BENCH [BENCH ...] [--positive positive.tsv] [--genome wolbachia]

For every BENCH it reads the main run (BENCH/{infected,tet}/results) and every variant made by
05_variants.py (BENCH/variants/<name>/{infected,tet}). Columns:
  infected / tet      passing loci of the focus genome (tet = false positives)
  structural          passing infected loci overlapping a structural RNA (column `structural`,
                      or the name regex of 04_compare.py for runs made before it existed)
  struct_reads        share of the infected passing reads in those loci
  known <name>        P = a passing locus overlaps the known small RNA on its strand,
                      - = loci there but none passing, 0 = no locus (only with --positive)
"""
from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("compare", HERE / "04_compare.py")
compare = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(compare)


def loci(results: Path, genome: str) -> list[dict] | None:
    p = results / "discovery/loci.tsv"
    if not p.exists():
        return None
    return [r for r in compare.read_tsv(p) if r["genome"] == genome]


def is_structural(r: dict) -> bool:
    if "structural" in r:
        return bool(r["structural"])
    return bool(compare.STRUCTURAL.search(r.get("annotation_types", "") + " " + r.get("annotation_names", "")))


def runs(bench: Path):
    yield "(main)", bench / "infected/results", bench / "tet/results"
    vdir = bench / "variants"
    if vdir.is_dir():
        for v in sorted(p for p in vdir.iterdir() if p.is_dir()):
            yield v.name, v / "infected", v / "tet"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bench", nargs="+", type=Path)
    ap.add_argument("--positive", type=Path)
    ap.add_argument("--genome", default="wolbachia")
    args = ap.parse_args()
    known = compare.read_tsv(args.positive) if args.positive else []
    header = ["bench", "variant", "infected", "tet", "structural", "struct_reads"] + [f"known {k['name']}" for k in known]
    rows = []
    for bench in args.bench:
        bench = bench.expanduser().resolve()
        for name, inf_dir, tet_dir in runs(bench):
            inf = loci(inf_dir, args.genome)
            if inf is None:
                continue
            tet = loci(tet_dir, args.genome)
            p = [r for r in inf if r["pass"] == "1"]
            s = [r for r in p if is_structural(r)]
            reads = sum(float(r["total_count"]) for r in p)
            sreads = sum(float(r["total_count"]) for r in s)
            row = [bench.name, name, len(p), "-" if tet is None else sum(r["pass"] == "1" for r in tet),
                   len(s), f"{sreads / reads:.0%}" if reads else "-"]
            for k in known:
                q = {"contig": k["contig"], "start": k["start"], "end": k["end"], "strand": k["strand"]}
                hits = [r for r in inf if r["contig"] == k["contig"] and compare.overlaps(r, q)]
                row.append("P" if any(r["pass"] == "1" for r in hits) else ("-" if hits else "0"))
            rows.append(row)
    print(compare.table(header, rows))


if __name__ == "__main__":
    main()
