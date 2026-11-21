#!/usr/bin/env python3
"""Where do the SmithHunter2 loci of one focus genome come from? Markdown report on stdout.

    python3 04_compare.py --infected infected/results [--tet tet/results] --genome wolbachia \
        [--genome-fasta data/wAlbB.fasta --v0-fasta v0_smithRNAs.fasta]

Sections: reads assigned to the genome; loci by type and filter; passing loci by annotation
(structural RNAs flagged); how many loci would pass at stricter min_rpm / min_count (approximate:
gross per-sample counts, the background is not subtracted); the negative control (tet: no
Wolbachia, every passing locus is a false positive); optional overlap with the v0 candidates
(sequences located on the genome by exact match). Standard library only.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from collections import Counter
from pathlib import Path

STRUCTURAL = re.compile(r"rRNA|tRNA|tmRNA|ssrA|ssrS|rnpB|ffs|\b(16S|23S|5S|6S)\b|RNase.?P|SRP|signal recognition",
                        re.IGNORECASE)
RPM_STEPS = (5, 10, 20, 50, 100)
COUNT_STEPS = (5, 10, 20)


def read_tsv(path: Path) -> list[dict]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_fasta(path: Path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:], []
            elif line:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def revcomp(s: str) -> str:
    return s[::-1].translate(str.maketrans("ACGTUN", "TGCAAN"))


class Run:
    def __init__(self, results: Path, genome: str):
        self.dir = results
        all_loci = read_tsv(results / "discovery/loci.tsv")
        self.loci = [r for r in all_loci if r["genome"] == genome]
        self.samples = [c[len("count_"):] for c in (all_loci[0] if all_loci else {}) if c.startswith("count_")]
        if not self.samples:
            self.samples = [r["sample"] for r in read_tsv(results / "collapsed/libraries.tsv")]
        self.library = {r["sample"]: float(r["reads_kept"]) for r in read_tsv(results / "collapsed/libraries.tsv")}
        self.origin = read_tsv(results / "origin/origin_summary.tsv")
        self.genome = genome

    def passing(self):
        return [r for r in self.loci if r["pass"] == "1"]

    def reads_on_genome(self) -> list[tuple[str, str, float]]:
        out = []
        for r in self.origin:
            if self.genome in re.split(r"[+,]", r["genomes"]):
                out.append((r["class"], r["genomes"], sum(float(r.get(s, 0) or 0) for s in self.samples)))
        return out


def fmt(x: float) -> str:
    return f"{x:,.0f}" if abs(x) >= 100 else f"{x:.3g}"


def table(header: list[str], rows: list[list]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return "\n".join(lines)


def section_reads(run: Run, label: str) -> str:
    total = sum(run.library.get(s, 0) for s in run.samples)
    rows = [[c, g, fmt(n), f"{n / total:.3%}" if total else "-"] for c, g, n in run.reads_on_genome()]
    return (f"### Reads touching `{run.genome}` ({label}; {len(run.samples)} libraries, "
            f"{fmt(total)} reads kept after trimming)\n\n"
            + (table(["class", "genomes", "reads", "share of library"], rows) if rows else "none") + "\n")


def section_loci(run: Run, label: str) -> str:
    loci = run.loci
    types = Counter(r["locus_type"] for r in loci)
    peak = [r for r in loci if r["locus_type"] == "peak"]
    expr = sum(r["expression_pass"] == "1" for r in peak)
    ends = sum(r["ends_pass"] == "1" for r in peak)
    both = len(run.passing())
    rows = [["loci in loci.tsv (all types)", len(loci)]]
    rows += [[f"  locus_type = {t}", n] for t, n in sorted(types.items())]
    rows += [["peak loci passing expression", expr], ["peak loci passing ends + isolation", ends],
             ["**pass = 1 (candidates.fasta)**", f"**{both}**"]]
    return f"### Loci of `{run.genome}` ({label})\n\n" + table(["", "n"], rows) + "\n"


def section_annotation(run: Run, label: str) -> str:
    p = run.passing()
    if not p:
        return f"### Passing loci by annotation ({label})\n\nnone\n"
    cls = Counter((r.get("annotation_class", "NA"), r.get("annotation_orientation", "")) for r in p)
    if "structural" in p[0]:   # column written by `discover` since the structural filter
        structural = [r for r in p if r["structural"]]
        sense = sum(r["structural"] == "sense" for r in structural)
        split = f" ({sense} sense, {len(structural) - sense} antisense)"
    else:
        structural = [r for r in p if STRUCTURAL.search(r.get("annotation_types", "") + " " + r.get("annotation_names", ""))]
        split = ""
    rows = [[c, o or "-", n] for (c, o), n in cls.most_common()]
    reads_all = sum(float(r["total_count"]) for r in p)
    reads_struct = sum(float(r["total_count"]) for r in structural)
    top = sorted(p, key=lambda r: -float(r["total_count"]))[:15]
    top_rows = [[r["locus_id"], f"{r['start']}-{r['end']}{r['strand']}", r["rep_length"], fmt(float(r["total_count"])),
                 r["samples_passing_rpm"], r.get("annotation_class", ""),
                 (r.get("annotation_names", "") or "")[:40], r["rep_sequence"]] for r in top]
    return (f"### Passing loci by annotation ({label})\n\n" + table(["class", "orientation", "loci"], rows)
            + f"\n\nOverlapping structural RNAs (rRNA, tRNA, tmRNA, RNase P, SRP, 6S): **{len(structural)} of {len(p)}** "
              f"loci{split}, {reads_struct / reads_all:.0%} of their reads.\n\nTop 15 by reads:\n\n"
            + table(["locus", "position", "len", "reads", "samples", "class", "names", "rep_sequence"], top_rows) + "\n")


def section_sweep(run: Run, label: str, min_samples: int | None) -> str:
    n = len(run.samples)
    need = min_samples if min_samples is not None else max(1, n - 1)
    cand = [r for r in run.loci if r["locus_type"] == "peak" and r["ends_pass"] == "1"]
    rows = []
    for mc in COUNT_STEPS:
        row = [mc]
        for rpm in RPM_STEPS:
            k = 0
            for r in cand:
                ok = sum(float(r[f"count_{s}"]) >= mc and float(r[f"rpm_{s}"]) >= rpm for s in run.samples)
                k += ok >= need
            row.append(k)
        rows.append(row)
    lib = sorted(run.library.get(s, 0) for s in run.samples)
    return (f"### Stricter expression filter ({label}; loci passing ends, in >= {need} of {n} libraries)\n\n"
            f"Approximate: gross counts, the local background is not subtracted, so the "
            f"(5 RPM, 5 reads) cell can be a little above the real pass count. Library sizes "
            f"{fmt(lib[0])}-{fmt(lib[-1])} reads: 5 RPM = {fmt(5e-6 * lib[0])}-{fmt(5e-6 * lib[-1])} reads.\n\n"
            + table(["min_count \\ min_rpm"] + [str(x) for x in RPM_STEPS], rows) + "\n")


def intervals(r: dict) -> list[tuple[int, int]]:
    """1-based closed intervals; a locus across the origin of a circular contig has end < start."""
    s, e = int(r["start"]), int(r["end"])
    return [(s, e)] if s <= e else [(s, 10 ** 12), (1, e)]


def overlaps(a: dict, b: dict) -> bool:
    if a["contig"] != b["contig"] or a["strand"] != b["strand"]:
        return False
    return any(s1 <= e2 and s2 <= e1 for s1, e1 in intervals(a) for s2, e2 in intervals(b))


def section_tet(inf: Run, tet: Run) -> str:
    p = tet.passing()
    shared = [r for r in p if any(overlaps(r, q) for q in inf.passing())]
    txt = (f"## Negative control (tet: no Wolbachia)\n\nEvery `{tet.genome}` locus passing in the tet line "
           f"is a false positive: host reads (or contamination) placed on the symbiont.\n\n"
           f"- passing loci in tet: **{len(p)}** (infected: {len(inf.passing())})\n"
           f"- of these, overlapping an infected passing locus: {len(shared)}\n\n")
    return txt + section_reads(tet, "tet") + "\n" + section_loci(tet, "tet") + "\n" + section_annotation(tet, "tet")


def section_v0(inf: Run, v0_fasta: Path, genome_fasta: Path) -> str:
    contigs = {name.split()[0]: seq.upper() for name, seq in read_fasta(genome_fasta)}
    located, unlocated = [], []
    for name, seq in read_fasta(v0_fasta):
        s = seq.upper().replace("U", "T").replace("-", "")
        hit = None
        for c, g in contigs.items():
            circ = g + g[:100]
            i = circ.find(s)
            for strand, query in (("+", s), ("-", revcomp(s))):
                i = circ.find(query)
                if 0 <= i < len(g):
                    end = (i + len(s) - 1) % len(g) + 1
                    hit = {"contig": c, "start": i + 1, "end": end, "strand": strand}
                    break
            if hit:
                break
        (located if hit else unlocated).append((name, s, hit))
    any_locus = sum(any(overlaps(h, r) for r in inf.loci) for _, _, h in located)
    passing = sum(any(overlaps(h, r) for r in inf.passing()) for _, _, h in located)
    p = inf.passing()
    v2_in_v0 = sum(any(overlaps(r, h) for _, _, h in located) for r in p)
    return (f"## Comparison with v0 (infected)\n\n"
            f"- v0 candidates: {len(located) + len(unlocated)}; located on `{inf.genome}` by exact match: {len(located)}"
            f" (not located: {len(unlocated)}, e.g. consensus sequences with mismatches)\n"
            f"- located v0 candidates overlapping any SmithHunter2 locus: {any_locus}; a passing locus: **{passing}**\n"
            f"- SmithHunter2 passing loci overlapping a v0 candidate: {v2_in_v0} of {len(p)}\n")


def why_failed(r: dict) -> str:
    if r["pass"] == "1":
        return "pass"
    reasons = []
    if r["locus_type"] != "peak":
        reasons.append("background locus")
    if r["expression_pass"] != "1":
        reasons.append(f"expression ({r['samples_passing_rpm']} samples)")
    if r["ends_pass"] != "1":
        reasons.append(f"ends (5' dominant {r['five_dominant_fraction']}, isolation {r['isolation']})")
    if r.get("structural_pass", "1") != "1":
        reasons.append(f"structural ({r['structural']})")
    return "; ".join(reasons) or "?"


def section_positive(runs: list[tuple[str, Run]], positive: Path) -> str:
    """Known small RNAs (TSV: name, contig, start, end, strand; 1-based): best locus per line."""
    rows = []
    for known in read_tsv(positive):
        q = {"contig": known["contig"], "start": known["start"], "end": known["end"],
             "strand": known["strand"]}
        for label, run in runs:
            hits = [r for r in run.loci if overlaps(r, q)]
            passing = [r for r in hits if r["pass"] == "1"]
            best = max(passing or hits, key=lambda r: float(r["total_count"]), default=None)
            reads = sum(float(r["total_count"]) for r in hits)
            if best is None:
                rows.append([known["name"], label, 0, 0, "-", "-", "no locus"])
            else:
                rows.append([known["name"], label, len(hits), fmt(reads), best["locus_id"],
                             f"{best['rep_sequence']} ({fmt(float(best['total_count']))})", why_failed(best)])
    return ("## Positive control (known small RNAs)\n\nLoci overlapping each region on its strand; "
            "the best is the passing locus with most reads, or the locus with most reads.\n\n"
            + table(["small RNA", "line", "loci", "reads", "best locus", "sequence (reads)", "status"], rows) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--positive", type=Path,
                    help="TSV name, contig, start, end, strand (1-based) of known small RNAs")
    ap.add_argument("--infected", required=True, type=Path, help="results dir of the infected line")
    ap.add_argument("--tet", type=Path, help="results dir of the tetracycline-cleared line")
    ap.add_argument("--genome", default="wolbachia")
    ap.add_argument("--genome-fasta", type=Path, help="needed with --v0-fasta")
    ap.add_argument("--v0-fasta", type=Path, help="v0 smithRNA candidates (FASTA)")
    ap.add_argument("--min-samples", type=int, help="as loci.min_samples (default n - 1)")
    args = ap.parse_args()

    inf = Run(args.infected, args.genome)
    out = [f"# SmithHunter2 loci of `{args.genome}`\n",
           "## Infected line\n", section_reads(inf, "infected"), section_loci(inf, "infected"),
           section_annotation(inf, "infected"), section_sweep(inf, "infected", args.min_samples)]
    tet = Run(args.tet, args.genome) if args.tet else None
    if tet:
        out.append(section_tet(inf, tet))
    if args.positive:
        out.append(section_positive([("infected", inf)] + ([("tet", tet)] if tet else []), args.positive))
    if args.v0_fasta:
        if not args.genome_fasta:
            sys.exit("--v0-fasta needs --genome-fasta")
        out.append(section_v0(inf, args.v0_fasta, args.genome_fasta))
    print("\n".join(out))


if __name__ == "__main__":
    main()
