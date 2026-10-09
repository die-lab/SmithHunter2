"""Command line entry points, one per workflow step."""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict

from . import __version__
from .assign import Contig, assign
from .collapse import collapse, write_outputs
from .loci import READ_COLUMNS, discover, locus_columns, write_fasta
from .reference import build_reference, parse_genomes, write_tables
from .seqio import open_text, read_fasta_names, read_tsv, write_tsv


def cmd_collapse(args):
    samples: dict[str, list[str]] = {}
    for item in args.sample:
        name, _, path = item.partition("=")
        if not path:
            sys.exit(f"--sample expects NAME=PATH, got {item!r}")
        samples.setdefault(name, []).append(path)
    names, ordered, stats = collapse(samples, args.min_length, args.max_length)
    write_outputs(names, ordered, stats, args.fasta, args.counts, args.libraries)
    print(f"{len(ordered)} unique sequences from {sum(s['reads_kept'] for s in stats)} reads")


def cmd_build_ref(args):
    genomes = parse_genomes(json.loads(args.spec), args.min_length, args.max_length)
    overhang = max(g.max_length for g in genomes)
    contigs = build_reference(genomes, args.fasta, overhang)
    write_tables(genomes, contigs, args.genomes, args.contigs)
    print(f"{len(contigs)} contigs from {len(genomes)} genomes")


def cmd_assign(args):
    contigs = {r["ref"]: Contig(r["ref"], r["genome"], r["contig"], int(r["length"]),
                                r["circular"] == "1") for r in read_tsv(args.contigs)}
    roles = {g["name"]: g["role"] for g in read_tsv(args.genomes)}
    with open_text(args.counts) as fh:
        samples = fh.readline().rstrip("\n").split("\t")[3:]
        table = {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            table[f[0]] = (int(f[2]), [int(x) for x in f[3:]])
    excess = set(read_fasta_names(args.excess)) if args.excess else set()

    summary: dict[tuple, list[int]] = defaultdict(lambda: [0] * len(samples))
    lengths: dict[tuple, list[int]] = defaultdict(lambda: [0] * len(samples))
    with open_text(args.alignments, "wt") as aln, open_text(args.origin, "wt") as org:
        aln.write("read_id\tgenome\tcontig\tstart\tend\tstrand\tweight\tmismatches\n")
        org.write("read_id\tclass\tgenomes\n")
        for rid, cls, genomes, placed in assign(args.sam, excess, list(table), contigs, roles,
                                                args.ambiguous, args.multimap):
            label = "+".join(genomes)
            org.write(f"{rid}\t{cls}\t{label}\n")
            for genome, hit, weight in placed:
                c = contigs[hit.ref]
                aln.write(f"{rid}\t{genome}\t{c.contig}\t{hit.start}\t{hit.end}\t{hit.strand}\t"
                          f"{weight:.6g}\t{hit.mismatches}\n")
            length, row = table[rid]
            acc = summary[(cls, label)]
            lacc = lengths[(cls, label, length)]
            for i, c in enumerate(row):
                acc[i] += c
                lacc[i] += c

    write_tsv(args.summary, [{"class": k[0], "genomes": k[1], **dict(zip(samples, v))}
                             for k, v in sorted(summary.items())], ["class", "genomes", *samples])
    write_tsv(args.lengths, [{"class": k[0], "genomes": k[1], "length": k[2], **dict(zip(samples, v))}
                             for k, v in sorted(lengths.items())],
              ["class", "genomes", "length", *samples])


def cmd_discover(args):
    samples, rows, reads = discover(
        args.alignments, args.counts, args.libraries, args.genomes, args.contigs,
        merge_gap=args.merge_gap, min_rpm=args.min_rpm, min_samples=args.min_samples,
        n_thre=args.n_thre, penalty=args.penalty, min_five=args.min_five, min_three=args.min_three)
    write_tsv(args.loci, rows, locus_columns(samples))
    write_tsv(args.reads, reads, READ_COLUMNS)
    write_fasta(args.fasta, rows)
    print(f"{len(rows)} loci, {sum(r['pass'] for r in rows)} passing filters")


def main(argv=None):
    p = argparse.ArgumentParser(prog="smithhunter", description=__doc__)
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("collapse", help="count unique sequences across samples")
    s.add_argument("--sample", nargs="+", required=True, metavar="NAME=FASTQ")
    s.add_argument("--min-length", type=int, required=True)
    s.add_argument("--max-length", type=int, required=True)
    s.add_argument("--fasta", required=True)
    s.add_argument("--counts", required=True)
    s.add_argument("--libraries", required=True)
    s.set_defaults(func=cmd_collapse)

    s = sub.add_parser("build-ref", help="combine all genomes into one reference")
    s.add_argument("--spec", required=True, help="JSON of the 'genomes' config block")
    s.add_argument("--min-length", type=int, required=True)
    s.add_argument("--max-length", type=int, required=True)
    s.add_argument("--fasta", required=True)
    s.add_argument("--genomes", required=True)
    s.add_argument("--contigs", required=True)
    s.set_defaults(func=cmd_build_ref)

    s = sub.add_parser("assign", help="assign reads to their genome of origin")
    s.add_argument("--sam", required=True)
    s.add_argument("--excess", default="")
    s.add_argument("--counts", required=True)
    s.add_argument("--genomes", required=True)
    s.add_argument("--contigs", required=True)
    s.add_argument("--ambiguous", default="report")
    s.add_argument("--multimap", default="fractional")
    s.add_argument("--alignments", required=True)
    s.add_argument("--origin", required=True)
    s.add_argument("--summary", required=True)
    s.add_argument("--lengths", required=True)
    s.set_defaults(func=cmd_assign)

    s = sub.add_parser("discover", help="define, quantify, score and annotate loci")
    s.add_argument("--alignments", required=True)
    s.add_argument("--counts", required=True)
    s.add_argument("--libraries", required=True)
    s.add_argument("--genomes", required=True)
    s.add_argument("--contigs", required=True)
    s.add_argument("--merge-gap", type=int, default=0)
    s.add_argument("--min-rpm", type=float, default=5.0)
    s.add_argument("--min-samples", type=int, default=None)
    s.add_argument("--n-thre", type=float, default=0.5)
    s.add_argument("--penalty", type=float, default=0.1)
    s.add_argument("--min-five", type=float, default=0.5)
    s.add_argument("--min-three", type=float, default=0.0)
    s.add_argument("--loci", required=True)
    s.add_argument("--reads", required=True)
    s.add_argument("--fasta", required=True)
    s.set_defaults(func=cmd_discover)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
