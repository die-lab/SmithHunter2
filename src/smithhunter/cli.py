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
        merge_gap=args.merge_gap, min_rpm=args.min_rpm, min_count=args.min_count,
        min_samples=args.min_samples,
        n_thre=args.n_thre, penalty=args.penalty, min_five=args.min_five, min_three=args.min_three,
        peaks=not args.no_peaks, peak_window=args.peak_window, peak_pvalue=args.peak_pvalue,
        background_flank=args.background_flank, satellite_distance=args.satellite_distance,
        satellite_fraction=args.satellite_fraction, min_isolation=args.min_isolation)
    write_tsv(args.loci, rows, locus_columns(samples))
    write_tsv(args.reads, reads, READ_COLUMNS)
    write_fasta(args.fasta, rows)
    print(f"{len(rows)} loci, {sum(r['pass'] for r in rows)} passing filters")


def cmd_simulate(args):
    from .simulate import Params, simulate
    params = Params(**{k: v for k, v in vars(args).items() if k in Params.__dataclass_fields__})
    sim = simulate(params)
    reads = sum(sum(s.counts) for s in sim.sources)
    print(f"{len(sim.samples)} samples, {reads} reads, {len(sim.sources)} sources written to "
          f"{args.outdir}")


def cmd_evaluate(args):
    from .evaluate import run
    print(run(args.truth, args.results, args.out), end="")


def cmd_regions(args):
    from .targets.regions import regions_from_genome, regions_from_transcripts, write_regions
    wanted = args.regions.split(",")
    if args.genome:
        if not args.annotation:
            sys.exit("--genome needs --annotation (GFF3)")
        regions = regions_from_genome(args.genome, args.annotation, wanted, args.whole_if_missing)
    elif args.transcripts:
        regions = regions_from_transcripts(args.transcripts, args.table, wanted,
                                           args.whole_if_missing)
    else:
        sys.exit("give --genome and --annotation, or --transcripts")
    write_regions(args.out, regions)
    print(f"{len(regions)} regions, {sum(len(r.sequence) for r in regions)} nt")


def parse_seed(text: str) -> tuple[int, int]:
    a, b = (int(x) for x in text.split("-"))
    if not 1 <= a < b:
        sys.exit(f"--seed must be START-END with 1 <= START < END, got {text}")
    return a, b


def cmd_decoys(args):
    from .targets.pipeline import read_candidates
    from .targets.shuffle import make_decoys
    candidates = [(cid, seq) for cid, _, seq in read_candidates(args.candidates)]
    decoys = make_decoys(candidates, args.n, parse_seed(args.seed), args.random_seed)
    with open_text(args.out, "wt") as fh:
        for did, _, seq in decoys:
            fh.write(f">{did}\n{seq}\n")
    print(f"{len(decoys)} decoys for {len(candidates)} candidates")


def cmd_sites(args):
    from .targets.pipeline import read_candidates, read_decoys, run_sites, write_outputs
    from .targets.regions import read_regions
    genomes = set(filter(None, args.genomes.split(","))) if args.genomes else set()
    candidates = [(cid, seq) for cid, genome, seq in read_candidates(args.candidates)
                  if not genomes or not genome or genome in genomes]
    ids = {cid for cid, _ in candidates}
    decoys = [d for d in read_decoys(args.decoys) if d[1] in ids] if args.decoys else []
    regions = read_regions(args.regions)
    classes = args.classes.split(",") if args.classes else None
    result = run_sites(candidates, decoys, regions, parse_seed(args.seed), classes)
    write_outputs(*result, args.sites, args.targets, args.fdr)
    print(f"{len(candidates)} candidates, {len(result[1])} sites, {len(result[2])} targets")


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
    s.add_argument("--min-count", type=float, default=1)
    s.add_argument("--min-samples", type=int, default=None)
    s.add_argument("--no-peaks", action="store_true",
                   help="loci are whole blocks of overlapping reads (no 5' peak splitting)")
    s.add_argument("--peak-window", type=int, default=2)
    s.add_argument("--peak-pvalue", type=float, default=0.001)
    s.add_argument("--background-flank", type=int, default=50)
    s.add_argument("--satellite-distance", type=int, default=10)
    s.add_argument("--satellite-fraction", type=float, default=0.1)
    s.add_argument("--n-thre", type=float, default=0.5)
    s.add_argument("--penalty", type=float, default=0.1)
    s.add_argument("--min-five", type=float, default=0.5)
    s.add_argument("--min-three", type=float, default=0.0)
    s.add_argument("--min-isolation", type=float, default=0.5)
    s.add_argument("--loci", required=True)
    s.add_argument("--reads", required=True)
    s.add_argument("--fasta", required=True)
    s.set_defaults(func=cmd_discover)

    s = sub.add_parser("simulate", help="simulate a small RNA library with known truth")
    s.add_argument("--outdir", default="sim")
    s.add_argument("--seed", type=int, default=1)
    s.add_argument("--samples", type=int, default=3)
    s.add_argument("--reads", type=int, default=200_000, help="mean reads per sample")
    s.add_argument("--read-length", type=int, default=75)
    s.add_argument("--paired", action="store_true")
    s.add_argument("--mito-length", type=int, default=16_000)
    s.add_argument("--nuclear-length", type=int, default=2_000_000)
    s.add_argument("--bacterium-length", type=int, default=300_000)
    s.add_argument("--srnas", type=int, default=30, help="ordinary mitochondrial smallRNAs")
    s.add_argument("--degradation", type=float, default=0.08,
                   help="fraction of reads from mitochondrial rRNA and mRNA degradation")
    s.add_argument("--error-rate", type=float, default=0.002)
    s.add_argument("--max-hits", type=int, default=50)
    s.add_argument("--ambiguous", default="report")
    s.set_defaults(func=cmd_simulate)

    s = sub.add_parser("evaluate", help="compare results of a simulated run with its truth")
    s.add_argument("--truth", default="truth")
    s.add_argument("--results", default="results")
    s.add_argument("--out", default="evaluation")
    s.set_defaults(func=cmd_evaluate)

    s = sub.add_parser("regions", help="target regions from a genome + GFF3 or transcripts")
    s.add_argument("--genome")
    s.add_argument("--annotation")
    s.add_argument("--transcripts")
    s.add_argument("--table", default="", help="transcript, region, start, end (1-based)")
    s.add_argument("--regions", default="three_prime_UTR",
                   help="comma-separated: three_prime_UTR, five_prime_UTR, CDS, transcript")
    s.add_argument("--whole-if-missing", action="store_true",
                   help="use the whole transcript when the wanted regions are missing")
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_regions)

    s = sub.add_parser("decoys", help="dinucleotide-shuffled decoys of the candidates")
    s.add_argument("--candidates", required=True)
    s.add_argument("--n", type=int, default=20)
    s.add_argument("--seed", default="2-8")
    s.add_argument("--random-seed", type=int, default=1)
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_decoys)

    s = sub.add_parser("sites", help="seed sites of candidates and decoys, site-count FDR")
    s.add_argument("--candidates", required=True)
    s.add_argument("--decoys", default="")
    s.add_argument("--regions", required=True)
    s.add_argument("--genomes", default="", help="only candidates of these genomes")
    s.add_argument("--seed", default="2-8")
    s.add_argument("--classes", default="")
    s.add_argument("--sites", required=True)
    s.add_argument("--targets", required=True)
    s.add_argument("--fdr", required=True)
    s.set_defaults(func=cmd_sites)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
