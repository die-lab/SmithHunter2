"""Command line entry points, one per workflow step."""

from __future__ import annotations

import argparse
import json
import os
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
                                                args.ambiguous, args.multimap, args.mode,
                                                args.evil_max_diff, args.max_hits):
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
        satellite_fraction=args.satellite_fraction, min_isolation=args.min_isolation,
        five_filter=args.five_filter, min_five_fraction=args.min_five_fraction,
        five_jitter=args.five_jitter)
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


def cmd_funnel(args):
    from .report import funnel, write_mqc_funnel, write_mqc_origin
    samples = [r["sample"] for r in read_tsv(args.libraries)]
    if len(args.fastp) != len(samples):
        sys.exit("--fastp needs one fastp JSON per sample, in the order of --libraries")
    rows = funnel(samples, args.fastp, args.libraries, args.summary, args.loci)
    write_tsv(args.out, rows, ["stage", "label"] + samples + ["total"])
    if args.mqc_funnel:
        write_mqc_funnel(args.mqc_funnel, samples, rows)
    if args.mqc_origin:
        write_mqc_origin(args.mqc_origin, samples, args.summary)
    print(" -> ".join(f"{r['stage']} {r['total']}" for r in rows))


def cmd_coverage(args):
    from .loci import load_counts
    from .report import coverage, coverage_columns
    samples, _ = load_counts(args.counts)
    n = 0
    with open_text(args.out, "wt") as fh:
        fh.write("\t".join(coverage_columns(samples)) + "\n")
        for row in coverage(args.alignments, args.counts, args.genomes, args.contigs):
            fh.write("\t".join(map(str, row)) + "\n")
            n += 1
    print(f"{n} covered positions on focus genomes")


def parse_window(text: str) -> tuple[int, int]:
    try:
        a, b = (int(x) for x in text.split(","))
    except ValueError:
        sys.exit(f"--window must be UPSTREAM_NEG,DOWNSTREAM, e.g. --window=-200,100, got {text}")
    if not a < 0 < b:
        sys.exit(f"--window must run from a negative to a positive position, got {text}")
    return a, b


def cmd_regions_bacterial(args):
    from .targets.regions import bacterial_windows, write_regions
    if not (args.genome and args.annotation):
        sys.exit("--mode bacterial needs --genome and --annotation (GFF3 with CDS features)")
    circular = (set(filter(None, args.circular_contigs.split(","))) if args.circular_contigs
                else args.topology == "circular")
    stats: dict = {}
    regions = bacterial_windows(args.genome, args.annotation, parse_window(args.window),
                                circular, args.min_window, stats)
    write_regions(args.out, regions)
    print(f"{stats['windows']} start-codon windows from {stats['cds']} CDS "
          f"({stats['clipped']} cut at contig ends, {stats['wrapped']} across the origin; "
          f"skipped: {stats['partial']} without start codon, {stats['pseudo']} pseudo, "
          f"{stats['too_short']} too short, {stats['no_contig']} contig not in the FASTA)")


def cmd_regions(args):
    if args.mode == "bacterial":
        return cmd_regions_bacterial(args)
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
    decoys = make_decoys(candidates, args.n, parse_seed(args.seed), args.random_seed,
                         keep_seed=args.keep_seed, mode=args.decoy_mode)
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


def cmd_score_sites(args):
    from .targets.pipeline import read_candidates, read_decoys
    from .targets.regions import read_regions
    from .targets.predictors import as_tuple, resolve_predictors
    from .targets.score import ScoreConfig, run_scoring, write_scores
    temperature = args.temperature
    if temperature is None:
        temperature = 37.0
        print("WARNING: --temperature not set, using 37 C. Set the temperature of the system you "
              "study (e.g. 25-28 C for many invertebrates and insect cells).", file=sys.stderr)
    genomes = set(filter(None, args.genomes.split(","))) if args.genomes else set()
    candidates = [(cid, seq) for cid, genome, seq in read_candidates(args.candidates)
                  if not genomes or not genome or genome in genomes]
    ids = {cid for cid, _ in candidates}
    decoys = [d for d in read_decoys(args.decoys) if d[1] in ids] if args.decoys else []
    classes = args.classes.split(",") if args.classes else None
    try:
        chosen = resolve_predictors(args.mode, args.predictors)
    except ValueError as exc:
        sys.exit(f"score-sites: {exc}")
    cfg = ScoreConfig(
        predictors=chosen, mode=args.mode,
        support_from=as_tuple(args.support_from) if args.support_from is not None else None,
        prefilter_by=args.prefilter_by,
        prefilter_keep_classes=as_tuple(args.prefilter_keep_classes) or (),
        report_dir=args.report_dir or (os.path.dirname(os.path.abspath(args.sites))
                                       if args.sites else ""),
        score_weak_sites=args.score_weak_sites, temperature=temperature,
        intarna_flank=args.intarna_flank, rnahybrid_flank=args.rnahybrid_flank,
        intarna_prefilter=args.intarna_prefilter, timeout=args.timeout,
        rnahybrid_require_site=not args.rnahybrid_any_duplex, resume=not args.no_resume,
        calibrate_samples=args.calibrate_samples, workers=args.workers, fdr=args.fdr,
        rnahybrid_batch=args.rnahybrid_batch, intarna_batch=args.intarna_batch)
    result = run_scoring(candidates, decoys, read_regions(args.regions), parse_seed(args.seed),
                         classes, cfg, args.workdir, args.dry_run)
    if result is None:
        return
    _, site_rows, target_rows, curves = result
    write_scores(site_rows, target_rows, curves, args.sites, args.targets, args.curves)
    print(f"{len(site_rows)} scored sites, {len(target_rows)} targets, "
          f"{sum(1 for r in site_rows if r['support'] > 0)} sites with support >= 1")


def cmd_score_windows(args):
    from .targets.bacterial import BacterialConfig, run_window_scoring, write_windows
    from .targets.pipeline import read_candidates, read_decoys
    from .targets.regions import read_regions
    temperature = args.temperature
    if temperature is None:
        temperature = 37.0
        print("WARNING: --temperature not set, using 37 C. Set the temperature of the system you "
              "study (the growth temperature of the bacterium).", file=sys.stderr)
    genomes = set(filter(None, args.genomes.split(","))) if args.genomes else set()
    candidates = [(cid, seq) for cid, genome, seq in read_candidates(args.candidates)
                  if not genomes or not genome or genome in genomes]
    ids = {cid for cid, _ in candidates}
    decoys = [d for d in read_decoys(args.decoys) if d[1] in ids]
    zone = tuple(int(x) for x in args.rbs_zone.split(","))
    cfg = BacterialConfig(
        temperature=temperature, seed_bp=args.seed_bp, no_seed=args.no_seed,
        acc_window=args.acc_window, acc_length=args.acc_length, rbs_zone=zone,
        workers=args.workers, fdr=args.fdr, timeout=args.timeout, resume=not args.no_resume)
    result = run_window_scoring(candidates, decoys, read_regions(args.regions), cfg, args.workdir,
                                args.dry_run)
    if result is None:
        return
    window_rows, target_rows, curves = result
    write_windows(window_rows, target_rows, curves, args.windows, args.targets, args.curves)
    print(f"{len(window_rows)} (candidate, window) rows, {len(target_rows)} targets, "
          f"{sum(r['pass_fdr'] for r in window_rows)} windows with q <= {args.fdr}")


def cmd_validate_targets(args):
    from collections import defaultdict as dd
    from .targets.regions import read_regions
    from .targets.validate import ROW_COLUMNS, validate
    targets = read_tsv(args.targets)
    truth = read_tsv(args.truth)
    universe = None
    if args.regions:
        genes = {r.gene for r in read_regions(args.regions)}
        universe = dd(set)
        for t in truth:
            universe[t["candidate"]] = genes
    rows, summary = validate(targets, truth, universe)
    write_tsv(args.out, rows, ROW_COLUMNS)
    print(f"{summary['found']}/{summary['pairs']} validated pairs have a scored site; "
          f"AUC by predictor: {summary['auc']}")
    for name, d in summary["detail"].items():
        print(f"  {name}: median percentile {d['median_pct']}, "
              f"in the top 5%: {d['top5']:.1%}, top 10%: {d['top10']:.1%}")


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
    s.add_argument("--mode", default="best", choices=["best", "severe", "evil"],
                   help="genome choice for reads fitting a focus and another genome")
    s.add_argument("--evil-max-diff", type=int, default=0,
                   help="evil mode: drop reads fitting an exclude genome up to this many "
                        "mismatches worse than their best focus alignment")
    s.add_argument("--max-hits", type=int, default=None,
                   help="reads with more best alignments are too_many_hits (needed when "
                        "bowtie does not apply -m, i.e. evil mode)")
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
    s.add_argument("--min-five", type=float, default=0.5, help="five-filter score: score > this")
    s.add_argument("--five-filter", default="dominant", choices=["dominant", "jitter", "score"],
                   help="5' end test: share of reads at the top position (dominant), within "
                        "--five-jitter nt of it (jitter), or the end score (score)")
    s.add_argument("--min-five-fraction", type=float, default=0.5,
                   help="dominant / jitter: minimum share of the reads")
    s.add_argument("--five-jitter", type=int, default=1)
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
    s.add_argument("--mode", choices=["ago", "bacterial"], default="ago",
                   help="bacterial: one window around the start codon of every CDS "
                        "(needs --genome and --annotation)")
    s.add_argument("--window", default="-200,100",
                   help="bacterial: UPSTREAM,DOWNSTREAM relative to the start codon (A = +1); "
                        "write it with '=' (--window=-200,100)")
    s.add_argument("--topology", choices=["linear", "circular"], default="linear",
                   help="bacterial: windows of circular contigs continue across the origin, "
                        "those of linear contigs are cut at the ends")
    s.add_argument("--circular-contigs", default="",
                   help="bacterial: comma-separated contigs that are circular (all others "
                        "linear); replaces --topology")
    s.add_argument("--min-window", type=int, default=30,
                   help="bacterial: drop windows shorter than this after cutting at contig ends")
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_regions)

    s = sub.add_parser("decoys", help="dinucleotide-shuffled decoys of the candidates")
    s.add_argument("--candidates", required=True)
    s.add_argument("--n", type=int, default=20)
    s.add_argument("--seed", default="2-8")
    s.add_argument("--random-seed", type=int, default=1)
    s.add_argument("--decoy-mode", choices=["shuffle", "keep_seed", "match_seed_gc"],
                   default="shuffle",
                   help="shuffle: whole sequence; match_seed_gc: whole sequence, seed with the "
                        "real seed's G+C count; keep_seed: only 3' of the seed")
    s.add_argument("--keep-seed", action="store_true",
                   help="shuffle only the part 3' of the seed (decoys share the real seed sites)")
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

    s = sub.add_parser("score-sites",
                       help="score seed sites with IntaRNA and RNAhybrid, q-values from decoys")
    s.add_argument("--candidates", required=True)
    s.add_argument("--decoys", required=True)
    s.add_argument("--regions", required=True)
    s.add_argument("--genomes", default="", help="only candidates of these genomes")
    s.add_argument("--seed", default="2-8")
    s.add_argument("--classes", default="")
    s.add_argument("--mode", default="ago", choices=["ago", "bacterial", "cross-kingdom"],
                   help="interaction model; sets the default predictors (ago: IntaRNA and "
                        "RNAhybrid, bacterial: IntaRNA only, cross-kingdom: both)")
    s.add_argument("--predictors", default=None,
                   help="comma list of predictors to COMPUTE (intarna,rnahybrid); a skipped one "
                        "costs nothing and its columns stay empty. Default: by --mode")
    s.add_argument("--support-from", default=None,
                   help="comma list of predictors that COUNT when calling sites (default: all "
                        "computed). Others are reported but ignored; calls are the union")
    s.add_argument("--prefilter-by", default="rnahybrid", choices=["rnahybrid"],
                   help="predictor that decides which sites IntaRNA scores when "
                        "--intarna-prefilter < 1")
    s.add_argument("--prefilter-keep-classes", default="",
                   help="comma list of site classes that IntaRNA always scores in full despite "
                        "the prefilter, e.g. 8mer")
    s.add_argument("--report-dir", default=None,
                   help="where predictor_report.tsv and prefilter_summary.tsv go "
                        "(default: the directory of --sites)")
    s.add_argument("--score-weak-sites", action="store_true",
                   help="also score 6mer and offset-6mer sites (many more runs)")
    s.add_argument("--temperature", type=float, default=None,
                   help="deg C of the interaction (IntaRNA); set it for the organism: 37 is "
                        "wrong for most invertebrates and bacteria outside the host")
    s.add_argument("--intarna-flank", type=int, default=100)
    s.add_argument("--rnahybrid-batch", type=int, default=2000,
                   help="sites per RNAhybrid process (0: one process per site); same numbers")
    s.add_argument("--intarna-batch", type=int, default=50,
                   help="windows per IntaRNA process (0: one process per site); same numbers")
    s.add_argument("--intarna-prefilter", type=float, default=1.0,
                   help="<1: IntaRNA only on this best share of sites by RNAhybrid MFE "
                        "(per class, same rule for real and decoy sites)")
    s.add_argument("--rnahybrid-flank", type=int, default=15)
    s.add_argument("--timeout", type=float, default=300,
                   help="seconds per IntaRNA / RNAhybrid call; a hang counts as a failed site")
    s.add_argument("--no-resume", action="store_true",
                   help="ignore checkpoints in the work directory and score everything again")
    s.add_argument("--rnahybrid-any-duplex", action="store_true",
                   help="keep RNAhybrid duplexes that do not cover the seed of the site")
    s.add_argument("--calibrate-samples", type=int, default=1000,
                   help="RNAcalibrate sample size per small RNA; 0 uses the 3utr_human preset")
    s.add_argument("--fdr", type=float, default=0.05)
    s.add_argument("--workers", type=int, default=1)
    s.add_argument("--workdir", default=".", help="scratch directory (must not be /tmp)")
    s.add_argument("--dry-run", action="store_true", help="only count the sites to score")
    s.add_argument("--sites", help="scored sites table")
    s.add_argument("--targets", help="scored targets table")
    s.add_argument("--curves", help="FDR curves table")
    s.set_defaults(func=cmd_score_sites)

    s = sub.add_parser("score-windows",
                       help="bacterial mode: IntaRNA of every small RNA and decoy against every "
                            "start-codon window, q-values per best interaction of a window")
    s.add_argument("--candidates", required=True)
    s.add_argument("--decoys", required=True)
    s.add_argument("--regions", required=True, help="regions.tsv made with --mode bacterial")
    s.add_argument("--genomes", default="", help="only candidates of these genomes")
    s.add_argument("--temperature", type=float, default=None,
                   help="deg C of the interaction: the growth temperature of the study")
    s.add_argument("--seed-bp", type=int, default=7,
                   help="IntaRNA seed: base pairs, at any position (6-7 usual)")
    s.add_argument("--no-seed", action="store_true", help="IntaRNA --noSeed")
    s.add_argument("--acc-window", type=int, default=150, help="IntaRNA --tAccW")
    s.add_argument("--acc-length", type=int, default=100, help="IntaRNA --tAccL")
    s.add_argument("--rbs-zone", default="-20,-1",
                   help="positions relative to the start codon counted as the RBS zone "
                        "(hit_rbs); use --rbs-zone=-20,-1")
    s.add_argument("--timeout", type=float, default=300,
                   help="seconds per IntaRNA call; a hang counts as a failed window")
    s.add_argument("--no-resume", action="store_true",
                   help="ignore the checkpoint in the work directory")
    s.add_argument("--fdr", type=float, default=0.05)
    s.add_argument("--workers", type=int, default=1)
    s.add_argument("--workdir", default=".", help="scratch directory (must not be /tmp)")
    s.add_argument("--dry-run", action="store_true", help="only count the IntaRNA runs")
    s.add_argument("--windows", help="scored windows table")
    s.add_argument("--targets", help="scored targets table")
    s.add_argument("--curves", help="FDR curves table")
    s.set_defaults(func=cmd_score_windows)

    s = sub.add_parser("validate-targets",
                       help="rank of validated small RNA -> gene pairs among scored targets")
    s.add_argument("--targets", required=True, help="targets_scored.tsv")
    s.add_argument("--truth", required=True, help="TSV with columns candidate, gene[, note]")
    s.add_argument("--regions", default="", help="regions.tsv: all its genes count as ranked")
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_validate_targets)

    s = sub.add_parser("funnel", help="reads left after each step, per sample")
    s.add_argument("--fastp", nargs="+", required=True, help="fastp JSON reports, one per sample")
    s.add_argument("--libraries", required=True)
    s.add_argument("--summary", required=True, help="origin_summary.tsv")
    s.add_argument("--loci", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--mqc-funnel", default="", help="MultiQC custom content table")
    s.add_argument("--mqc-origin", default="", help="MultiQC custom content bar graph")
    s.set_defaults(func=cmd_funnel)

    s = sub.add_parser("coverage", help="depth and 5' ends per position of focus genomes")
    s.add_argument("--alignments", required=True)
    s.add_argument("--counts", required=True)
    s.add_argument("--genomes", required=True)
    s.add_argument("--contigs", required=True)
    s.add_argument("--out", required=True)
    s.set_defaults(func=cmd_coverage)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
