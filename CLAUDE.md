# SmithHunter2 – development notes

Context for anyone (human or Claude) continuing this project. Read this first.

## Working with the maintainer

- Use English in conversation (the maintainer's preference), as in code, comments and docs.
- The maintainer co-authored the original SmithHunter (ESZlab).
- `die-lab/SmithHunter` is an earlier experimental rewrite: it is NOT the original and
  must not be used as a model. The original is https://github.com/ESZlab/SmithHunter.
- The pipeline must run on a Linux server (university machine, VPN + SSH). Development happens
  there through the Claude Code extension for VS Code (Remote-SSH). Direct calls from the
  server to api.anthropic.com are blocked, GitHub works. Ask before launching heavy jobs:
  the machine is shared.
- Development so far was done on a Windows PC without Python; module A was verified only
  through GitHub Actions. First thing on the server: create the env and run `.test`.
- Plan with literature review: https://claude.ai/artifact/D5movNUPfevGUTpjd3qiPP

## Goal

Find small RNAs produced by any genome present in a sample. The user declares every genome
with a role: `focus` (studied) or `exclude` (competes for reads, then dropped). Reads are
assigned to their best genome; ties between genomes are `ambiguous` and handled by policy.

## Design decisions (agreed)

- Snakemake; conda env per rule plus a single `environment.yml`; Python core with no
  third-party dependencies (parses bowtie SAM directly, testable without pysam).
- fastp for SE and PE. PE: keep R1, R2 only for overlap correction (`--merge` needs >= 30 bp
  overlap and lost 84% of the ~31 nt inserts in the test data). Collapse identical reads,
  align unique sequences once.
- bowtie 1, `-v 1 -k 50 -m 50 --best --strata` (no `--reorder`: bowtie 1.3.1 deadlocks with it and `-p` > 1), on one combined reference
  (contigs named `genome|contig`). Replaces v0's bowtie2 cascade mito → nuclear → mito, which
  dropped any read touching the nuclear genome even when it fit the mitochondrion better.
- Loci by genomic coordinates and strand (replaces v0 vsearch clustering by sequence identity).
- Circular contigs: overhang appended to the reference, folded back after alignment, loci may
  wrap the origin.
- Expression filter: RPM on library and raw count per sample (`min_rpm`, `min_count`);
  `min_samples` defaults to n - 1 (as in the 2025 GBE paper).
- End score: cleaned version of v0 `sharp_smith.R` (see `src/smithhunter/ends.py`).
- Module B (to do): PITA is dropped. Seed 2–8 by default (canonical seed; Bartel lab,
  McGeary et al. 2019; Plazzi et al. 2024 for smithRNAs). Site classes 8mer / 7mer-m8 /
  7mer-A1 / 6mer / offset 6mer. IntaRNA + RNAhybrid (+ miRanda optional), all scores kept,
  empirical FDR from shuffled decoy smallRNAs. Per-genome interaction mode: `ago` (seed, UTRs)
  or `bacterial` (IntaRNA around the start codon, no seed rule; CopraRNA optional).
- Precursor folding (module B): RNAfold, MFEI, p-value against shuffled sequences.

## Known problems of v0 to avoid (found reading ESZlab code)

Directory checks with `[ -f dir ]` (never cleaned, `>>` accumulates), sample counting by
substring grep (micro1 matches micro10), quadratic grep loops, PITA column 7 is dGduplex not
ddG, seed options not passed to RNAhybrid/PITA, BED 0-based coordinates passed to 1-based
seqret, `getfasta` without `-s`, sharp_smith.R parsing numbers as strings and using the loop
variable `three_score` in the filter.

## Status

- Module A implemented: `collapse`, `build-ref`, `assign`, `discover` subcommands
  (`python -m smithhunter`), Snakemake rules in `workflow/rules/`.
- Verified on the server (2026-10-09): env `smithhunter2` in `~/miniforge3`, `.test` runs end to end.
  The `.test` reads are random reads from a real sample: they only show that the workflow runs.
- `simulate` / `evaluate` subcommands: simulated library with known truth (isomiRs, non-templated
  tails, errors, degradation, tRFs, NUMTs, origin-spanning locus, close/overlapping pairs, SNP,
  heteroplasmy). Locus definition by overlap alone chained a third of the smallRNAs with
  degradation into gene-sized loci, so loci are now split at 5' end peaks (Poisson test against the
  local background, `split_peaks` in `loci.py`); the expression filter subtracts each sample's
  background. Simulations: 44/44 outcomes as expected from 0% to 30% degradation, several seeds.
- Open: on `.test` (real reads) weak satellite peaks 3-15 nt from strong loci pass the filters
  (12 -> 48 passing loci). To discuss: merge peaks close to a much stronger one.

## Next steps

1. Satellite peaks on real data (see Status), then legacy comparison with v0 output on the same data.
2. Module B (targets), then HTML report, then optional tools (MINTmap, Kraken2 on unmapped
   reads, ShortStack cross-check, AGO-CLIP chimeras).

## Commands

```bash
pytest                                            # unit tests (repo root)
cd .test && snakemake -s ../workflow/Snakefile --cores 2   # integration test
```
