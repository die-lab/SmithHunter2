# SmithHunter2 – development notes

Context for anyone (human or Claude) continuing this project. Read this first.

## Working with the maintainer

- The maintainer writes in Italian; answer in Italian. Code, comments and docs are in English.
- The maintainer co-authored the original SmithHunter (ESZlab).
- `die-lab/SmithHunter` is an earlier experimental rewrite: it is NOT the original and
  must not be used as a model. The original is https://github.com/ESZlab/SmithHunter.
- The pipeline must run on a Linux server (university machine, reached over SSH).
- Plan with literature review: https://claude.ai/artifact/D5movNUPfevGUTpjd3qiPP

## Goal

Find small RNAs produced by any genome present in a sample. The user declares every genome
with a role: `focus` (studied) or `exclude` (competes for reads, then dropped). Reads are
assigned to their best genome; ties between genomes are `ambiguous` and handled by policy.

## Design decisions (agreed)

- Snakemake; conda env per rule plus a single `environment.yml`; Python core with no
  third-party dependencies (parses bowtie SAM directly, testable without pysam).
- fastp for SE and PE (PE: merge pairs). Collapse identical reads, align unique sequences once.
- bowtie 1, `-v 1 -k 50 -m 50 --best --strata --reorder`, on one combined reference
  (contigs named `genome|contig`). Replaces v0's bowtie2 cascade mito → nuclear → mito, which
  dropped any read touching the nuclear genome even when it fit the mitochondrion better.
- Loci by genomic coordinates and strand (replaces v0 vsearch clustering by sequence identity).
- Circular contigs: overhang appended to the reference, folded back after alignment, loci may
  wrap the origin.
- Expression filter: RPM on library, `min_samples` defaults to n - 1 (as in the 2025 GBE paper).
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
- Written on a machine without Python or Linux tools: verified only through GitHub Actions
  (`.github/workflows/ci.yml`: unit tests + full run on `.test`). Check the latest CI run.

## Next steps

1. Run `.test` on the server; look at `origin_summary.tsv` and `loci.tsv` with the maintainer.
2. Simulated reads with known smallRNAs, tRNA fragments and NUMTs inserted, to measure recovery.
3. Legacy comparison with v0 output on the same data.
4. Module B (targets), then HTML report, then optional tools (MINTmap, Kraken2 on unmapped
   reads, ShortStack cross-check, AGO-CLIP chimeras).

## Commands

```bash
pytest                                            # unit tests (repo root)
cd .test && snakemake -s ../workflow/Snakefile --cores 2   # integration test
```
