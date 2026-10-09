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
- On `.test` (real reads) peaks produced weak satellites 3-15 nt from strong loci and clusters of
  comparable peaks (imprecise 5' ends), and five_score inside a +-2 nt peak is always ~1.
  Fixed with satellite merging (<= 10 nt, < 10% of the reads) and the `isolation` filter (share of
  reads starting within 10 nt of the peak, >= 0.5). `.test`: 25 passing loci, 9 of the 12 found by
  overlap blocks; the 3 lost have 0-3 reads per sample within +-2 nt of their 5' end.

## Where we stopped (2026-10-09)

- Module B step 1 is done and pushed (`834a2b9`, CI green): `regions`, `decoys`, `sites`
  subcommands, rule `targets` (run separately from module A, after reviewing candidates).
  Target sets use the general seed logic (`ago`) by default; `mode: bacterial` is opt-in and
  only prints "not implemented" for now. The simulator plants 15 seed sites in nuclear 3'UTRs;
  `evaluate` finds 15/15. Site-count FDR is ~1 by nature (a seed match alone is not a target):
  ranking needs step 2.
- RNA tools are installed in a separate env `smithhunter2-rna` (IntaRNA 3.4.1, RNAhybrid 2.1.2,
  ViennaRNA 2.7.2; recipe `workflow/envs/rna.yaml`). Installing them into `smithhunter2`
  (Python 3.14) silently gave IntaRNA 1.2.5 and ViennaRNA 2.4.7 (py27), so they were removed.
- Pushing: `git push git@github.com:die-lab/SmithHunter2.git main` (SSH; the HTTPS remote has no
  credentials). Repo-local `user.email` is thepokemonmuia@gmail.com (the other address is blocked
  by GitHub email privacy). `gh` is not installed: CI status via
  `curl https://api.github.com/repos/die-lab/SmithHunter2/actions/runs?per_page=3`.

## Next steps

1. Single environment again: pin Python in `environment.yml` (e.g. `>=3.10,<3.14`) and add
   `intarna>=3.4`, `rnahybrid`, `viennarna>=2.6` with minimum versions; check with a dry solve
   (`mamba env create --dry-run -f environment.yml`). If it solves, drop `smithhunter2-rna`.
2. Module B step 2 (`docs/module_b_design.md` sections 4.2-4.3): wrap IntaRNA (whole region for
   accessibility, `--tRegion` on the site, seed pinned to guide 2-8) and RNAhybrid (window around
   the site, `-f` matching the seed, `-d xi,theta` from RNAcalibrate on the same regions); score
   real and decoy sites; q-values with `targets/fdr.py:qvalues`; `support` column. Add stored
   tool outputs as parser test fixtures; extend the simulation check.
3. Then: precursor folding (step 3), bacterial mode (step 4), report.
4. Module A: run on a full-size real library to tune thresholds, then legacy comparison with v0.
   Optional tools later (MINTmap, Kraken2 on unmapped reads, ShortStack cross-check, AGO-CLIP).

## Commands

```bash
source ~/miniforge3/etc/profile.d/conda.sh && conda activate smithhunter2
pytest                                            # unit tests (repo root)
cd .test && snakemake -s ../workflow/Snakefile --cores 2   # integration test (module A)
# simulation benchmark, modules A and B (outside the repo, e.g. in a scratch dir):
PYTHONPATH=src python -m smithhunter simulate --outdir /tmp/sim
cd /tmp/sim && snakemake -s <repo>/workflow/Snakefile --cores 2 && \
  snakemake -s <repo>/workflow/Snakefile --cores 2 targets && \
  PYTHONPATH=<repo>/src python -m smithhunter evaluate
```
