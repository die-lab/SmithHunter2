# SmithHunter2 – development notes

> **Resuming? Start with `docs/HANDOFF.md`** (state, running jobs, ordered next steps, server rules),
> then `docs/module_b_development_log.md`. The notes below are older.

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

## Report and v0 comparison (2026-10-28)

- Reporting layer adapted from the colleague's `~/14_SmithHunther2/smithhunter24` (a faithful
  Snakemake port of v0 module A, not under git): `report.py` + `funnel` / `coverage` subcommands,
  `workflow/rules/report.smk`, R scripts in `workflow/scripts/`. `read_funnel.tsv` and
  `provenance.yaml` are in rule `all` (stdlib only); plots and MultiQC are rule `report`
  (envs `r.yaml`, `multiqc.yaml`; tested on the server by putting the colleague's existing R and
  MultiQC envs on PATH instead of creating new ones). The "Found 1 samples" line of MultiQC for
  the origin bar graph counts datasets: the plot has every sample.
- v0 comparison on `.test` in `~/THE_CAPTAIN/compare_v0_20261024` (a copy of the colleague's
  pipeline with the adapters set: with `adapter_r1: null` fastp autodetection left the adapter
  on most reads). Their 30 candidates: 15 match or overlap our 25 passing loci, 15 fail our
  filters (mostly ends/isolation). 15 of ours are not theirs, e.g. `mito_00356` (812 reads,
  1 mismatch on mito, 2 on a nuclear contig) removed by their nuclear filter; their sequence
  clustering splits one locus into several candidates.
- `five_score` takes few values (1, 0.5, 0.476 = two ends 2 nt apart, 0.333, ...) and the filter is
  `> min_five_score` (0.5), so it means "one 5' position holds >= n_thre of the reads". Loci with
  two main 5' ends 1-2 nt apart fail (e.g. mito_00060: 39 reads at +2, 32 at +4).
  Decided by the maintainer: a user flag, as for `assignment.mode`. `ends.five_filter`:
  `dominant` (default; share of reads at the top 5' position >= `min_five_fraction`, 0.5),
  `jitter` (share within `five_jitter` nt, default 1, of the top position), `score` (the old
  `five_score > min_five_score`). `dominant` at 0.5 passes exactly the loci the old rule passed
  (test `test_score_above_half_means_dominant_fraction_reaches_n_thre`), so defaults keep results;
  new column `five_jitter_fraction`. The 3' end still uses `three_score`.

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

## Genome-choice modes (added 2026-10-27) - READ BEFORE TOUCHING `assign.py`

`assignment.mode` in the config (`best` | `severe` | `evil`, plus `evil_max_diff`), implemented in
`classify()` (`assign.py`), `align`/`assign_origin` rules (`mapping.smk`), `common.smk` checks,
tests in `tests/test_assign.py`. Requested by the maintainer; v0 dropped any read touching the
nuclear genome even when it fit the focus genome better.
- `best` (default): fewer mismatches wins, tie = double-mapping read, `ambiguous` policy applies.
- `severe`: ties between genomes always discarded (class `ambiguous`, nothing placed).
- `evil`: as severe, and discard (class `evil_rejected`) when
  0 < (best mismatches on an exclude genome) - (best mismatches on the focus genome) <= `evil_max_diff`.
  The maintainer phrased the parameter as "focus minus filtering-out"; that difference is negative when
  the focus fits better, so the config uses its absolute value (how much worse the exclude genome fits).
  `evil_max_diff` 0 = severe; must be <= `mapping.mismatches` (checked in `common.smk`).
- Evil needs the worse strata: `align` then runs `bowtie --best -k max_hits+1` (no `--strata`, no `-m`)
  and `assign --max-hits` sets aside reads with too many best alignments. Caveat: with `-k` truncation, a
  read with more than max_hits+1 alignments in total may lose exclude hits of a worse stratum (rare).
  I assumed bowtie reports best-to-worst with `--best`; `classify` does not rely on the order, only on the
  set of hits it receives.
- Checked on `.test` copy (best/severe/evil, `evil_max_diff` 1): best = severe there because the policy is
  `report` (ties are not used either way); evil rejected 251 of 3283 mapped unique reads.
- To do / ideas: not yet in the simulator (`simulate.py` has no NUMT with mismatches; add one and test
  `evil` in `evaluate`); `severe`/`evil` also discard ties between two focus genomes (decide if that is wanted);
  the origin summary has `evil_rejected` as a new class (check downstream code and plots that list classes);
  `evil` does not look at hit counts per genome, only at the best mismatch count.

## Module B step 2 (2026-11-09)

- Done: `targets/intarna.py`, `rnahybrid.py`, `score.py`, `validate.py`; subcommands `score-sites`
  (with `--dry-run`) and `validate-targets`; rule `score_sites` (in `targets`; `predictors: []` skips it);
  tests `tests/test_score.py` (parsers on captured real output, geometry, q-values with the tools
  faked, real-tool tests skipped when the tools are not on PATH) and `tests/test_validate.py`.
- Single environment is back: dry solve of `environment.yml` (python>=3.10,<3.14 + intarna,
  rnahybrid, viennarna) gives python 3.13, IntaRNA 3.4.1, RNAhybrid 2.1.2, ViennaRNA 2.7.2.
  NOT yet created for real (it would install under ~/miniforge3, outside THE_CAPTAIN); `smithhunter2-rna`
  still works and can be dropped after creating it.
- Tool facts checked on the server: IntaRNA ranges are `from-to`; `--tRegion` must contain the
  whole interaction (guide reach), so the seed is pinned with `--seedTRange`; IntaRNA output has
  `id1` = target, `id2` = query. RNAhybrid `-f` counts from the 5' end of the guide; it returns
  `mfe 0.0 p 1.0` when the helix constraint cannot be met (treated as no interaction). RNAcalibrate
  prints `name n xi theta`; xi ~1.5-1.7, theta ~0.2, stable from k=1000.
- Simulation check (CI asserts it): planted 8mers srna_001-005 with extended pairing to guide 9-21
  (4 of them long enough) are supported (q <= 0.05 in both predictors) 4/4; bare seed sites 0/11;
  other sites 0/170 (1.1% in another seed). `evaluate` prints this ("Scored sites" block).
- Code review by a subagent (2026-11-10) found, and these are fixed with tests in
  `tests/test_review_fixes.py`: candidates without decoys got q = 0 (now no q; pooled decoys
  weighted 1/n of their candidate; no decoys at all raises); `t1` of non-canonical seeds
  (`sites.py`) was off by one, so IntaRNA silently scored nothing for `--seed 4-10`;
  `validate` gave AUC 1 for one-gene universes; no tool timeouts; a stage-2 error wiped a
  site's RNAhybrid score; `best_site_class`. Two design points it raised: the full-shuffle
  decoy null differs from the real sequence in seed composition (added `decoy_mode:
  keep_seed`), and the RNAhybrid duplex may sit on another seed match of the window (added
  `rnahybrid_require_site`; my first `at_site` flag compared the duplex start with a fully paired
  guide and was wrong for seed-only sites).
- Literature check by two subagents (2026-11-10), points acted on: p-value calibration: window-length calibration (~55 nt)
  was tried and fails (RNAcalibrate gives xi < 0 below ~100 nt), so it uses RNAcalibrate's default 500+-300 nt
  (`CALIBRATION_LENGTH`); the region-length version was replaced; `fdr_pseudocount` (+1) in q-values;
  decoy modes `shuffle` / `match_seed_gc` / `keep_seed` (`--decoy-mode`); README section
  "Statistical notes and limits". NOT done: seed-composition-matched decoys beyond GC count,
  measuring the sensitivity to the IntaRNA flank, recall of two-stage mode on real data. Citations
  that could not be verified were left out (Fu & Qian 2014, Min & Yoon 2010, Cheng/Yu; PITA flank sizes).
- Checkpoints (`score.py`): `checkpoint_stage1.pkl` / `checkpoint_scores.pkl` in `--workdir`, keyed by
  a signature of guides + decoys + target sequences + settings. Tests run `run_scoring` in a tmp cwd
  (autouse fixture) so checkpoints never leak between tests. Developing while a long job runs: copy
  `src` and `tests` to `.scratch_solve/dev` and give it its own `pytest.ini` (the repo's
  `pyproject.toml` otherwise puts the repo `src` first on the path); deploy after the tests pass.
- Costs measured on the shared server: IntaRNA 0.14 s/site (flank 150), 0.05 s (flank 60);
  RNAhybrid 0.009 s/site. Hence `intarna_flank` default 100 and the two-stage `intarna_prefilter`.
- Calibration data for bacterial -> host small RNAs: `~/THE_CAPTAIN/calib_bacterial_rnai`, with
  `README.md` there (what was chosen and why, accessions, download state).
- Rule for this server: never write outside THE_CAPTAIN, also not to /tmp (pytest: use
  `--basetemp` inside THE_CAPTAIN; conda: `CONDA_PKGS_DIRS` inside THE_CAPTAIN; `pkill -f` can
  kill your own ssh shell: kill by pid after checking `/proc/<pid>/cwd`).

## Next steps

1. Single environment again: pin Python in `environment.yml` (e.g. `>=3.10,<3.14`) and add
   `intarna>=3.4`, `rnahybrid`, `viennarna>=2.6` with minimum versions; check with a dry solve
   (`mamba env create --dry-run -f environment.yml`). If it solves, drop `smithhunter2-rna`.
2. (done 2026-11-09) Module B step 2. Left: run the calibration on the bacterial -> host data
   (`calib_bacterial_rnai`), extend `evaluate` to check q-values on the simulation, miRanda.
3. Then: precursor folding (step 3), bacterial mode (step 4), plots of module B in the report.
4. Module A: run on a full-size real library to tune thresholds, then legacy comparison with v0.
   Optional tools later (MINTmap, Kraken2 on unmapped reads, ShortStack cross-check, AGO-CLIP).

## Commands

```bash
source ~/miniforge3/etc/profile.d/conda.sh && conda activate smithhunter2
pytest                                            # unit tests (repo root)
cd .test && snakemake -s ../workflow/Snakefile --cores 2   # integration test (module A)
snakemake -s ../workflow/Snakefile --cores 2 --sdm conda report   # plots + MultiQC (from .test)
# simulation benchmark, modules A and B (outside the repo, e.g. in a scratch dir):
PYTHONPATH=src python -m smithhunter simulate --outdir /tmp/sim
cd /tmp/sim && snakemake -s <repo>/workflow/Snakefile --cores 2 && \
  snakemake -s <repo>/workflow/Snakefile --cores 2 targets && \
  PYTHONPATH=<repo>/src python -m smithhunter evaluate
```
