# Improvement plan for the SmithHunter2 workflow (module B and its validation)

Written 2026-11-10 by the manager (the main Claude session). Each task is independent: it has
its own working copy and data directory, a fixed deliverable and an acceptance test. Nobody but
the manager edits the repository; agents deliver patches, the manager merges them in the order
below and reruns the whole test suite after every merge.

## Why

What the first real positive control (40 human miRNAs, 2,673 reporter-validated pairs) showed:

1. The best-MFE ranking of validated genes is good (AUC 0.79) but equal to counting canonical seed
   sites (0.797): the site energy adds nothing we can show at gene level.
2. The site-level decoy q-values almost never fire on true targets (2 supported sites in 1.7M).
3. The run is slow: 98 CPU-hours for 3.5M sites without the prefilter, ~5 h per 40-miRNA control.
4. The bacterial side has no positive control yet (Pseudomonas: none recoverable by a seed model;
   Wolbachia WsnRNA-46: nearly absent from the deposited reads).
5. The bacterium-to-bacterium mode (`bacterial`) is not implemented, and the interaction modes
   differ by organism (development log, section 6).

## Tasks

| id | task | deliverable | acceptance | CPU |
|---|---|---|---|---|
| `jobs` | resume the running work: finish and interpret the Wolbachia CDS run and the miRNA keep_seed run; compare decoy modes on the positive control | report + results added to `calib_bacterial_rnai/README.md` | numbers for both runs, comparison table shuffle vs keep_seed, reproducible commands | read-only |
| `genelevel` | test a gene-level statistic (aggregate of the site scores of a gene) against the same decoys, so UTR length and site number are in the null | `targets/genelevel.py` + tests + experiment report, or a documented negative result | beats the site-count baseline within a UTR-length decile by >= 0.03 AUC, otherwise not merged | 4 workers |
| `speed` | batch RNAhybrid (many windows per process) and investigate IntaRNA batching; cut Python overhead | patch to `rnahybrid.py` / `score.py`, benchmark | >= 5x faster RNAhybrid stage, results identical to the current code on 2,000 random sites | 4 workers |
| `bactmode` | implement the `bacterial` mode end to end (start-codon windows, IntaRNA only, window-level q-values) | patch + tests + docs section | unit tests and a synthetic planted-pair check with the real tool | 2 workers |
| `bactdata` | build the bacterial validation set (E. coli K-12 trans-encoded sRNAs with validated targets, sources for every pair) | `calib_bacterial_rnai/bact_ecoli/` (genome, candidates, truth, SOURCES.md) | >= 30 pairs from >= 10 sRNAs, each with a verified source, uncertain ones flagged | none |
| `transfect` | a second positive control WITH negatives: miRNA transfection expression changes; evaluation script | data + `evaluate_repression.py` + a first result | >= 8 miRNAs with per-gene log2FC from a cited public dataset; result vs the site-count baseline | 4 workers |
| `verify` | reproducibility: build the single environment, run the whole test suite, the CI simulation and a workflow dry run from a clean state; check README commands | report with exact commands and discrepancies | environment solves and installs inside THE_CAPTAIN; tests and CI simulation pass | 2 workers |

## Rules for every agent

Write only inside THE_CAPTAIN (never `/tmp`, never the home config, no installs into
`~/miniforge3`); the repository is read-only for agents (work in
`SmithHunter2/.scratch_solve/agents/<id>/`); `nice -n 10`; check `uptime` before heavy jobs
(wait when the load is above 40); kill only your own processes, never with `pkill -f`; do not
touch the running jobs `run_mirna` and `run_wolb_cds`; report negative results as results.

## Merge order and conflicts

1. `speed` (touches `rnahybrid.py`, `score.py`), 2. `bactmode` (new files + `regions.py`,
`targets.smk`, `cli.py`), 3. `genelevel` (new module, small hook in `score.py`), 4. data tasks
(`bactdata`, `transfect`) feed the next validation round, 5. `verify` last, on the merged code.
The manager rewrites `docs/module_b_development_log.md` after each merge; agents never edit it.

## After the merge

Re-run the bacterial and the miRNA controls with the merged code; fill the final documentation
from the development log (section 7 checklist).
