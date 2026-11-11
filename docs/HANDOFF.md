# SmithHunter2 / module B: handoff (written 2026-11-10 12:30 CET)

Read this first when resuming in a new chat. It states where everything is, what is finished,
what is still running or half done, and what to do next, in order. Companion documents:
`docs/module_b_development_log.md` (the logic: decisions with evidence, experiments,
corrections, open questions), `docs/improvement_plan.md` (the agent plan), `CLAUDE.md`
(older development notes), `calib_bacterial_rnai/README.md` (data and calibration runs).

## 0. Rules that apply to everything (from the user's global instructions)

- Server: connect with the ssh command in your global instructions (the `^` in the user name is
  part of the name); all work happens in `~/THE_CAPTAIN` on the server. (The connection details are
  deliberately not in this file because the repository is public; the copy
  `~/THE_CAPTAIN/HANDOFF_SmithHunter2.md` outside the repository has them.)
- **Write only inside `THE_CAPTAIN`** (also not to `/tmp`, `~/.cache`, `~/.conda`, `~/miniforge3`);
  outside it is read-only. Anything outside needs the user's explicit confirmation THREE times,
  in three separate messages. Use absolute paths. `pytest` needs `--basetemp` inside THE_CAPTAIN
  and `-p no:cacheprovider`; conda needs `CONDA_PKGS_DIRS` inside THE_CAPTAIN; export `TMPDIR`
  inside THE_CAPTAIN.
- The server is shared (64 cores, load 45-60 most of the time): `nice -n 10`, ask before heavy
  jobs, never `pkill -f` (it kills your own ssh shell); kill by PID after checking `/proc/<pid>/cwd`.
- ssh quoting: apostrophes inside a single-quoted ssh command break it; put scripts in a file
  and pipe it (`ssh ... 'bash -s' < script.sh`). Backslash-n inside text sent through some tools
  becomes a real newline: build such strings with `chr(10)` or edit with an editor tool.
- Do not commit or push without being asked. The user's commit e-mail is the repo-local
  `thepokemonmuia@gmail.com`; commit trailer: `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`.
- The focus genome of a study is the organism the paper is about (not always mitochondrial);
  test the bacterial side as much as the mitochondrial one.

## 1. Where things are

| what | path (under `~/THE_CAPTAIN`) |
|---|---|
| repository (git) | `SmithHunter2/` branch **`module-b-step2`**, commit `ceeb274` (24 files), **not pushed**; `main` untouched (`bf4c098`) |
| calibration data, runs, truth files | `calib_bacterial_rnai/` (24 GB; see its README) |
| clean single environment (built by an agent, Python 3.13, IntaRNA 3.4.1, RNAhybrid 2.1.2, ViennaRNA 2.7.2) | `envs/sh2_clean` (package cache `envs/pkgs_clean`; an agent used `HOME=envs/fakehome` to keep mamba from writing to `~/.conda`) |
| existing user envs | `~/miniforge3/envs/smithhunter2` (Python 3.14, no RNA tools) and `smithhunter2-rna` (the tools); activate both: `PATH=$HOME/miniforge3/envs/smithhunter2/bin:$HOME/miniforge3/envs/smithhunter2-rna/bin:$PATH` |
| agent working copies and deliveries | `SmithHunter2/.scratch_solve/agents/<id>/` (git-ignored, ~900 MB with dev copies) |
| simulations | `sim_ci_clean/`, `sim_step2_ext/`, `sim_step2_test/` |

Test command (100 pass; 97 pass + 3 skip when the RNA tools are not on PATH):
`cd SmithHunter2; export PATH=...(both envs)...; PYTHONPATH=src python -m pytest -q -p no:cacheprovider --basetemp=$PWD/.scratch_solve/pytmp/t tests`

## 2. What is finished (in the commit)

- Module A (discovery) as before; module B step 1 (regions, decoys, seed sites).
- **Module B step 2**: `score-sites` (IntaRNA anchored on the seed with `--seedQRange/--seedTRange/--tRegion`;
  RNAhybrid with `-f` helix constraint, duplex must cover the seed; RNAcalibrate p-values at the
  tool's default lengths; optional two-stage prefilter `intarna_prefilter`; decoy modes `shuffle`,
  `match_seed_gc`, `keep_seed`; q-values per predictor and class, pooled and per candidate, decoy
  weights 1/n, +1 pseudocount, `support`), checkpoint/resume (`checkpoint_stage1.pkl`,
  `checkpoint_scores.pkl`), timeouts, per-site error tolerance (stop above 5%), progress/ETA,
  `validate-targets` (rank, percentile, AUC, top-5/10%), simulator with extended-pairing targets,
  CI assertions. Temperature defaults to 37 C with a loud warning.
- Single `environment.yml` again (verified by a clean build: 82 s).
- Docs: README "Scoring sites" and "Statistical notes and limits", design doc, development log.

## 3. Key results (all numbers are reproducible from `calib_bacterial_rnai/README.md` and the log)

1. **Simulation (CI)**: planted extended-pairing sites supported 4/4, bare seed matches 0/11, other 0/170.
2. **miRNA positive control** (40 human miRNAs, 2,673 miRTarBase reporter pairs, 19,634 genes, shuffle decoys):
   2,157/2,673 pairs have a strong seed site; gene-level AUC of best RNAhybrid MFE 0.791; **number of strong
   sites alone 0.797**, UTR length alone 0.644, MFE within UTR-length decile 0.747. Only 2 of 1.7M sites get
   `support` (the site-level decoy q-values have almost no power on true targets). IntaRNA covered only the
   best 15% of sites (prefilter), so its ranks are not comparable.
3. **Pseudomonas PA14 OMV sRNAs vs human 3'UTRs**: the published sRNA52320 is recovered by module A;
   0/5 "validated" genes have a strong site (predicted by BLASTN/miRanda, not seed-based); one supported
   candidate (a paralog, RNAhybrid MFE -59..-61, probably repeat-like complementarity, not checked).
4. **Wolbachia**: module A on the 4 Aag2 libraries works (wMel reads only in infected samples; 570 loci pass).
   The published WsnRNA-46 is nearly absent from the deposited reads (4 reads in 37M). Its Dhc site is in the
   CDS (per the paper), so the 3'UTR run was uninformative; a CDS run was started (section 4).
5. **Bacterial mode (agent `bactmode`)**: works end to end on synthetic data; 16 bp planted pairing supported,
   12 bp not (q 0.25), 14 bp at the floor 0.05 with 20 decoys; 0/200 random windows called. ~0.15-0.2 s per
   IntaRNA run -> ~100 CPU-h for a 4,000-gene genome with 30 sRNAs and 20 decoys.
6. **Bacterial validation set (agent `bactdata`)**: *E. coli* K-12, 69 pairs from 26 sRNAs (41 strong, 20 moderate,
   8 weak evidence; 33 verified in full text), `calib_bacterial_rnai/bact_ecoli/`; *Salmonella* 13 pairs (weaker).
7. **Speed (agent `speed`, its own measurements under load ~52)**: batching RNAhybrid is 24.9x faster per site
   (11.79 -> 0.474 ms/site, 5,000 real windows) with bit-identical results (mfe, p, pos, span, at_site) and
   identical end-to-end outputs on the simulation; the 6.7M-site stage would drop from ~2 h to ~5 min on
   12 workers (extrapolated). IntaRNA: 93% of the ~100-135 ms/site is the accessibility of the target window;
   batching helps only when several small RNAs share a window (7x with `keep_seed` decoys, because IntaRNA
   computes the accessibility once per target), 1.1x with `shuffle` decoys. Padding windows to a fixed layout
   is not exact; a smaller flank is faster but changes energies (flank 100/75/50/30 -> 135/89/56/30 ms,
   mean |dE| up to 0.8 kcal/mol). Side findings: `make_decoys` is not reproducible across processes unless
   `PYTHONHASHSEED` is fixed (2,533 vs 3,106 sites on identical code: **fix: use a deterministic order
   instead of set iteration; do this before comparing runs**); RNAcalibrate aborts ("Error in nextSQ") on a
   subset containing a 17,983-nt sequence.
8. **Literature checks** (subagents): seed rules supported; IntaRNA has no miRNA benchmark; FDR is
   target-decoy-style, conservative on average, not proven; `support` is a count of correlated votes;
   smithRNA evidence for AGO and the seed rule (2-8 vs 4-10) is weak; RNAcalibrate fails below ~100 nt windows
   (kept at default 500+-300 nt).

## 4. Running or unfinished at the time of writing (check these first)

| item | state at 12:20 | where / how to check |
|---|---|---|
| miRNA control, mode 2 `keepseed` (8.2M sites) | RNAhybrid stage 40%, then IntaRNA (~3 h more); has checkpoints | `calib_bacterial_rnai/run_mirna/run_d.log` ends with `MIRNA_keepseed_DONE`, `MIRNA_ALL_DONE`; outputs in `run_mirna/keepseed/` |
| Wolbachia WsnRNA-46 vs *Aedes* CDS at 27 C | ~75%, ~30 min left | `calib_bacterial_rnai/run_wolb_cds/run.log` ends with `WOLB_CDS_DONE` and prints `validation.tsv` (rank of Dhc = `LOC5568788`) |
| agent `jobs` | following the two jobs, will append Results to `calib_bacterial_rnai/README.md` | `.scratch_solve/agents/jobs/REPORT.md` |
| agent `genelevel` | gene-level statistic experiment (merge only if it beats the site-count baseline within UTR-length decile by >= 0.03 AUC) | `.scratch_solve/agents/genelevel/DELIVERY/` |
| agent `speed` | **delivered and reported** (no REPORT.md file; the report is summarised in section 3, point 8): `.scratch_solve/agents/speed/DELIVERY/{speed.patch,docs_readme.patch,bench/}`; `patch -p1` from the repo root (score.py is CRLF), 125 tests pass in its copy (31 new in `tests/test_batch.py`). `rnahybrid_batch`/`intarna_batch` are ScoreConfig fields only, NOT wired into the CLI/config yet (do that at merge). Re-check the identity claim before trusting it | same folder |
| agent `bactmode` | **delivered**: `.scratch_solve/agents/bactmode/DELIVERY/bactmode.patch` (130 tests in its copy); its report is in the chat only, summarised in section 3 | apply with `patch -p1` on a copy |
| agent `predictors` | just started: per-set predictor policy, `support_from`, `called_by`/`best_q`, `prefilter_by`/keep-classes, `predictor_report.tsv`, README paragraph | `.scratch_solve/agents/predictors/DELIVERY/` |
| agent `transfect` | building a miRNA-transfection positive control WITH negatives | `calib_bacterial_rnai/transfect/` (has `expression.tsv`, `evaluate_repression.py`, `run_scoring.sh` already) |
| agents `verify`, `bactdata` | finished | `.scratch_solve/agents/verify/REPORT.md`, `calib_bacterial_rnai/bact_ecoli/SOURCES.md` |

If the usage limit cuts the agents off, their partial files stay where listed; nothing is lost
except unreported conclusions. Re-launch them from the prompts in this chat's history or redo the
task from `docs/improvement_plan.md` (the plan lists goal, deliverable, acceptance test per task).

## 5. What to do next (in this order)

0. **Make `make_decoys` deterministic** (it depends on `PYTHONHASHSEED`: sets are iterated, see section 3 point 7) and add a
   test that two processes give the same decoys; until then export `PYTHONHASHSEED=0` for any run you want to compare.
1. **Read results of the finished runs**: Wolbachia CDS (`run.log`), then miRNA keepseed (when done):
   `python calib_bacterial_rnai/scripts/analyze_scored.py <dir>`, `baseline_analysis.py <targets_scored> <truth> <regions>`,
   and the printed `validate-targets` summary. Decide the default decoy mode (shuffle vs keep_seed vs match_seed_gc)
   on the positive control and update `docs/module_b_development_log.md` (D6, section 3).
2. **Merge agent patches into a fresh branch from `module-b-step2`**, one at a time, running the full
   tests after each (`patch -p1 --dry-run` first):
   speed -> bactmode -> predictors -> genelevel (only if accepted). Known coupling: `bacterial.py` imports
   `rna`, `CSV_COLUMNS`, `parse_csv`, `run_tool` from `intarna.py` and `qvalues` from `fdr.py`; `speed`
   rewrites `rnahybrid.py/intarna.py` and the scoring loop of `score.py`; `predictors` touches the post-scoring
   part of `score.py`. Verify `speed`'s identity claim yourself (old vs new on 2,000 windows).
3. **Fix the discrepancies found by `verify`** (report in `.scratch_solve/agents/verify/REPORT.md`):
   `validate-targets` snippet lacks `--out` and the truth-file columns; `fdr_pseudocount` has no flag/config key;
   resume re-runs RNAcalibrate before checking the checkpoint (move the checkpoint check first); README says
   150 nt window (code 100), "first five" extended sites (four), `loci.min_count` default; provenance
   `conda_pins` collide on `python`; add missing config keys/columns to the README.
4. **Bacterial validation**: run the *E. coli* set (`bact_ecoli/`) through the merged `bacterial` mode, reduced
   (5-8 sRNAs x all start-codon windows, ~20 CPU-h; ask the user about CPU first); report the rank/percentile of the
   69 validated pairs with `validate-targets`, count pairs outside the window separately.
5. **Second positive control**: evaluate `transfect/` (repression vs score; baselines site count and UTR length).
6. **Decide and implement the product consequences of finding 2** (site q-values nearly never fire on true targets):
   gene-level statistic (if `genelevel` accepted), or document that ranking, not calling, is the supported use
   of seed sites; reconsider the meaning of `support`.
7. **Decisions the user must take** (ask): default decoy mode; whether to push `module-b-step2` / open a PR to `die-lab/SmithHunter2`
   (push: `git push git@github.com:die-lab/SmithHunter2.git module-b-step2`); keep or delete the 24 GB calibration data,
   simulations and `.scratch_solve`; whether the `bacterial` and `cross-kingdom` presets go into the default config;
   smithRNA seed (2-8 vs 4-10) default.
8. **Final documentation** from `docs/module_b_development_log.md` section 7 checklist (modes and presets with
   evidence, validation chapter including negative results, reproduction commands, clean install, CI green), then cleanup.

Not started by decision: precursor folding (step 3), report plots for module B, miRanda.

## 6. Cheat sheet

```bash
# environment for the user envs
export PATH=$HOME/miniforge3/envs/smithhunter2/bin:$HOME/miniforge3/envs/smithhunter2-rna/bin:$PATH
cd ~/THE_CAPTAIN/SmithHunter2; export PYTHONPATH=src
python -m smithhunter score-sites --help          # --dry-run counts the sites first
python -m smithhunter validate-targets --targets T --truth PAIRS --regions R --out OUT
python -m smithhunter decoys --decoy-mode shuffle|match_seed_gc|keep_seed ...
# positive-control files
calib_bacterial_rnai/mirna/control40/{candidates.fasta,truth.tsv}      # 40 miRNAs, 2,673 pairs
calib_bacterial_rnai/run_pa14/results/targets/human/regions.tsv        # human 3'UTRs (19,634 genes)
calib_bacterial_rnai/bact_ecoli/{candidates.fasta,truth.tsv,SOURCES.md}
calib_bacterial_rnai/truth/*.tsv                                       # PA14 sRNA52320, WsnRNA-46 pairs
```

## 7. Message to paste in a new chat

"Resume SmithHunter2 module B. Read ~/THE_CAPTAIN/SmithHunter2/docs/HANDOFF.md first, then
docs/module_b_development_log.md. Check the running jobs and agent deliveries in section 4, then
follow section 5 in order. Follow the server rules in section 0 (write only inside THE_CAPTAIN)."

## 8. Status update 2026-11-11 (later session; supersedes section 4/5 where they differ)

- Git rule changed by the user: commit and push to branches without asking, never to `main`
  (merge to main only when the user says so). Audit what goes out before every push: the remote
  is public, no server address, login or personal paths in the tree.
- Branches pushed to `die-lab/SmithHunter2`: `module-b-step2` (ceeb274, 50d05ac),
  `module-b-step3` (+ deterministic decoys, test over 4 hash seeds), `module-b-merge`
  (= step3 + speed patch + bactmode patch; 156 tests pass; local branch name `merge-speed`).
- Section 5 step 0 done. Step 2 done for `speed` (batch == per-site on 1,027 real windows, full
  dict equality) and `bactmode`. NOT done: `predictors` and `genelevel` (both agents died on the
  usage limit; DELIVERY folders are empty, partial work in `.scratch_solve/agents/<id>/`).
  `rnahybrid_batch` / `intarna_batch` are still not CLI/config keys.
- Finished jobs: Wolbachia CDS (WsnRNA-46 -> Dhc `LOC5568788`: IntaRNA percentile 0.075, RNAhybrid
  0.124, AUC 0.93 / 0.88, support 0, best q 0.64, one validated pair) and miRNA keep_seed (AUC
  RNAhybrid 0.7908, IntaRNA 0.557 on the 15% prefiltered, still only 2 supported sites): the decoy
  mode does not rescue site-level power, so finding 2 stands. Not yet written into the log (D6, section 3).
- Tests need a SHORT `TMPDIR` (e.g. `~/THE_CAPTAIN/tt`), otherwise AF_UNIX forkserver paths are too long.
- Next: step 1 write-up in the log, step 3 (verify discrepancies), re-run `predictors` and `genelevel`,
  then step 4 (ask the user about CPU first).
