# Module B development log: the logic, the evidence, the open points

Living document. It records *why* module B is built the way it is, what has been checked and
how, what failed, and what is still open, so that the final documentation can be written from
it and every claim can be traced to a test, a run or a source. Started 2026-11-10 (step 2).
Reproduction details for each experiment are in `~/THE_CAPTAIN/calib_bacterial_rnai/README.md`.

## 0. Principles

1. **The focus genome is whatever the study is about.** A mitochondrion, a bacterium, an
   endosymbiont (Wolbachia), a plastid... Roles (`focus` / `exclude`) are declared per genome;
   nothing in module B assumes a mitochondrial origin. Target sets (what the small RNAs are
   tested against) are declared separately from the genomes that make the small RNAs.
2. **Keep every score; trust is measured, not assumed.** No fixed cut-off. Each predictor's
   score is kept; q-values come from decoy small RNAs that go through the identical pipeline.
3. **The site comes first** ('ago' mode): predictions are anchored on canonical seed sites; the
   energy predictors score sites, they do not search the transcriptome freely.
4. **Say what a number means.** The FDR here is an empirical decoy estimate relative to a
   composition-matched random null; it is not a probability that a site is functional.
5. **Negative results are results.** A positive control that cannot be recovered by the model
   is reported as such, not tuned away.

## 1. Data flow of step 2

```
candidates (module A or published)  ->  decoys (N per candidate, a decoy mode)
target regions (genome+GFF3 or transcripts) -> unique sequences (isoforms scanned once)
scan(sites)  ->  RNAhybrid on every site (window +-15 nt, -f = seed range, -d from RNAcalibrate)
             ->  [optional prefilter: best X% per class by MFE, one pooled cutoff, real+decoy]
             ->  IntaRNA on the survivors (window +-100 nt, seed pinned, interaction confined)
             ->  per predictor and class: q-values (pooled and per candidate), support
checkpoints after stage 1 and at the end (resumable); tables + FDR curves
validate-targets: rank / percentile / AUC of validated pairs among all genes
```

## 2. Decisions, reasons and evidence

| # | decision | reason | evidence / status |
|---|---|---|---|
| D1 | seed classes at guide 2-8 (8mer, 7mer-m8, 7mer-A1, 6mer, offset-6mer) | canonical miRNA seed rules | literature SUPPORTED (Agarwal 2015; Bartel-lab papers); unit tests on hand-made targets |
| D2 | IntaRNA per site, seed pinned (`--seedQRange`, `--seedTRange`), `--tRegion` = guide reach, `seedBP` per class | energy of the interaction anchored on the site, not the best duplex anywhere | PARTLY supported: no IntaRNA benchmark on miRNA data exists in what we found (Mann 2017; Raden 2020 are bacterial); tested on planted sites with the real tool |
| D3 | IntaRNA window = site +-100 nt | accessibility in context at a bounded cost (0.14 s/site at 150 nt, 0.05 s at 60 nt) | same scale as IntaRNA defaults (tAccW 150 / tAccL 100); sensitivity to the flank NOT measured |
| D4 | RNAhybrid per site with `-f` = seed range | forces the seed helix | Rehmsmeier 2004 uses the same constraint; the duplex must also cover the seed of THIS site (`rnahybrid_require_site`) |
| D5 | p-values from RNAcalibrate on random sequences (default 500+-300 nt, target dinucleotide frequencies, same `-f`) | the paper advises per-query calibration with the same constraint | informational only; q-values use MFE / energy. Calibrating at the real window length (~55 nt) was tried and FAILS (xi < 0 below ~100 nt, probed with RNAcalibrate on the server), so the extrapolation to short windows through RNAhybrid's length normalisation is kept and stated |
| D6 | decoys: dinucleotide shuffles of the whole small RNA (`shuffle`); alternatives `match_seed_gc`, `keep_seed` | standard null in RNA work (Clote 2005) | CONTESTED: decoy seeds differ from the real seed (GC confounder). On the Pseudomonas panel `shuffle` and `keep_seed` agreed. Not yet compared on a positive control (miRNA run in progress) |
| D7 | q-value = min FDR over thresholds, decoy sites weighted 1/N of their candidate, +1 pseudo decoy, per predictor and class, pooled and per candidate | target-decoy estimate (Elias & Gygi 2007; Storey & Tibshirani 2003; Barber & Candes 2015) | SUPPORTED; conservative on average, no proof of control; candidates without decoys get no q |
| D8 | `support` = number of predictors with pooled q <= FDR | simple, transparent | a count of correlated votes, NOT a confidence (Spearman IntaRNA/RNAhybrid 0.3-0.4) |
| D9 | two-stage mode: RNAhybrid on all sites, IntaRNA on the best share per class | 98 CPU-hours vs ~19 for 3.5M sites | valid with one pooled label-blind cutoff; IntaRNA q is conditional on passing; recall on the simulation 4/4 at 15% and 30%; recall on real data NOT measured |
| D10 | checkpoints, timeouts, per-site error tolerance (stop above 5%) | hours-long shared-server runs | tested (fakes + end-to-end identical outputs on resume) |
| D11 | strong classes only by default (`score_weak_sites: false`) | 6mer / offset-6mer are ~10x more frequent | cost decision, not a biological one |
| D12 | `bacterial` mode: IntaRNA over the whole start-codon window (-200..+100), seedBP 7 at any position, local accessibility (tAccW 150, tAccL 100), no RNAhybrid and no prefilter | trusted sRNA-mRNA pairs lie in -150..+100 (101 of 102); the IntaRNA 2.0 benchmark used -200..+100 and seedBP 7; RNAhybrid alone is much worse than IntaRNA for bacteria (section 6) | literature, second-hand (section 6); implemented 2026-11-10; behaviour checked on synthetic data only (section 8), no real sRNA positive control yet |
| D13 | unit of the FDR in bacterial mode = best interaction per (small RNA, window), not per seed match | the miRNA positive control showed that the unit matters (section 3b): many site-level hits per gene gave no power. One interaction per window keeps the number of tries equal to the number of genes for real and decoy alike | by construction; power on planted pairs in section 8 |
| D14 | windows of neighbouring operon genes overlap and may contain the end of the upstream CDS: kept, `upstream_cds` reports the overlap | sRNAs also pair inside upstream coding sequence; the 5'UTR is not annotated in most genomes; removing the overlap would hide real pairs | design decision; not evaluated on real data |
| D15 | circular contigs: the window continues across the origin (`topology: circular` / `circular_contigs`); linear contigs: cut, the start codon keeps its place | otherwise the windows near the origin of a circular genome are wrong | unit tests on both strands (genome mapping checked base by base); end-to-end test with a planted site straddling the origin |
| D16 | IntaRNA energy only: `intarna_q` (pooled) and `intarna_q_cand`; no `support` | one predictor; same weighted target-decoy logic as D7 (decoy weights, +1) | unit tests with a faked IntaRNA; real tool in section 8 |

## 3. Experiments (what was run, what came out)

| run | what | result | where |
|---|---|---|---|
| simulation, CI parameters | planted sites: extended pairing (4 usable), bare seed, 170 others | extended supported 4/4, bare seed 0/11, other 0/170 | `sim_ci_clean/`; asserted in CI |
| Pseudomonas PA14 OMV vs human 3'UTRs, shuffle decoys | 30 candidates incl. the published sRNA52320 | 0/5 "validated" genes have a scored site; one supported candidate (paralog, RNAhybrid MFE -59..-61) | `calib_bacterial_rnai/run_pa14/results/targets/human` |
| same, keep_seed decoys + duplex-on-site rule | same candidates | 14 vs 18 supported sites, same candidate | `.../human_keepseed` |
| Wolbachia WsnRNA-46 vs Aedes 3'UTRs | one published sRNA, 20 decoys | no strong site in the Dhc 3'UTR (220 nt). SUPERSEDED: the literature check shows the paper's Dhc site is in the CDS (RNAhybrid MFE -30.3), so a 3'UTR-only run could not test it | `run_wolb/` |
| Wolbachia WsnRNA-46 vs Aedes CDS, 27 C (the paper's cell temperature), 40 decoys | the same sRNA, coding sequences | IN PROGRESS | `run_wolb_cds/` |
| Wolbachia small RNAs, module A (focus = wMel, exclude = Aedes), 4 libraries | candidates for the bacterial side | works as designed: wMel reads 1.10M (Pop cytoplasm) and 1.85M (Pop nucleus) vs 206 and 370 in the uninfected controls; 41,213 loci, 570 pass (433 CDS fragments, 69 intergenic, 40 tRNA, 20 rRNA...). **The published WsnRNA-46 (27 nt) is nearly absent from the deposited reads**: 1 exact read in Pop-cytoplasm, 3 reads with a different 5' extension (TTTTCTGGAT...) in Pop-nucleus, 0 in controls, so it fails the expression filter (4.25 reads) | `run_wolb_A/` |
| 40 human miRNAs, 2,673 miRTarBase reporter pairs (11.0), shuffle decoys | positive control, 19,634 genes | 2,157/2,673 pairs have a strong seed site in the 3'UTR. Gene-level AUC of the best RNAhybrid MFE 0.79, top 5%: 20.7%, top 10%: 37.9% of pairs. But only 2 of 1.7M sites get `support`. Baselines (same tie handling): 3'UTR length 0.644, **number of strong sites 0.797**, best MFE within its length decile 0.747. IntaRNA values exist for 15% of sites only (prefilter), so its percentiles are not comparable | `run_mirna/shuffle/`, `scripts/baseline_analysis.py` |
| same, keep_seed decoys | the decoy-mode comparison on a positive control | IN PROGRESS (8.2M sites; saves checkpoints incl. decoy scores) | `run_mirna/keepseed/` |

## 3b. What the positive control says (2026-11-10)

1. A strong canonical seed site exists in the 3'UTR for 81% of the reporter-validated pairs; the
   energy (RNAhybrid MFE) ranks validated genes well above chance at gene level (AUC 0.79).
2. **That ranking is mostly the number of seed sites**: counting strong sites alone gives 0.797,
   UTR length alone 0.644, and the MFE within a length decile 0.747. The energy carries signal
   beyond length but, at gene level, none beyond the site count that we could show.
3. **The site-level decoy q-values have almost no power on known true targets** (2 supported
   sites in 1.7M): a canonical seed match scores about as well as the seed matches of shuffled
   decoys, so "better than random-seed sites" does not separate it. The statistics are valid,
   the unit (one site) is not where the signal is.
4. Reporter-assay positives in miRTarBase are biased toward long, well-studied 3'UTRs and include
   non-canonical sites and CLASH rows; the set has no negatives. The AUCs are a ranking sanity
   check, not a performance estimate.
5. Proposed next step: a **gene-level statistic** (aggregate of the site scores of a gene) tested
   against the same decoys, so that UTR length and site number are in the null. To be evaluated
   offline on the checkpointed scores of the keep_seed run before any code is written; if it does
   not beat the site-count baseline within a length decile it is not worth keeping.

## 4. Corrections made along the way (so they are not repeated)

- `rnahybrid_at_site` first compared the duplex *start* with a fully paired guide: wrong for
  seed-only sites (25% flagged). Now: the duplex must cover the seed.
- `t1` of non-canonical seeds (`--seed 4-10`) was off by one, so IntaRNA scored nothing silently.
- Candidates without decoys got q = 0; unequal decoy counts were not weighted.
- The calibration was first done at region length, then at window length on a reviewer's
  suggestion; the second fails (xi < 0 for ~55 nt), so it is back at RNAcalibrate's defaults.
  Docs had said "shuffled copies" (RNAcalibrate samples random sequences from the dinucleotide
  distribution of the target set).
- The Wolbachia/Dhc case was called "uninformative" from a 3'UTR-only run; the paper's site is in
  the CDS, so the verdict was premature. The target region is part of the hypothesis.
- The default temperature (37 C) was silent; it now warns when not set.
- My first baseline put unscored genes strictly last while `validate-targets` gives them the mean
  tail rank (0.715 vs 0.79 for the same score). Both now use the tie rule; they agree to 3 decimals.
- Worker CPU was measured at the wrong process level (Python 3.14 uses forkserver).
- Scratch files went to `/tmp` once (deleted); everything since stays in THE_CAPTAIN.

## 5. Open questions (to close before the tool is called finished)

1. A **seed-based positive control with real numbers** (miRNA run) and the decoy-mode comparison on it.
2. **Bacterial side:** the Wolbachia small-RNA libraries through module A (focus = Wolbachia,
   exclude = Aedes), then candidates against Aedes 3'UTRs; and the `bacterial` mode (now implemented, section 8)
   (start-codon windows, no fixed seed) for bacterium -> bacterium targets.
3. **Which mechanism applies to which genome** (AGO-guided seed rules for mitochondrial and for
   bacterial small RNAs in a host vs Hfq-style pairing inside the bacterium): a literature
   comparison is in progress; its recommendations go into section 6.
4. Sensitivity to the IntaRNA flank and to temperature; recall of the two-stage mode on real data.
5. miRanda as an optional third predictor; weak sites; per-organism temperature.
6. Precursor folding (step 3) and the report plots (not started, by decision).

## 6. Parameter presets per interaction mode

Source: a literature comparison by a subagent (2026-11-10; WebFetch summaries are second-hand,
so treat as leads to confirm in the primary text). Evidence: V = read in a fetched source,
R = recalled, I = inferred. Where the literature is thin, the preset is marked as a default
to be tested, not as a finding.

| parameter | `ago` (mitochondrial smithRNAs, miRNA-like) | `bacterial` (sRNA -> bacterial mRNA) | `cross-kingdom` (bacterial sRNA -> host) |
|---|---|---|---|
| machinery | AGO binding of smithRNAs is only weakly supported (re-analysis of public AGO2 CLIP, Pozzi & Dowling 2022; no chimeras) (V) | Hfq / RNase E / RNase III, no AGO, no guide strand (R) | host AGO only if tested; Legionella: AGO2 inhibitor reduced the effect, AGO2-IP not significant; Pseudomonas, Wolbachia: AGO not tested (V) |
| seed | canonical 2-8 (Plazzi 2024); Pozzi 2017 / SmithHunter use a shifted 4-10 (V): keep 2-8 default, `seed: [4,10]` as an exploratory option | no fixed position, ~6-12 nt, GU common (R); IntaRNA seedBP 6-7 over a window (V) | inherit the host rule; Wolbachia used 2-8 with seed-mutant reporters (V) |
| target region | 3'UTR primary (Pozzi 2017, Plazzi 2024); 5'UTR/CDS reported apart, never in one q-value | -200..+100 around the start codon (V: 101 of 102 trusted pairs in -150..+100) | as the study proposes: Wolbachia Dhc site is in the CDS (V), others 3'UTR |
| predictors | RNAhybrid + IntaRNA; 3'UTR accessibility | IntaRNA with local accessibility does far better than RNAhybrid alone (Pain 2015, scores as reported by the subagent: CopraRNA 46, IntaRNA 27, RNAplex 21, RNAhybrid 4; scale not checked) (V); do NOT prefilter with RNAhybrid; CopraRNA optional | RNAhybrid + IntaRNA, then ask whether the host AGO is plausible |
| temperature | organism-specific (clams 12-22 C in the papers, computed at 37 C "for software limits") (V): do not default to 37 | 37 C or the growth temperature of the study | the host cell temperature (Aag2 27 C, THP-1 37 C) (V) |
| decoys | dinucleotide shuffles of the small RNA; keep length and 5' base; compare `shuffle` / `match_seed_gc` / `keep_seed` | shuffled sRNA and shuffled mRNA windows | as `ago` |
| do not assume | AGO loading; 3'UTR-only targeting; a 5' U | a seed position; 3'UTR targeting | that the bacterial RNA behaves as a miRNA |

Consequences for the code (status):
- done: `seed` configurable (non-canonical classes work, `t1` bug fixed); `temperature` warns when unset;
  `regions` accepts CDS / 5'UTR / 3'UTR per target set.
- done (2026-11-10): `bacterial` mode (start-codon windows, IntaRNA only, no RNAhybrid prefilter; section 8).
- to do: keep one q-value per region type; a `cross-kingdom` preset in the config; make the unit of "unvalidated mechanism" explicit in the output.
- noted: the smithRNA papers quote RNAhybrid thresholds in kJ although RNAhybrid reports kcal/mol (probable labelling slip, I).

## 7. Checklist for the final documentation

- [ ] README: modes, parameters, outputs, limits (statistical notes section exists)
- [ ] Per-organism parameter presets (section 6) with the evidence for each
- [ ] Validation chapter: simulation, miRNA positive control, bacterial/cross-kingdom cases, negatives included
- [ ] Reproduction commands for every number quoted
- [x] Bacterial mode implemented and documented in the README (2026-11-10); a real positive control is still open
- [ ] Environment created and `pytest` + CI green on a clean machine
- [ ] Cleanup: `.scratch_solve/`, simulation and scratch directories; git commit and push (not done: nothing is committed)

## 8. Bacterial mode: synthetic checks (2026-11-10)

Code: `targets/bacterial.py` (scoring), `targets/regions.py` (`bacterial_windows`). Tests:
`tests/test_bacterial.py` (windows on both strands, clipping, circular wrap, GFF flavours, q-values with
a faked IntaRNA, checkpoint/resume, one test with the real IntaRNA and one end-to-end run).
These are checks of the machinery on planted data, not a validation on real sRNA-mRNA pairs.

**Power by length of the planted pairing.** One sRNA (22 nt), 200 random windows of 300 nt (GC
0.4-0.6, start codon at 201), 5 planted windows (random, GC 0.5, complementary to the 5' 10/12/14/16/18 nt
of the sRNA, ending at -9), 20 decoys (dinucleotide shuffles), IntaRNA seedBP 7, tAccW 150, 37 C:
4,305 IntaRNA runs. IntaRNA found no interaction in 601 of the 4,305 pairs (14%; the seed rule
excludes them, for real and decoy alike). Best E of the 162 random windows with an interaction:
median -5.9, minimum -14.8. Decoys (3,449 pairs with an interaction): median -6.2, minimum
-18.7, 30 pairs <= -14, 114 <= -12.

| planted pairing | E (kcal/mol) | hit | q alone among the 200 random windows |
|---|---|---|---|
| 10 bp | -10.0 | 5UTR, RBS | 1.00 |
| 12 bp | -16.6 | 5UTR, RBS | 0.25 |
| 14 bp | -21.2 | 5UTR, RBS | 0.05 |
| 16 bp | -25.1 | 5UTR, RBS | 0.05 |
| 18 bp | -25.2 | 5UTR, RBS | 0.05 |

No random window passed q <= 0.05 (0 of 200). q = 0.05 is the floor with 20 decoys and one passing
window (the +1 pseudo decoy), so "supported" at an FDR of 0.05 means that no decoy window in 20 x 200
tries scored better. **A perfect 12-bp pairing is not separable from chance** at this scale (the
best chance interactions of shuffled decoys reach -14..-19); 14 bp or more is. Real sRNA-mRNA
duplexes are usually longer than 12 bp but often imperfect, so the sensitivity for short or
interrupted pairing is limited by the null, not by the code; more decoys lower the floor of q but
not the ratio. The task asked for a pairing of about 12 bp; the end-to-end test uses 16 bp for this reason.

**End to end.** A synthetic circular genome (30 genes on both strands, GFF3, two windows crossing the
origin), one candidate with a 16-bp target planted at -34..-19 of the last gene (minus strand, the
target straddling the origin), one random candidate, 20 decoys each, run through `regions --mode
bacterial`, `decoys` and `score-windows`: the planted window is the top row (E -22.6, q 0.05, `5UTR`,
`hit_rbs` 1), the reported genomic span maps back to the planted positions in the genome (start > end,
across the origin), the random candidate has no window with q <= 0.05. A 12-gene version went through
the Snakemake rules (`targets`) with the same tables.

Cost measured: 0.15-0.2 s per IntaRNA run (300-nt window, 22-nt query) on the loaded shared server;
batching several windows into one call did not make a run faster (the time is the dynamic
programming, not the start-up).

Open points: no real positive control (a published sRNA with known targets in a genome with a GFF3
would test the window and the seed settings); seedBP 6 vs 7 and `no_seed` not compared; sensitivity
to the window size and to the temperature not measured; the query is a fragment when only reads are
available (section 6); decoys are shuffles of the query only (shuffled mRNA windows are not used as a
second null); windows with `upstream_cds` > 0 are not treated separately.
