# Bacterial internal mode (sRNA and its own host's mRNAs): parked results

Status: alternative, additional branch. The main target of the tool is host-parasite RNAi; this branch keeps the
bacterial-internal evaluation (E. coli K-12) for a possible internal-detection tool. Scripts: `analysis/bacterial_internal/`
(they expect the run directory layout of the E. coli run: windows_scored.tsv, regions.tsv, work/windows_checkpoint, truth tables).

## Setup
8 sRNAs (Spot42, RyhB, GcvB, OmrA, OmrB, FnrS, CyaR, RprA) against 4,285 start-codon windows (-200/+100), IntaRNA with
accessibility and seed 7, 37 C, 20 shuffle decoys per sRNA. 36 validated pairs (tier 1).

## Results
* Ranking: median percentile of validated pairs 4.1%, AUC about 0.89; 25/36 in the top 10%, 9 in the top 1%.
* Temperature (27 C, 15 C): ranking practically unchanged.
* Pooled shuffled-sRNA FDR (hit counts): no validated pair passes. It measures sRNA stickiness, not specificity.
* Per-sRNA null (Gumbel on each sRNA's own energies, BH over 34,185 tests): no window passes q <= 0.25.
* Paired window test (real energy vs the 20 decoys at the SAME window, z with normal approximation, BH): 16/36 validated
  pairs beat all 20 decoys (3.6% of all genes do); 312 windows pass q <= 0.05 (5 validated, 12 curated known). Ranking is not
  better than raw energy (AUC 0.856 vs 0.876). The two strongest calls overlap the sRNA's own locus (to be labelled, not removed).
* Truth tiers: curated known targets (Melamed 2016 table S3) AUC 0.85; RIL-seq pairs in >= 2 libraries AUC 0.67, only 4% of
  them have the chimeric fragment inside the window. RIL-seq is a soft positive set only.
* Truth set is incomplete (GcvB 5 of about 54 known targets).

## Open items
Fold-stability z-score of real sRNAs against many shuffles; paired test with 50-100 decoys per window; label (not remove)
windows overlapping the sRNA's own locus; conservation scoring; contiguous >= 16 nt score.
