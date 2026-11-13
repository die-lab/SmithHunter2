# Module B – revision after the literature review (2026-11-13)

This document revises `module_b_design.md` in the light of a three-strand literature
review (mechanism; validated pairs; prediction tools). The full reviews are in
`research/A_mechanisms.md`, `research/B_validated_pairs.md`, `research/C_tools_benchmarks.md`.

It does not replace the design. Most of the design survives: the architecture that was
chosen for taste reasons in October ("every score is kept", "trust is measured",
"consensus is a column, not a filter") turns out to be the only architecture the
evidence supports. What changes is the framing, several defaults, and the claim the
tool is allowed to make.

---

## 1. The three findings that change the design

### 1.1 A seed match carries no information. This is not a bug in our FDR.

Step 1 gave a site-count FDR of ~1 and that was read as "ranking needs step 2". The
literature says it is the central result of the field:

- Qin 2023: ~10⁶ plant sRNAs could target **97% of all fungal genes** by complementarity.
- Friedman's signal:background for endogenous miRNA sites is **2.6 for conserved 8mers,
  1.2 for 6mers** — a third of conserved 8mers are chance, in the best-characterised
  system that exists.
- Qin 2022: 3,185 predicted tomato targets, 163 responsive, **none** surviving genetic
  test; `dcl1 dcl2` mutants fully virulent.

Our FDR of ~1 is the correct measurement. The published cross-kingdom work that reports
thousands of confident targets simply never computed one. **That is the contribution of
this tool, and it should be the headline rather than a caveat.**

### 1.2 There is no benchmark set, and there cannot be a ROC curve.

49 named pairs across the whole cross-kingdom literature; 4 at grade A, 17 B. Every
single A/B pair was discovered by running psRNATarget / TAPIR / miRanda and validating
the top hits. **The positive set is the output of the tool class we would be
benchmarking.** Published non-targets: one.

Consequences:

- No AUC, no ROC, no precision/recall against literature pairs. Any such number would be
  circular and we should refuse to produce one.
- The defensible metrics are **median rank of the true target** (needs no negatives — the
  standard in the bacterial sRNA field, cf. Pain 2015) and **calibration of the
  decoy-derived q-values on simulation**, which we already have.
- The simulator becomes the primary validation instrument, not a secondary one. It is
  the only place where ground truth is actually known.

### 1.3 For *Wolbachia* specifically, the mechanism is unestablished.

- **No bacterial sRNA has ever been shown by AGO-IP or CLIP to enter a host Argonaute.**
  Not one, in any system.
- No degradome or 5'-RACE evidence for any bacterial→animal pair.
- Three *Wolbachia* sRNA papers exist in total. One reporter-level host pair
  (WsnRNA-46 → *Dhc*), and its effect is **up**-regulation.
- Mayoral 2014, the foundational paper, used BLAST + RNAhybrid + RNA22 with no reported
  cutoffs, no FDR and no decoys.
- The better-supported mechanism in *Wolbachia* is modulation of **host** miRNAs — a
  confounder almost never controlled.
- ≥60% of named bacterial "effectors" in the literature are tRF/rRF/tmRNA fragments: the
  same species that dominate degradation and contamination.

So the `ago` mode applied to a bacterial endosymbiont is an **extrapolation from
eukaryote-to-eukaryote systems**, and the tool must say so in its own output. Realistic
prior for a sequence-only bacteria→animal prediction: **≲1%**.

---

## 2. What this means for the tool's claim

Replace "target prediction" with **calibrated hypothesis generation**. Concretely, module
B answers: *given this candidate, how much complementarity signal is there against a
composition-matched null, and how does that compare to the chance expectation for this
search space?* It does not answer *which genes does this sRNA regulate*.

This is not modesty. It is the difference between this tool and the
3,185-predicted-targets-no-phenotype literature, and it is defensible at review in a way
that a ranked target list is not.

**Action:** the output header of `targets.tsv` and the report must state the computed base
rate for the actual run: `search space L nt; expected chance 7mer matches ≈ L/16384 = N`.
A reader seeing N = 73 next to 91 reported hits calibrates correctly with no further text.

---

## 3. Revisions to the design, item by item

### 3.1 Keep as designed (now with evidence behind it)

| Design decision | Status |
|---|---|
| Every score kept, consensus a column not a filter | **Confirmed.** Seed and thermodynamic scoring fail in opposite directions; their agreement *is* the information. Never sum them. |
| Empirical FDR from decoys through the identical pipeline | **Confirmed, and it is the novel contribution.** The field's own methods review recommends no randomised controls at all. |
| Precursor folding reported beside targets, never a filter | **Confirmed.** tRNA z = −1.59 (p ≈ 0.12), mRNA ≈ 0. MFE-based hairpin criteria are known high-FP. Annotation only. |
| Seed 2–8 as default, site classes, best-class-only assignment | **Keep** for `ago` mode, with the caveat of §3.3. |
| Bacterial mode: −200/+100 around start codons, no seed rule | **Confirmed.** Pain 2015: replacing this window with real RNA-seq UTRs changes AUC by 4.9–17.5%, none significant. Do not spend effort on exact UTR boundaries. TargetRNA3 SHAP: seed length contributes little in bacteria. |
| Decoys preserve dinucleotide composition | **Confirmed.** Mononucleotide shuffling makes real sequences look spuriously significant. |

### 3.2 Change: search region is 3'UTR **+ CDS**, not 3'UTR only

Decision 3 of the October list (3'UTR only, CDS opt-in) should be reversed.

- CLASH puts ~50% of miRNA sites in CDS.
- Every validated **plant**-host cross-kingdom pair is CDS, near-complementary, cleavage-type.
- The 3'UTR consensus in the cross-kingdom literature is a **methodological artefact of
  the tools used**, not a measured enrichment. Nobody has tested 3'UTR vs CDS vs 5'UTR on
  equal footing.

Restricting to 3'UTR would make a true effector invisible by construction. Report the
region as a column and let the analysis decide.

Also: **normalise for 3'UTR length.** Predicted targeting multiplicity scales with UTR
length at r = 0.63. Without normalisation the output is a list of long-UTR genes.

### 3.3 Change: three modes, not two

The plant/animal split is mechanistically real and both reviews found it independently:

| mode | pairing model | justification |
|---|---|---|
| `ago_animal` (current `ago`) | contiguous WC 7-mer at guide 2–8; 6mers unscored | every validated animal-host pair is seed-driven miRNA **mimicry** |
| `ago_plant` (new) | psRNATarget-style penalty scoring over the full guide, ≤4–5 mismatch-equivalents, positions 2–13 weighted, **pairing required across 10–11** | every validated plant-host pair is cleavage-type |
| `bacterial` | no seed rule, IntaRNA over start-codon windows | unchanged |

`ago_plant` is low priority for Zyginidia but costs little and makes the tool general, which
was the stated goal. Benchmarking the two `ago` modes together would average two
incompatible mechanisms.

### 3.4 Change: A1 anchor becomes a flag, not a class requirement

The A1 anchor is a property of the AGO MID pocket, not of base pairing. For a bacterium
with no established AGO loading, enforcing it assumes the conclusion. Emit the
A1-agnostic call **and** a boolean `a1_anchor_present`.

### 3.5 Change: drop the TargetScan context features

Section 4.1 lists AU content, distance from region ends, and position within the UTR.
**Do not implement these as scores.** They are quantitative AGO/UTR-biology constants
fitted in human, mouse and fly cells. Importing their coefficients into a
bacterium-to-insect setting manufactures false confidence.

**Keep** 3' supplementary pairing at guide 13–16 as a reported integer (`supp_bp_13_16`) —
it has a generic base-pairing rationale and was the most repressive CLASH class — but
**do not weight it into any score.**

Note: McGeary 2019 found flanking dinucleotides move 8mer affinity ~100-fold, often more
than site type does. We cannot use that model here, and that is itself worth stating: the
single largest known determinant of site efficacy is unavailable to us.

### 3.6 Change: two nulls, not one

Current design: 20 dinucleotide-shuffled decoys per candidate. Revise to:

- **Null A (primary):** 500 Altschul–Erickson dinucleotide-preserving shuffles per
  candidate, scored against the **unshuffled** target set. 500 rather than 20 because
  q ≈ 0.01 needs decoy mass in the tail. Shuffle the **whole** sRNA including positions
  2–8 — holding the seed fixed collapses the null for the feature under test.
- **Null B (orthogonal):** ≥200 real sRNAs of matched length from a phylogenetically and
  ecologically distant bacterium with no shared host. Controls for everything a shuffle
  destroys: real folding, real k-mer vocabulary, real 5'-end bias. Document the donor.
- **Never shuffle the targets.** It alters the search space rather than the query,
  destroys real UTR/CDS composition structure, and gives an anti-conservative FDR.

Report both q-values. **Disagreement between them is a demotion signal** and should be a
column.

Cost: 500 decoys × 30 candidates is 25× the current plan. Mitigate by keeping 6mer and
offset-6mer out of the energy predictors by default (already decided) and by scattering in
Snakemake (already designed). If cost is prohibitive, 200 decoys resolves q to ~10⁻⁵,
which is enough.

### 3.7 Change: stratify the q-values further

Current: per site class, pooled over candidates. Add **transcript-length decile × GC
quartile**. Long GC-rich transcripts win on energy for reasons that have nothing to do
with biology, and a pooled q hides that the 6mer stratum has ~293 expected chance hits per
sRNA in a bacterial window versus ~18 for 8mer.

### 3.8 New: hard filters before any pairing score

These matter more than the scoring function. None are currently in the design.

1. **Contamination filter.** Reject candidates with exact/near-exact identity to host
   rRNA, tRNA or abundant host ncRNA (the peu-MIR2910 lesson). Flag candidates whose
   abundance correlates with batch rather than biology. Where library metadata permits,
   model index hopping at 0.003–0.009 per read. Context: >7% of public miRNA datasets are
   cross-contaminated; xenomiRs in 824 datasets sit at 0.001% and 81% are rodent.
2. **Alternative-hypothesis flag.** Any tRF / rRF / tmRNA / Y-RNA fragment carries a flag
   that its host effect may be **TLR7/TLR3 innate sensing** rather than base-pairing. This
   is the majority class among published bacterial "effectors", so it is not an edge case.
   Flag, do not exclude — selective export is the discriminating evidence, not sequence class.
3. **Confounder flag.** For symbiont target sets, state in the report that an observed host
   phenotype may be mediated by symbiont-induced changes in **host** miRNAs. This is the
   better-supported *Wolbachia* mechanism and the control that is almost always missing.
4. **Export/enrichment field (optional input).** If the user can supply a
   vesicular/host-cytoplasmic fraction, enrichment over whole-symbiont RNA is the single
   strongest available prior. Accept it as an optional per-candidate column.

Module A already annotates tRNA/rRNA overlap via `annotation_class`; items 2 and 3 are
mostly reporting work on top of what exists.

### 3.9 New: bidirectional effects

Two catalogued pairs — including the only *Wolbachia* one — report **up**-regulation. A
silencing-only output model is known to be incomplete. At minimum, do not word the output
as "repression"; "predicted interaction" throughout.

---

## 4. Concrete parameters

### 4.1 IntaRNA (bacterial mode)

From the only systematic parameter sweep (Raden 2020), which recovered ~30% more verified
interactions at ~20% of baseline runtime:

```
--seedBP 7
--seedNoGU
--seedMaxE -2
--seedMinPu 0.001
--intLenMax 60
--intLoopMax 8
--outMaxE -10
--qAcc C --qAccW 0 --qAccL 0        # fold the whole sRNA: its structure is real
--tAcc C --tAccW 150 --tAccL 100    # local window on target; global folding of long mRNA
                                    # is slow and unreliable
--mode H
--outMode C --outNumber 1 --outOverlap N   # one best interaction per transcript, for FDR
```

Optional second run with the **Andronescu-2007** energy model (+5 verified pairs over
Turner-2004 in the same sweep); report hits found by both parameterisations.

Opt-in refinement: restrict `--qRegion` to predicted single-stranded regions of the sRNA
(Beisel's Spot 42 result). Cheap, and the closest thing the bacterial field has to a free
lunch after parameter tuning.

### 4.2 RNAhybrid

- `-f 2,8` in seed mode (this is how a seed is pinned in RNAhybrid — fixes the v0 bug).
- **`-n` must be raised above the sRNA length: the default is 30 and long sRNAs fail
  silently.** This is a real trap for bacterial-mode queries.
- `-e` with a generous `-b` so the hit cap does not truncate.
- `-d ξ,θ` from `RNAcalibrate` on **our own** target set. Never `-s 3utr_human/fly/worm`.
  Already in the design; the review confirms it matters.

### 4.3 Seed pinning is mandatory in `ago` mode

Scanning all sRNA windows multiplies the chance-hit base rate by (L−6): ~94× for a 100-nt
sRNA. Keep an `--all-windows` escape hatch, **off by default**, with a loud log warning
stating the expected chance-hit count.

### 4.4 Base-rate arithmetic to print

- Pinned 7mer over a 1.2 Mnt bacterial window: ~73 chance hits.
- Unpinned over a 100-nt sRNA: ~6,900.

---

## 5. Revised implementation order

Supersedes §10 of the design. Ordered by evidence-weighted value, not by code dependency.

1. **Base-rate line in the output header, and seed pinning enforced.** Two orders of
   magnitude of false positives removed by a default. Nearly free.
2. **Tune the IntaRNA parameters (§4.1).** Largest single performance gain available;
   costs a config file.
3. **Null A: 500 dinucleotide-shuffle q-values, stratified (§3.6, §3.7).** The core
   scientific contribution over existing cross-kingdom practice.
4. **RNAcalibrate wiring for RNAhybrid (§4.2).** Small, and it makes the second opinion
   mean something.
5. **CDS added to the default search region + `supp_bp_13_16` (§3.2, §3.5).** Cheap.
6. **Hard filters and flags (§3.8).** Mostly reporting on top of module A's annotation.
7. **Null B: unrelated-organism decoys (§3.6).** The orthogonal check.
8. **Simulation extension** — planted sites of known class, calibration check that
   unplanted candidates give ≈ `fdr` of their calls. This is now the *primary* validation.
9. **Beisel single-stranded `--qRegion`** (opt-in).
10. **Folding annotation fields.** Last, explicitly non-filtering.
11. **`ago_plant` mode (§3.3).**

**Deliberately not on this list:** TargetScan context++ features (no mechanistic warrant
here); a trained ML model (the entire bacterial field has ~4,386 pairs, none cross-kingdom);
conservation/comparative scoring.

On the last: comparative genomics is **the** signal that makes CopraRNA the best tool in its
class — median best-target rank 4, versus IntaRNA 27, RNAup 22.8, RNAplex 20.8, TargetRNA2
18.5, RNAhybrid 4.1 (Pain 2015); top-10% recall in *E. coli* 74% for CopraRNA versus 14%
for IntaRNA, and 0% for IntaRNA in *P. multocida*. It is structurally unavailable to us
because one bacterial genome is present per library by design.

**Module B should state openly that it works without the field's strongest feature**
rather than hoping no one notices. Whether a multi-strain *Wolbachia* comparative mode is
worth adding later is a genuine open question — it would be the single biggest available
improvement, and it fits the sRNA atlas project.

---

## 6. Validation plan

Since no benchmark set exists (§1.2):

1. **Simulation is primary.** Extend `simulate` with planted sites per class and per mode;
   `evaluate` reports recovery per class and **q-value calibration** — unplanted candidates
   should yield ≈ `fdr` of their calls. A miscalibrated q is a bug; a low recovery is not
   necessarily one.
2. **Median rank of the true target** on simulated data, reported as the headline metric,
   following bacterial-field practice.
3. **Sanity/regression set** from the 4 grade-A literature pairs (Bc-siR3.2→*AtMPK1/MPK2*;
   Fol-milR1→*SlyFRG4*; *M. incognita*→tomato; ghr-miR166/miR159→*V. dahliae*). These test
   that the plant mode recovers known cleavage-type pairs at all. **They are not a
   benchmark** — they were found by the tool class being tested — and must be labelled as
   regression fixtures in the code.
4. **Negative control:** WsnRNA-46 vs mouse *Il1rl1* 3'UTR is the one published clean
   non-target (hpo-miR-71 system); include as a fixture.
5. **No AUC/ROC figure in the paper.** If a reviewer asks for one, the answer is §1.2.

---

## 7. Open questions for the maintainer

1. **Does the claim change?** §2 proposes reframing module B from "target prediction" to
   "calibrated hypothesis generation". This affects the paper as much as the code.
2. **500 decoys vs 20** — acceptable cost on the shared server, or settle for 200?
3. **Null B donor organism** — needs a bacterium with no shared host. A free-living
   alphaproteobacterium? *Rickettsia* shares arthropod hosts and would be conservative.
4. **Is `ago_plant` worth building** given the group's focus, or is generality not worth
   the maintenance?
5. **Comparative mode** (§5, last paragraph) — the biggest available win, unavailable in
   the current one-genome design. Worth a schema change to allow multiple related
   symbiont genomes?

---

## 8. Data

The atlas project gains support from this review. *Wolbachia* has **one** purpose-built
WsnRNA dataset (GSE55210: one host, cell lines). A multi-host reanalysis across *Aedes*,
*Drosophila*, *Laodelphax*, *Tetranychus* and *Zyginidia*, with a calibrated null, would be
genuinely novel rather than confirmatory — and would be the first *Wolbachia* sRNA work
with a stated false discovery rate.

Accessions captured for reuse: GSE55210, GSE45320/45321/45323, GSE71598, GSE80421,
GSE55978/55941, GSE157226, PXD021844.

Note for scoping: no *M. tuberculosis*, *Salmonella* or *Listeria* sRNA → specific human
mRNA pair exists above prediction stage, in case those come up as extensions.
