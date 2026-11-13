# Strand C — Prediction tools for small-RNA → mRNA targeting: what exists, what it was validated on, and how well it actually works

Survey for **SmithHunter2 module B** (bacterial sRNA → host mRNA target prediction).
Compiled 2026-11-13. All claims paraphrased; reference numbers in square brackets map to §8.

---

## 0. Executive summary (read this if nothing else)

1. **No tool in any of the three families was built for the problem module B actually has** — a bacterial sRNA acting on a eukaryotic host transcriptome, with unknown mechanism, no Hfq, no AGO guarantee, and only one bacterial genome available.
2. **The bacterial-sRNA field's best tool (CopraRNA) wins almost entirely because of comparative genomics**, which SmithHunter2 structurally cannot use (one genome per library, by design). Strip conservation away and the remaining tools (IntaRNA, RNAup, RNAplex) plateau at a *median best-target rank in the tens-to-hundreds* on a 4,000-transcript search space [1,2].
3. **The miRNA field's best models are quantitative and well-calibrated but only because they were fitted to ~10⁵–10⁶ measurements of repression and binding affinity for *specific* miRNAs in *specific* cell types** [16,17,19]. Their parameters (3'UTR AU content, distance from stop codon, AGO occupancy, flanking-dinucleotide K_D) are AGO-biochemistry constants. Transplanting them to a bacterial sRNA in a host cell is unjustified; transplanting the *site-class scaffold* (8mer/7mer-m8/7mer-A1/6mer) is defensible, because that is just a k-mer taxonomy.
4. **Cross-kingdom tooling is essentially unvalidated.** The field's standard is "run psRNATarget or RNAhybrid, report thousands of hits, validate 1–2 by reporter assay." Multiple well-powered studies that went looking for the phenomenon found no functional evidence at all [31,32].
5. **Therefore: combine, don't choose.** Score seed-class and thermodynamic accessibility as *independent* lines of evidence, calibrate each against its own empirical null, and report the agreement explicitly. Do not multiply them into one opaque score.
6. **Build the FDR from dinucleotide-shuffled sRNA decoys plus a second, orthogonal null (real sRNAs from an unrelated bacterium).** Mononucleotide shuffling is known to be the wrong null for anything thermodynamic because it destroys stacking statistics [27,28].
7. **Do not trust RNAhybrid's or IntaRNA's built-in p-values on this problem.** RNAhybrid's EVD parameters are fitted on fly/worm/human 3'UTRs [25]; IntaRNA/CopraRNA fit a GEV per whole-genome run on a *bacterial* transcript population [4]. Both must be re-derived empirically from your own decoys.
8. **Precursor folding is close to worthless as evidence for a bacterial sRNA.** MFE z-scores do not separate mRNA from random sequence at all, and separate tRNA only weakly [29]. Implement it as an annotation field, never as a filter.

---

## 1. Family 1 — Bacterial sRNA → bacterial mRNA predictors

### 1.1 The tools

| Tool | Year | Interaction model | Seed constraint | Accessibility | Comparative genomics | Scale-out |
|---|---|---|---|---|---|---|
| **TargetRNA** [8] / **TargetRNA2** [9] | 2008 / 2014 | Hybridisation energy; TargetRNA2 adds sRNA & mRNA secondary structure, sRNA conservation, optional RNA-seq co-differential expression | implicit | partial (structure terms) | yes (conservation of sRNA) | ~1 min / transcriptome [1] |
| **IntaRNA** [3,4,5] | 2008 → 3.x | ΔG_hybrid + ΔG_unfold(query) + ΔG_unfold(target), exact DP | explicit, tunable (`--seedBP`, `--seedMaxE`, `--seedNoGU`, `--seedMinPu`) | full ED model, local windows | no (single sequence) | minutes–hours; 38 h for 93,280 *E. coli* pairs at defaults [6] |
| **RNAplex** [10] | 2008/2011 | simplified duplex energy + accessibility profile (RNAplfold) | no hard seed | yes, precomputed | no | fastest exact-ish tool; ~1 min / transcriptome [1]; 3 orders of magnitude faster than full methods [10] |
| **RNAup** | 2006/2008 | full opening-energy + duplex, most rigorous | no | yes (exact) | no | slow |
| **RNApredator** [10] | 2011 | RNAplex backend + accessibility + GO post-processing; 2,155 genomes preloaded | no | yes | no | fast web server |
| **sTarPicker** [11] | 2011 | two-step: stable duplex seeding → extension; ensemble ML classifier on seed + binding-region energy and accessibility | yes (two-step) | yes | no | moderate |
| **CopraRNA** [2,4] | 2013/2014 | IntaRNA per organism → GEV-fitted p-values → **Hartung combination for dependent p-values**, phylogenetically weighted → q-values | inherits IntaRNA | inherits IntaRNA | **yes — this is the whole point** | slow; "prohibitively slow" for many sRNAs [7] |
| **sRNARFTarget** [6] | 2021 | Random Forest on trinucleotide-frequency *differences* between sRNA and mRNA — **no thermodynamics at all** | none | none | none | ~100× faster than IntaRNA |
| **TargetRNA3** [7] | 2023 | Gradient boosting over 9 selected features from 111 candidates; top SHAP contributors = RNAplex hybridisation energy + number of sRNA:target homologs | seed lengths 8–10 had **little predictive power** | indirect | yes (homolog count is a feature) | seconds |
| **GraphRNA / kGraphRNA / sInterRF / sInterXGB** [12] | 2025 | GNN over sRNA–mRNA interaction graph; ensembles with CopraRNA | learned | learned (no external folding needed) | optional | fast |
| **sInterBase-trained / expression-augmented** [13] | 2025 | ML + co-expression modules from 2,143 *E. coli* / 177 *Salmonella* RNA-seq samples | — | — | — | — |

### 1.2 How small the training/benchmark sets really are

This is the single most important fact about this family.

| Benchmark | Size | Organisms |
|---|---|---|
| Richter & Backofen 2012 feature study [14] | **71** verified sRNA–target pairs | *E. coli*, *S. enterica* |
| Wright 2014 CopraRNA benchmark [4] | **101** verified pairs, 18 sRNAs | enterobacteria |
| Pain 2015 [1] | **102** trusted pairs, **22 sRNAs**, all *E. coli*; 82 with direct interaction evidence, 55 with compensatory mutations | *E. coli* only |
| Lai & Meyer 2015 [15] | 102 bacterial sRNA–mRNA + 54 fungal snoRNA–rRNA | mixed |
| IntaRNA 2.0 test set [6] | **119** pairs | mixed |
| sRNARFTarget training [6] | **745** positives (+745 decoys made by randomly swapping sRNAs) | mixed; from sRNATarBase 3.0 + MAPS/GRIL-seq/CLASH/RIL-seq |
| Raden 2020 parameter sweep [5] | **149** verified pairs, **15 sRNAs** | *E. coli* + *S. typhimurium* (4,319 / 4,552 targets) |
| TargetRNA3 [7] | **4,386** interactions from 77 sRNAs in 13 genomes, 4 phyla (largest to date; mostly from RIL-seq/CLASH/GRIL-seq), vs 325,162 labelled non-interactions | multi-phylum |

Two decades of tool development were tuned against **~100 *E. coli* pairs**. TargetRNA3's 4,386 is the first set big enough to train on, and it comes from high-throughput ligation methods whose interactions are *physical*, not necessarily *regulatory* — see §1.5.

### 1.3 Head-to-head results

**Pain et al. 2015** [1] — 22 *E. coli* sRNAs vs 4,317 mRNAs (window −200/+100 around start codon):

| Tool | Median rank of best true target | AUC (default UTRs / RNA-seq UTRs) | Base-pair recall within predicted hybrid |
|---|---|---|---|
| CopraRNA | **4** | 46.41 / — | 76.7% |
| IntaRNA | — | 27.02 / 28.83 | 76.7% |
| RNAup | — | 22.75 / 27.59 | 78.9% |
| RNAplex | — | 20.81 / 21.88 | 73.6% |
| TargetRNA2 | — | 18.52 / — | 55.9% |
| RNAhybrid | — | 4.14 / 4.38 | — |
| RNAduplex | — | 3.86 / 3.55 | — |
| RNAcofold | 1340 | 0.57 / 1.42 | — |
| Pairfold | — | 0.25 / 1.22 | — |

Runtime per sRNA against 4,317 mRNAs: ~1 min (TargetRNA2, RNAplex) to ~120 min (Pairfold) [1].

Three conclusions that matter for us:
- **CopraRNA wins because it adds conservation on top of accessibility** — and it needs homologs (min 3 organisms, ≥5 recommended, max 8) [4]. Without them it is just IntaRNA.
- **Fixed windows are fine.** Replacing the −200/+100 window with RNA-seq-derived real UTRs gave only 4.9–17.5% AUC gains, none significant [1]. This retires a whole class of worry.
- **RNAhybrid is near the bottom** (AUC 4.14 vs IntaRNA's 27) on bacterial sRNA–mRNA pairs, because it models neither accessibility nor the unfolding cost of a long structured sRNA [1].

**Umu & Gardner 2016** [16] — across all domains of life: the MFE tools that add accessibility (RNAup, IntaRNA, RNAplex) had the best true-positive rates at high PPV in every dataset, and could separate roughly half of native interactions from background. Alignment-based and internal-binding-energy methods had high TPR with low PPV.

**Lai & Meyer 2015** [15] — of 14 general RNA–RNA interaction methods, the best were **non-comparative, energy-based, accessibility-using tools predicting short interactions**. They explicitly warn that accuracy degrades with input length, "causing implications for de novo transcriptome-wide searches."

**sRNARFTarget 2021** [6] — recall of true targets in the top N%:

| Organism | cutoff | CopraRNA | sRNARFTarget | IntaRNA |
|---|---|---|---|---|
| *E. coli* | top 10% | 74% | 21% | 14% |
| *Synechocystis* | top 50% | 100% | 70% | 55% |
| *P. multocida* | top 20% | 18% | 10% | **0%** |

AUROC on the IntaRNA 2.0 119-pair set: sRNARFTarget 0.61 vs IntaRNA 0.59; *E. coli* subset 0.61 vs 0.62; *Salmonella* 0.58 vs 0.51 [6]. These are barely-above-chance numbers. Runtime: *P. multocida* 1,804 pairs in 31 s vs IntaRNA's 6,196 s; *E. coli* 93,280 pairs, 15 min vs >38 h (146×) [6].

**Raden et al. 2020** [5] — parameter sweep on IntaRNA, 149 pairs, 15 sRNAs. **Just changing parameters recovered ~30% more verified interactions at ~20% of baseline runtime.** Seed constraints alone: +26% recovery at half runtime. Best individual settings:

- `--seedBP 8` best for recovery; performance **collapses at ≥9**; 7 used as reference
- `--seedNoGU` — forbidding GU in seeds significantly improved recovery *and* cut runtime
- `--seedMaxE ≈ -5` kcal/mol with GU allowed, `≈ -2` without GU
- `--seedMinPu` best in 0.001–0.1
- interaction length: little gain above 50 nt; **~60 is the runtime/recovery sweet spot**
- interaction energy threshold: no effect above −10 kcal/mol, drops below
- `--seedMaxEhybrid`: no significant gain
- SHAPE data on targets (ShaKer): +4 verified pairs — marginal
- Energy model matters: Andronescu-2007 +5 pairs, Langdon-2018 +4, Turner-1999 +3 vs Turner-2004

**Published recommended set** [5]: 7-bp canonical seeds, no GU in seed, min unpaired probability 0.001 for both seed sites and both interaction sites, max interaction length 60, max interior/bulge loop 8.

**TargetRNA3 2023** [7]: at FPR 0.05 it had the highest sensitivity of six tools (CopraRNA, RNAup, IntaRNA, RNAplex, sRNARFTarget, TargetRNA3). Its SHAP analysis is the headline result for us: the two dominant features were **RNAplex hybridisation energy** and **number of sRNA:target homologs**. Seed length 8–10 contributed little. Again: energy + conservation, not seed.

### 1.4 The "binding ≠ regulation" problem

- Beisel et al. 2012 [17]: for Spot 42, restricting prediction to the **three unstructured regions** of the sRNA substantially improved recovery of known and novel targets. Predicted non-targets failed for identifiable reasons: no Hfq site, occluded target site, or overlapping Hfq/targeting sites.
- Richter & Backofen 2012 [14]: on 71 pairs, **only the sRNA-side interaction sites are sequence-conserved, not the target-side sites**; base-pairing complementarity is generally *not* conserved across more distant species. Accessibility, by contrast, is significantly higher at functional sites on **both** partners.
- Peer & Margalit 2011 [18]: experimentally determined sRNA target-binding regions are both more conserved and more accessible than random regions.
- Safari et al. 2026 [19]: integrating Hfq-CLASH with matched transcriptomics/proteomics in *E. coli*, an interpretable ML model predicted functional outcome at **AUC 0.78 (transcriptome) / 0.74 (proteome)**. Physical pairing alone was insufficient; target-side Hfq occupancy predicted functional regulation, AR2-domain occupancy on the sRNA predicted non-responsiveness.
- Reißer et al. 2021 [20] — the most uncomfortable finding: across sRNA–mRNA and snoRNA–rRNA datasets, **prediction performance was independent of the quality of the estimated accessibility profiles**. Estimated accessibility tracked sRNA structure reasonably but disagreed with structure annotations of longer transcripts. Their conclusion: accessibility's main practical effect is to *stop the thermodynamic model predicting implausibly long interactions*. It is acting partly as a length regulariser, not purely as biophysics.

### 1.5 What this family says to SmithHunter2

- Conservation is the strongest single signal and is **unavailable** to module B by design (one bacterium per library, host genome usually not a comparative set).
- Of the signals that remain, **accessibility-aware hybridisation energy** is the best validated, and **RNAplex-class energy is as informative as anything fancier** [7,10].
- Absolute performance without conservation is poor: AUC ~0.6, top-10% recall ~14–21%. **Any claim module B makes must be calibrated, not absolute.**

---

## 2. Family 2 — Animal miRNA → mRNA predictors

### 2.1 Tools and models

| Tool | Model | Training / benchmark data | Reported performance |
|---|---|---|---|
| **TargetScan 1–5** [21,22] | Seed classes (8mer > 7mer-m8 > 7mer-A1 > 6mer > offset-6mer) + conservation (branch length) + P_CT | Conserved-site counting across vertebrate UTR alignments with GC/length/PUF-matched control k-mers | Signal:background at branch-length cutoff 1.0: **8mer 2.6, 6mer 1.2**; ~534 ± 25 conserved sites and ~497 ± 49 conserved targets per miRNA family; **57.8 ± 3.0%** of human RefSeq genes are conserved targets at cutoff 0.6, authors estimate ~two-thirds overall [22] |
| **TargetScan context / context+** [23] | 5 context features: local AU content, co-targeting proximity, **3' supplementary pairing at miRNA nt 13–16**, ≥15 nt from stop codon, position away from centre of long UTRs | Transfection microarrays + endogenous miRNA data | Quantitatively predicts site efficacy without conservation; identifies non-conserved functional sites and siRNA off-targets [23] |
| **+ SPS / TA** [24] | Seed-pairing stability and target-site abundance as *miRNA-level* proficiency modifiers | HeLa reporters + arrays; lsy-6 and miR-23 as weak-proficiency cases | Weak SPS + high TA → low proficiency; improved prediction when modelled [24] |
| **TargetScan 7 (context++)** [25] | Site type + **14** additional features | Compendium of transfection datasets, pre-processed to remove confounding bias | Outperformed previous models; "as informative as the best high-throughput in vivo crosslinking approaches"; **showed that reported non-canonical CLIP sites bind but do not repress** [25] |
| **miRanda-mirSVR** [26] | Alignment + SVM regression on continuous down-regulation | ~Transfection expression datasets | Regression over classification; scores correlate with repression magnitude |
| **PITA** [27] | **ΔΔG = ΔG_duplex − ΔG_unfolding of the site + flank** | 190 *Drosophila* literature pairs (102 functional, 88 non-functional) + >60 luciferase interactions in S2 cells, 4–8 replicates | **AUC 0.79** (3 nt up / 15 nt down flank) vs Stark 0.74, PicTar 0.71, miRanda 0.64, STarMir ~0.49. ΔΔG vs measured repression **r = 0.70** (P < 4×10⁻⁴); **ΔG_duplex alone only r = 0.36** |
| **RNA22** | Pattern/motif-based, no seed requirement, target-first | — | Widely used in cross-kingdom work precisely because it has no seed rule |
| **DIANA-microT** | Seed + conservation + weighted site score | — | — |
| **RNAhybrid** [28] | Pure MFE duplex, no accessibility; EVD p-value | 3'UTR parameter sets for fly/worm/human | Bottom of the bacterial benchmark (AUC 4.14) [1] |
| **McGeary et al. 2019 — biochemical + CNN** [29] | AGO-RBNS relative K_D for all k-mers ≤12 nt; biochemical occupancy model; CNN generalises K_D to any miRNA sequence | 6 miRNAs by AGO-RBNS (miR-1, let-7a, miR-7, miR-124, miR-155, lsy-6); >1.5M relative K_D values; 68,112 expression estimates over 4,257 transcripts in 16 transfections; tested on 12 held-out miRNAs | **Biochemical model r² = 0.34** for mRNA changes vs **r² ≤ 0.14 for 30 prior algorithms**; TargetScan7 0.14 → 0.25 on same data, 0.28 retrained; biochemical+ (accessibility, supplementary pairing, conservation) **0.37**; CNN-predicted K_D vs repression **r² = 0.76** (RNAduplex 0.21, MIRZA 0.56); on held-out miRNAs biochemical+ with CNN inputs **0.20** vs TargetScan **0.13** |
| **scanMiR** [30] | McGeary k-mer affinity, as lightweight linear models + aggregated transcript repression; flexible 3'-supplementary alignment | Bioconductor | Practical reimplementation of [29] |
| **TargetNet** [31], **miRBind2**, **MiRformer**, **DeepExoMir** | Deep nets on miRNA–CTS pairs with relaxed seed criteria | CLASH/CLIP-derived pairs | TargetNet improves on prior classifiers; but see [32] |
| **miRBench** [32] | Bias-controlled negative generation | New curated datasets | Exposed a **miRNA frequency class bias** in standard negatives; a *simple CNN* retrained on debiased data reached **average precision 0.81–0.86**, beating state-of-the-art — i.e. much reported deep-learning gain was dataset artefact |

### 2.2 The single most important number: the flanking-dinucleotide effect

McGeary et al. [29] measured a **~100-fold range in binding affinity for the miR-1 8mer purely from the two flanking dinucleotides** — frequently exceeding the difference between *site types*. For each miRNA tested, the best-flanked 6-nt canonical site bound more strongly than the 8mer in its worst flanking context. A linear model of flanking effects alone reached r² = 0.63.

This is the empirical core of what "3'UTR AU content" was always a proxy for, and it is the reason a pure site-class call (8mer vs 7mer) is a weak predictor of anything quantitative.

Site occupancy figures from the same work: at 7.3 pM AGO2–miR-1, the **8mer occupied only 3.8–17%** of the complex; all marginal sites together 6.1–9.8%; **reads with no identifiable site occupied 32–53%**.

### 2.3 What fraction of predictions are functional?

- Friedman et al. [22]: at the 8mer level with conservation filtering, S/B = 2.6 implies **P_CT ≈ (2.6−1)/2.6 ≈ 0.62**, i.e. ~38% of even *conserved* 8mer sites are chance. For 6mers, S/B = 1.2 → **P_CT ≈ 0.17**: five in six conserved 6mers are noise.
- Pinzón et al. 2016 [33]: miR-223-guided repression is **often smaller than inter-individual variation in gene expression among wild-type mice** — i.e. most predicted targets are functionally insensitive to the miRNA. Highly conserved sites cluster in haplo-insufficient (dose-sensitive) genes. Many conserved sites are conserved for **miRNA-independent reasons** and are only incidentally complementary. Their conclusion is blunt: the role of miRNAs has been overestimated because false-positive rates are routinely ignored.
- Helwak et al. CLASH [34]: 18,514 miRNA–mRNA interactions, 399 miRNAs, 6,959 genes, but chimeras were only ~2% of AGO1 reads. **Only ~37% of seed interactions had uninterrupted Watson–Crick pairing**; ~60% of seed interactions were non-canonical; a seed-less class (nt 13–21 pairing) made up ~16%. Location: ~50% CDS, ~42% 3'UTR, ~4–5% 5'UTR after junction filtering. The class combining seed with **nt 13–16 supplementary pairing** was both the most conserved and the most strongly repressive. Non-seed effects were statistically significant but modest.
- Agarwal et al. 2015 [25] directly contradicted the CLASH expansion: the non-canonical sites bind but **do not mediate repression**; the vast majority of functional sites are canonical.
- Target abundance / AGO competition: Arvey et al. [35] showed across 178 transfections in two cell lines that per-gene down-regulation varies with total available target-transcript concentration; they argue targeting must be treated quantitatively, not as a binary label. Denzler et al. and Khan et al. [36] showed transfection itself globally perturbs endogenous miRNA regulation — a confounder in most of the benchmark data.

### 2.4 What this family says to SmithHunter2

- The **site-class taxonomy (8mer / 7mer-m8 / 7mer-A1 / 6mer / offset-6mer) is a sound organising scheme** — it is a k-mer hierarchy with a monotone affinity ordering that was reproduced biochemically [29].
- The **quantitative context model is not portable**: AU content, UTR position and stop-codon distance are AGO-loading and UTR-biology constants, measured in human/mouse/fly cells, with an A1 anchor that is a property of the **AGO MID domain**, not of RNA. A bacterial sRNA has no reason to obey them.
- **An A1 requirement is the most suspect portable feature.** 7mer-A1 and 8mer both require an adenosine opposite guide position 1, enforced by AGO's nucleotide-specificity pocket. If the host AGO is not loading the bacterial sRNA, the A1 classes lose their mechanistic basis entirely.
- The **3'-supplementary pairing at guide 13–16** is the one context feature with a generic base-pairing rationale, and it is cheap to compute. Worth carrying as an annotation.

---

## 3. Family 3 — Cross-kingdom and plant tools

### 3.1 The tools

| Tool | Model | Validation status |
|---|---|---|
| **psRNATarget** (2011, 2017) [37,38] | Penalty-score complementarity (plant-style near-perfect pairing) + **target-site accessibility (UPE)**; customisable for canonical and non-canonical | Trained/tuned on plant miRNA–mRNA pairs; 2017 release claims higher recall without inflating total output. **No cross-kingdom validation.** Used in cross-kingdom work at expectation cutoffs as loose as **E ≤ 2** [39], which is permissive |
| **TAPIR** [40] | Two modes: fast (FASTA-like) and precise (RNAhybrid-based); supports target mimics with large loops | Plant miRNAs; no cross-kingdom benchmark |
| **TargetFinder** | Penalty score | In Arabidopsis: **89% precision, 97% recall**; in non-Arabidopsis species only **46%** of true interactions detected; score optimisation lifted recall to 70% at 65% precision [41] |
| **P-TarPmiR** [42] | TarPmiR retrained on plant-only TarDB data | Key finding: **excluding animal training data entirely produced the most accurate plant predictor** — animal-derived features actively degraded plant prediction |
| **RNAhybrid / miRanda / RNA22** | Generic duplex | The de facto cross-kingdom tools, used with no organism-appropriate calibration |
| **Zanini et al. pipeline** [43] | A practical review-style procedure for finding/characterising/validating sRNA effectors in plant–microbe interaction | A guideline paper, not a benchmarked tool |

### 3.2 What validation actually exists

Honest answer: **almost none, and the negative results are as strong as the positive ones.**

Positive / partially validated:
- Weiberg et al. 2013 [44]: *Botrytis cinerea* sRNAs hijack host AGO1 — the founding case, with genetic validation.
- Cai et al. 2018 [45]: plants send sRNAs in extracellular vesicles to silence fungal virulence genes.
- Wong-Bajracharya et al. 2022 [46]: *Pisolithus microcarpus* miRNA in cross-kingdom silencing during symbiosis.
- Silvestri et al. 2024 [47]: a fungal sRNA silences a host plant transcription factor, with degradome support.
- Johnson et al. 2019 [48] — methodologically the most interesting: *Cuscuta* trans-species sRNAs target host mRNAs almost exclusively in **highly conserved protein-coding regions**, and vary in a **three-nucleotide period** to compensate for synonymous-site variation across hosts. If module B is looking at a real cross-kingdom effector, this is the signature to look for: CDS-biased, codon-wobble-tolerant targeting — *not* 3'UTR targeting.

Negative / failed-to-replicate:
- Kettles et al. 2019 [31]: *Zymoseptoria tritici* — 389 in-planta-induced sRNA loci predicted to target wheat defence mRNAs; **molecular validation failed for every tested target**; DCL and AGO deletions were dispensable for full virulence; the fungus cannot even take up dsRNA.
- Qin et al. 2022 [32]: *Botrytis*–tomato — 7,042 fungal sRNAs predicted to target 3,185 tomato mRNAs; 163 of the predicted targets were indeed down-regulated during infection; but deleting the transposon locus producing ~10% of the sRNAs, and deleting both Dicer-like genes (>99% reduction in transposon-derived sRNA), **did not reduce virulence on any plant tested**. "No evidence for any detectable role."

The gap between "7,042 sRNAs predicted to hit 3,185 transcripts" and "no detectable phenotype" is the exact failure mode module B must be engineered against.

### 3.3 Tools in the Wolbachia / parasite sRNA literature

Mayoral et al. 2014 (PNAS) [49] is the reference case and is instructive about the field's standards:
- sRNAs **WsnRNA-46 and WsnRNA-59** (~30 nt, from ~65–70 nt stem-loop precursors), detected by Northern blot in *A. aegypti*, *D. melanogaster*, *D. simulans*.
- Target prediction: **BLAST, RNAhybrid and RNA22**, against the *A. aegypti* and *w*MelPop-CLA genomes.
- **No cutoffs are reported** for any of the three. The paper describes a "stringent target prediction" without defining it. The only number given is MFE −30.3 kcal/mol for the host *Dhc* interaction.
- Host target *Dhc* was selected partly for **complete complementarity to guide positions 2–8** — i.e. a miRNA seed rule applied to a bacterial sRNA, with no justification offered.
- Validation: mimic/inhibitor RT-qPCR plus GFP reporters with 177 bp (*murD*) and 181 nt (*Dhc*) fragments; seed-mutated reporters lost the effect. The *Dhc* effect was an **increase**, not repression.
- **No FDR, no decoy control, no quantitative threshold.**

Barnard/Hussain-lineage work on Wolbachia-modified mosquito miRNA profiles is a different claim (host miRNA modulation, not bacterial sRNA → host mRNA) and does not supply target-prediction validation either.

**Plain statement for the manuscript:** there is at present no validated computational method for predicting bacterial-sRNA → animal-host-mRNA targets. SmithHunter2 module B would be building the first one that is explicitly calibrated, and should say so rather than implying inherited validity from IntaRNA or TargetScan.

---

## 4. Design question 1 — Which scoring model is least wrong across kingdoms?

**Answer: both, scored separately, reported as independent evidence, never collapsed into one number.**

Reasoning:

1. **The seed-centric model's mechanistic warrant is AGO.** The 5'-seed dominance, the A1 pocket, the 8mer > 7mer-m8 > 7mer-A1 > 6mer ordering — these are properties of Argonaute's binding channel, re-measured directly by AGO-RBNS [29]. If the bacterial sRNA is loaded into a host AGO, the model transfers. If it acts by direct base-pairing (the bacterial mechanism), it does not, and the bacterial field's own evidence says seed length 8–10 carried **little predictive power** in the largest ML feature study [7].
2. **The thermodynamic/accessibility model's warrant is physics**, and it is the only one that is kingdom-agnostic. It is also the best-validated signal in the bacterial benchmarks [1,15,16]. But Reißer et al. [20] showed its contribution is partly as a length regulariser, so it should not be oversold either.
3. **They fail in different directions**, which is the whole argument for combining them. Seed scoring is blind to structure and to interaction length; energy scoring systematically favours long, GC-rich, highly structured transcripts regardless of biology. A hit that both nominate is informative in a way that neither alone is.
4. **What people facing this choice actually do:** Mayoral et al. used *both* (RNAhybrid for energy, a seed-complementarity requirement for site selection) [49]. The plant cross-kingdom community uses psRNATarget, which is itself a complementarity-penalty score **plus** an accessibility (UPE) term [37,38]. PITA's own data say why: ΔΔG correlated with repression at r = 0.70 while duplex energy alone gave r = 0.36 [27] — accessibility is not optional once you are doing thermodynamics at all.
5. **Concretely, the opt-in "bacterial mode" is the right architecture.** Default mode = seed-constrained (the host might be using AGO); bacterial mode = no seed rule, start-codon windows (the sRNA might be acting as a bacterial sRNA). Run both, and report the union with the mode flag attached, plus a flag for hits recovered by both. Do not try to pick the right biology computationally — you cannot, and the output should say so.

**One concrete refinement:** given Johnson et al. [48], add a CDS-targeting evidence channel. If the real mechanism is cross-kingdom effector targeting, the signal is in conserved coding regions, and a 3'UTR-only or 5'UTR-only search will miss it by construction. CLASH found ~50% of miRNA sites in CDS anyway [34].

---

## 5. Design question 2 — How to build an empirical FDR

### 5.1 The candidate nulls, and what each preserves

| Null | Preserves | Destroys | Verdict |
|---|---|---|---|
| **Mononucleotide shuffle of the sRNA** | base composition | dinucleotide/stacking statistics, folding propensity | **Do not use for anything thermodynamic.** Stacking energies are dinucleotide quantities; mononucleotide shuffles systematically *destabilise*, so real sequences look spuriously significant [28,29-folding]. The folding literature is explicit that mononucleotide shuffling should not be used for z-scores [29-folding] |
| **Dinucleotide shuffle of the sRNA** (Altschul–Erikson) | base and dinucleotide composition, hence approximately the stacking free-energy budget and GC content | the actual seed k-mer, the real secondary structure, any evolved accessibility pattern | **The right primary null.** It is the standard control for MFE-based inference and the one Gesell & Washietl argue for [28] |
| **Shuffling the targets** | nothing useful | target 5'UTR/CDS composition structure, codon bias, start-codon context, local GC landscape | **Biased.** It changes the search space, not the query. Any composition bias in real UTRs (AU-rich 5'UTRs in AT-rich hosts; codon-driven periodicity in CDS) is destroyed, so the null is systematically *easier* than reality and the FDR is anti-conservative. Use only as a secondary diagnostic |
| **Real sRNAs from an unrelated organism as decoys** | all realistic sequence statistics: real length distribution, real folding, real k-mer vocabulary, real 5'-end nucleotide bias | nothing — which is also the risk | **The best orthogonal null, with one caveat.** It controls for everything a shuffle cannot. The caveat: if the unrelated organism shares host range or GC content, decoys may be genuine targeters, inflating the null and making the FDR conservative. Pick a phylogenetically and ecologically distant donor and say which |
| **Permuting sRNA–target assignments** | both marginal distributions exactly | the pairing itself | **Correct null only when you already have a set of candidate pairs** and are asking whether a *specific assignment* is special. This is what sRNARFTarget used for negatives [6], and what miRBench showed can go wrong: random swapping introduced a **miRNA frequency class bias** that models learned instead of biology [32]. Safe for evaluation; unsafe as ML training negatives |

### 5.2 The known critique

- **Rivas & Eddy (2000)** established the general result that MFE-based significance against shuffled sequence is not sufficient to detect ncRNA at genome scale; the z-score literature has been living with this ever since [29-folding cites it as unresolved].
- **Gesell & Washietl (2008)** [28] built dinucleotide-controlled null models precisely because naive shuffling gives the wrong background for RNA structure inference.
- **Empirically** [29-folding]: against dinucleotide-preserving shuffles, mean z-scores were SECIS −4.74, U2 snRNA −4.22, srpRNA −3.56, hammerhead III −3.19, U1 −1.75, **tRNA only −1.59 (p ≈ 0.12)**, and **mRNA ≈ −0.18 (p ≈ 0.48)**, 5'UTR **+0.18**, 3'UTR −0.11, CDS −0.13. Under a *weaker* randomisation the same ordering appeared with uniformly smaller magnitudes — i.e. the choice of shuffle directly scales your apparent significance.
- **miRBench** [32] is the modern restatement for target prediction: in-silico negatives are "challenging and prone to introducing biases," and models learn dataset artefacts rather than biology.

### 5.3 Recommended construction for module B

Use **two nulls and require agreement**, because a single shuffle-based FDR is exactly the thing the literature warns about.

- **Null A (primary, per-sRNA):** 200–500 Altschul–Erikson dinucleotide-preserving shuffles of the sRNA, each scored against the **identical, unshuffled** target set with the identical pipeline. Empirical q-value = (number of decoys achieving a score ≥ observed, across all decoys and targets) / (number of decoys × targets), Storey–Tibshirani corrected.
- **Null B (orthogonal, per-run):** 200+ real sRNAs of matched length drawn from a phylogenetically distant bacterium with no shared host, scored identically.
- **Report both q-values.** Discrepancy between them is itself the diagnostic: if q_A ≪ q_B, the hit is explained by realistic sequence statistics the shuffle destroyed, and should be demoted.
- **Shuffle preserving guide-position constraints:** if you shuffle positions 2–8 along with everything else, the decoy seed k-mers will have the dinucleotide-conditioned composition of the real seed, which is what you want. Do **not** hold the seed fixed and shuffle the rest — that collapses the null for the exact feature you are testing.
- **Never shuffle the targets** as the primary null.
- **Stratify the FDR by site class and by transcript length decile.** A single pooled FDR across 8mers and 6mers conceals that the 6mer stratum is almost pure noise (§6).
- **Sanity check with permuted assignment** across your own candidate table, as an evaluation-only control, with the miRBench frequency-bias caveat stated.

---

## 6. Design question 3 — Realistic precision at transcriptome scale

### 6.1 The arithmetic

For a uniform-composition search space of total length L nucleotides, the expected count of exact matches to a single fixed k-mer is ≈ L / 4^k.

**Bacterial-mode window** (4,000 transcripts × 300 nt for −200/+100 around start codon ≈ **1.2 Mnt**), single seed drawn from guide positions 2–8:

| Site class | k | Expected chance hits in 1.2 Mnt |
|---|---|---|
| 8mer | 8 | **≈ 18** |
| 7mer-m8 / 7mer-A1 | 7 | **≈ 73** each |
| 6mer / offset-6mer | 6 | **≈ 293** each |

**Host-transcriptome mode** (20,000 transcripts × ~1,200 nt of 3'UTR ≈ **24 Mnt**):

| Site class | Expected chance hits |
|---|---|
| 8mer | **≈ 366** |
| 7mer | **≈ 1,465** each |
| 6mer | **≈ 5,860** each |

**If you scan every 7-nt window of the sRNA rather than pinning positions 2–8**, multiply by (sRNA length − 6). For a 100-nt bacterial sRNA that is ~94 seeds: **~6,900 chance 7mer hits in the bacterial window, ~138,000 in a host 3'UTR space.** This is the single biggest decision in the design. Seed pinning is not a nicety; it is a ~100× change in the base rate.

Corrections that make it worse: AT-rich hosts and AU-rich 3'UTRs raise the hit rate for AU-containing seeds well above the uniform estimate; and A1 "anchoring" is a single-base constraint that only costs a factor of ~4 while being mechanistically unjustified outside AGO.

### 6.2 How much do energy and accessibility move it?

Not nearly as much as one would like:

- **Conservation is the biggest lever, and even it leaves ~38% noise.** Friedman's 8mer S/B of 2.6 with branch-length ≥1.0 [22] means that after *both* seed matching *and* deep conservation filtering, over a third of surviving sites are chance.
- **Accessibility-aware energy lifts AUC from ~0.5 to ~0.6** in the bacterial setting [6] and from AUC 0.64 (miRanda) to 0.79 (PITA) in the *Drosophila* setting [27] — real, but far from a precision filter.
- **Top-10% recall of true targets is 14% (IntaRNA) to 74% (CopraRNA)** [6]; without conservation you are at the low end.
- **The ceiling on *functional* prediction is lower than the ceiling on *binding* prediction.** The best quantitative animal model explains r² ≈ 0.34–0.37 of repression [29]; the best bacterial functional model reaches AUC 0.74–0.78 and only by adding Hfq-CLASH occupancy data [19].
- **At cross-kingdom scale with permissive settings the output is uninterpretable by construction**: psRNATarget at E ≤ 2 yielded 1,677 predicted host transcripts for one fungal milRNA set [39]; Qin et al. got 3,185 predicted tomato targets and no phenotype [32].

### 6.3 The honest expectation to state in the paper

With a seed-pinned 7mer scan over a few thousand transcripts, **expect on the order of 10²–10³ chance hits per sRNA before any energy filtering**. Accessibility-aware energy plus an empirical q-value can plausibly reduce the reported set to tens of candidates at q < 0.1, but the *precision* of that set is unknown and, on all available evidence from better-resourced fields, is likely well under 50%. Module B should be positioned as a **ranked hypothesis generator with calibrated relative evidence**, and the documentation should say that an absolute precision claim is not supportable.

---

## 7. Design question 4 and 5 — Practical notes

### 7.1 Running IntaRNA well

Parameter choices that materially matter (from [5], the only systematic sweep):

| Setting | Recommended | Why |
|---|---|---|
| `--seedBP` | **7** (8 scored best on recovery but 7 is the safer reference; **≥9 collapses**) | Longer seeds over-constrain; the recovery cliff at 9 is sharp |
| `--seedNoGU` | **on** | Improved recovery *and* reduced runtime — a rare free win |
| `--seedMaxE` | **−5** kcal/mol with GU, **−2** without | Tuned optimum |
| `--seedMinPu` | **0.001** (published rec.) — sweep showed 0.001–0.1 all good | Seed-site accessibility |
| `--intLenMax` | **60** | No recovery gain above ~50; 60 is the runtime sweet spot |
| `--intLoopMax` | **8** | No effect alone; matters in combination |
| interaction accessibility | min unpaired prob **0.001–0.01** | Sweep peaked at 0.01; published rec. uses 0.001 |
| `--outMaxE` | no threshold above **−10** kcal/mol | Nothing above −10 contributes |
| energy model | consider **Andronescu-2007** as a secondary run | +5 verified pairs over Turner-2004 |

**The headline practical fact: default IntaRNA leaves ~30% of recoverable verified interactions on the table while running ~5× slower than the tuned configuration** [5]. Tuning is not optimisation theatre; it is the single largest performance change available to module B without adding data.

Further notes:
- **`--tRegion` vs whole-transcript accessibility.** The fixed −200/+100 window is empirically justified: RNA-seq-derived real UTRs gave only 4.9–17.5% non-significant AUC gains [1]. Use the window as the *search* region, but compute accessibility with a **local sliding window (`--tAccW`/`--tAccL`)**, not over the whole transcript — global folding of a long mRNA is both slow and, per Reißer et al. [20], disagrees with real structure annotations of long transcripts. Window ~150, max base-pair span ~100 is the conventional RNAplfold-style setting.
- **Query side:** for a long structured bacterial sRNA, compute `--qAcc` over the **full sRNA**, not a window — the sRNA's own fold is real and is what Beisel et al. exploited by restricting to unstructured regions [17]. Consider adding Beisel's trick explicitly: restrict `--qRegion` to predicted single-stranded regions of the sRNA, which substantially improved recovery for Spot 42.
- **Mode:** exact (`--mode H`) is affordable at a few thousand transcripts with the tuned parameters; use `--mode S`/seed-only heuristics only if profiling says otherwise. `--threads` scales near-linearly.
- **Pitfalls.** (i) Energies are not comparable across transcripts of different lengths — longer targets win spuriously; rank within length strata or regress out length before thresholding. (ii) GC-rich transcripts win spuriously for the same reason — this is exactly what the dinucleotide-preserving null is for. (iii) Overlapping suboptimals (`--outNumber`, `--outOverlap`) will inflate hit counts if not deduplicated per transcript; take one best interaction per transcript for FDR purposes.

### 7.2 Running RNAhybrid well

- **The EVD is the whole issue.** RNAhybrid normalises energy as `en = e / log(m·n)` for target length m and query length n, and models duplex energies as an extreme-value distribution with position parameter ξ and shape θ: `P[X ≤ en] = 1 − exp(−exp(−(−en − ξ)/θ))` [25]. If `-d ξ,θ` is not given, the parameters are **estimated from the query's maximal duplex energy assuming a linear dependence**, using one of three built-in sets: `3utr_fly`, `3utr_worm`, `3utr_human` [25].
- **None of those three sets is appropriate for a bacterial sRNA or for a non-model host.** Using `-s 3utr_human` on a Zyginidia or mosquito transcriptome with a 30-nt bacterial sRNA gives p-values that are decorative.
- **Therefore: always run `RNAcalibrate` on your actual target population and pass the resulting ξ,θ via `-d`.** This is the documented purpose of the companion tool and it is the only way RNAhybrid's p-values mean anything here.
- **Even then, treat the p-value as a within-run ranking statistic, not an error rate.** The EVD is fitted to a *sequence-composition* null, not to the biology; your empirical decoy q-value (§5) is the number to report.
- Other options: `-f from,to` **forces a helix between given 1-based query positions** — this is how you pin the seed in RNAhybrid (e.g. `-f 2,8`), and per §6.1 it is worth ~100× in base rate. `-e` energy cutoff interacts with `-b` (hit count) — `-b` stops reporting once reached, so a tight `-b` with a loose `-e` silently truncates. `-m` defaults to max target length 2000 and `-n` to max query length 30; **long bacterial sRNAs exceed `-n` 30 and will be rejected or truncated** — a classic silent failure. `-u`/`-v` cap unpaired nucleotides in internal loops and bulges.
- **RNAhybrid has no accessibility model at all**, which is why it sat near the bottom of the bacterial benchmark (AUC 4.14 vs IntaRNA 27.02) [1]. Use it as a fast, independent, structure-naive second opinion — which is genuinely useful for the "both models agree" channel — not as a primary scorer.

### 7.3 Are the reported p-values trustworthy?

**No, in both cases, and for different reasons.**

- RNAhybrid: parameters fitted on the wrong organisms' 3'UTRs by default [25].
- IntaRNA/CopraRNA: p-values come from fitting a **generalized extreme value distribution to the IntaRNA energy scores within each whole-genome prediction** [4]. That is a sensible within-run calibration for a bacterial transcriptome, but the fitted population is your target set — so the p-value answers "is this energy extreme for this transcriptome?", not "is this a target?". It inherits every composition bias in the target set.
- CopraRNA's q-values come from combining per-organism p-values by **Hartung's method for dependent p-values**, phylogenetically weighted with a root function to damp outliers, with homolog clusters from DomClust [4]. That machinery is unavailable without homologs.
- **Recalibrate empirically. Always.** Then report the tool p-value and the empirical q-value side by side.

### 7.4 Precursor / hairpin folding as evidence

The numbers (dinucleotide-preserving shuffles, mean z-scores) [29-folding]:

| Class | mean z | p |
|---|---|---|
| SECIS | −4.74 | 0.0 |
| U2 snRNA | −4.22 | 0.002 |
| srpRNA | −3.56 | 0.046 |
| Hammerhead III | −3.19 | 0.008 |
| U1 snRNA | −1.75 | 0.10 |
| **tRNA** | **−1.59** | **0.12** |
| **mRNA (full)** | **−0.18** | **0.48** |
| mRNA 5'UTR | **+0.18** | 0.46 |
| mRNA 3'UTR / CDS | −0.11 / −0.13 | 0.53 / 0.51 |

Read that table carefully: **tRNA — a molecule with a famous, universal, crystallographically determined fold — does not reach significance** on MFE z-score. mRNA is indistinguishable from random. Only the most extreme structural families separate cleanly.

Add the classical critique (Rivas & Eddy): structure alone is generally not statistically significant for ncRNA detection at genome scale, and the Carleton authors note their own work does not resolve that concern [29-folding]. The randfold tool exists precisely because raw MFE is not interpretable without a shuffled-background p-value — and even the randfold p-value is only as good as the shuffle (§5).

**Is folding worth implementing for SmithHunter2 module B?**

**No, not as evidence, and not as a filter. Yes, as a cheap annotation field.** Specifically:

- A bacterial trans-acting sRNA need not have a hairpin precursor at all. Many functional bacterial sRNAs are processed from 3'UTRs or are tRNA fragments — the Koeppen/Xie OMV sRNAs that started this whole line of inquiry are **tRNA and ncRNA fragments**, not hairpin products.
- Mayoral's WsnRNAs *did* come from ~65–70 nt stem-loops [49], so hairpin structure is sometimes present — but that is an observation to record, not a criterion to filter on.
- Filtering candidates on MFE or MFEI would import the known high false-positive rate of hairpin-based miRNA prediction into a setting where the premise (a Dicer-processed hairpin precursor) may be false by construction.

**Concrete recommendation:** compute RNAfold MFE, MFEI (MFE per nt normalised by GC content), and a randfold-style z-score **against dinucleotide-preserving shuffles** for every sRNA locus, store all three, and surface them in the output table with an explicit note that z > −2 does not argue against the sRNA. Never gate module B's target predictions on them. If later data shows hairpin-derived and non-hairpin-derived sRNAs behave differently, the field is already in the database and the analysis is a one-liner.

---

## 8. References

1. Pain A. *et al.* (2015) An assessment of bacterial small RNA target prediction programs. *RNA Biology*. https://consensus.app/papers/details/w:39fa655d-2c45-486d-bf04-e1cac8f36b8e/ (PMC4615726)
2. Wright P.R. *et al.* (2013) Comparative genomics boosts target prediction for bacterial small RNAs. *PNAS*. https://consensus.app/papers/details/w:dd131c39-814c-49d6-a446-64b3ed3e3ae6/
3. Busch A., Richter A.S., Backofen R. (2008) IntaRNA: efficient prediction of bacterial sRNA targets incorporating target site accessibility and seed regions. *Bioinformatics*. https://consensus.app/papers/details/w:... (see Freiburg RNA tools)
4. Wright P.R. *et al.* (2014) CopraRNA and IntaRNA: predicting small RNA targets, networks and interaction domains. *Nucleic Acids Research*. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4086077/
5. Raden M. *et al.* (2020) The impact of various seed, accessibility and interaction constraints on sRNA target prediction — a systematic assessment. *BMC Bioinformatics*. https://pmc.ncbi.nlm.nih.gov/articles/PMC6956497
6. Naskulwar K., Peña-Castillo L. (2021) sRNARFTarget: a fast machine-learning-based approach for transcriptome-wide sRNA target prediction. *RNA Biology*. https://www.biorxiv.org/content/10.1101/2021.03.05.433963v1.full.pdf
7. Tjaden B. (2023) TargetRNA3: predicting prokaryotic RNA regulatory targets with machine learning. *Genome Biology*. https://pmc.ncbi.nlm.nih.gov/articles/PMC10691042/
8. Tjaden B. (2008) TargetRNA: a tool for predicting targets of small RNA action in bacteria. *NAR*. https://consensus.app/papers/details/w:... (NAR 2008)
9. Kery M.B. *et al.* (2014) TargetRNA2: identifying targets of small regulatory RNAs in bacteria. *NAR*. https://consensus.app/papers/details/w:aa51f83c-414a-4069-b5f4-9c453ef4e5f8/
10. Eggenhofer F., Tafer H., Stadler P.F., Hofacker I.L. (2011) RNApredator: fast accessibility-based prediction of sRNA targets. *NAR*. https://consensus.app/papers/details/w:9f2732bd-7df9-49ef-9523-c3dbb2be73d4/
11. Ying X. *et al.* (2011) sTarPicker: a method for efficient prediction of bacterial sRNA targets based on a two-step model for hybridization. *PLoS ONE*. https://consensus.app/papers/details/w:fedd6d0a-1357-487e-a398-d2be8b03dff0/
12. Attias Cohen S. *et al.* (2025) GNNs and ensemble models enhance the prediction of new sRNA-mRNA interactions in unseen conditions. *BMC Bioinformatics*. https://consensus.app/papers/details/w:9fb2ec84-98e3-4ef9-ac92-8db18078b4c5/
13. Derinkok Y. *et al.* (2025) Improving prediction of bacterial sRNA regulatory targets with expression data. *NAR Genomics and Bioinformatics*. https://consensus.app/papers/details/w:47a44796-3fc9-482c-94be-cd88d850a51b/
14. Richter A.S., Backofen R. (2012) Accessibility and conservation: general features of bacterial small RNA–mRNA interactions? *RNA Biology*. https://consensus.app/papers/details/w:ba65c16a-7fa2-49cf-a4da-a90dcfbc1e66/
15. Lai D., Meyer I.M. (2015) A comprehensive comparison of general RNA–RNA interaction prediction methods. *NAR*. https://consensus.app/papers/details/w:24936be1-a3b8-4756-90b0-3cf6804677dd/
16. Umu S.U., Gardner P.P. (2016) A comprehensive benchmark of RNA–RNA interaction prediction tools for all domains of life. *Bioinformatics*. https://consensus.app/papers/details/w:bebdecd3-91cc-481f-b84e-14ab22062838/
17. Beisel C.L., Updegrove T.B., Janson B.J., Storz G. (2012) Multiple factors dictate target selection by Hfq-binding small RNAs. *EMBO J*. https://consensus.app/papers/details/w:650c0ed8-6a53-4dc2-9ba1-6a39226980dd/
18. Peer A., Margalit H. (2011) Accessibility and evolutionary conservation mark bacterial small-RNA target-binding regions. *J Bacteriol*. https://consensus.app/papers/details/w:ad47c655-6d39-4edb-9b11-93d7ae9c60ef/
19. Safari F. *et al.* (2026) Molecular determinants of functional bacterial sRNA–mRNA interactions revealed by integrating RNA interactomes and interpretable machine learning. *bioRxiv*. https://consensus.app/papers/details/w:c0dfd171-b806-44d1-8941-5d7ab0a05542/
20. Reißer S. *et al.* (2021) Investigating the concept of accessibility for predicting novel RNA-RNA interactions. *bioRxiv*. https://consensus.app/papers/details/w:201ca7ee-42b8-417a-93e0-29a18227fad0/
21. Lewis B.P. *et al.* (2003) Prediction of mammalian microRNA targets. *Cell*. https://consensus.app/papers/details/w:... (Cell 2003)
22. Friedman R.C., Farh K.K.-H., Burge C.B., Bartel D.P. (2009) Most mammalian mRNAs are conserved targets of microRNAs. *Genome Research*. https://genome.cshlp.org/content/19/1/92
23. Grimson A. *et al.* (2007) MicroRNA targeting specificity in mammals: determinants beyond seed pairing. *Molecular Cell*. https://consensus.app/papers/details/w:69ce10d5-403f-4e8b-a881-e93fee235109/
24. Garcia D.M. *et al.* (2011) Weak seed-pairing stability and high target-site abundance decrease the proficiency of lsy-6 and other microRNAs. *NSMB*. https://consensus.app/papers/details/w:66920de7-8bc1-474a-84e5-18ca1874a578/
25. Agarwal V., Bell G.W., Nam J.-W., Bartel D.P. (2015) Predicting effective microRNA target sites in mammalian mRNAs. *eLife*. https://consensus.app/papers/details/w:2cad58dc-2ca0-41fc-8e7a-3470200d57de/
26. Betel D. *et al.* (2010) Comprehensive modeling of microRNA targets predicts functional non-conserved and non-canonical sites (mirSVR). *Genome Biology*. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2945792/
27. Kertesz M. *et al.* (2007) The role of site accessibility in microRNA target recognition (PITA). *Nature Genetics*. https://members.cbio.mines-paristech.fr/~jvert/svn/bibli/local/Kertesz2007role.pdf
28. Gesell T., Washietl S. (2008) Dinucleotide controlled null models for comparative RNA gene prediction. *BMC Bioinformatics*. https://link.springer.com/doi/10.1186/1471-2105-9-248
29. McGeary S.E. *et al.* (2019) The biochemical basis of microRNA targeting efficacy. *Science*. https://pmc.ncbi.nlm.nih.gov/articles/PMC7051167
29-folding. Clote P. *et al.* Structural RNA has lower folding energy than random RNA of the same dinucleotide frequency. https://people.scs.carleton.ca/~kranakis/Papers/Zscores.pdf (and Rivas E., Eddy S.R. (2000) Secondary structure alone is generally not statistically significant for the detection of noncoding RNAs, *Bioinformatics*, cited therein)
30. Soutschek M. *et al.* (2022) scanMiR: a biochemically based toolkit for versatile and efficient microRNA target prediction. *Bioinformatics*. https://consensus.app/papers/details/w:4055f856-3ee7-433d-b4e3-fac8406fb81d/
31. Kettles G.J. *et al.* (2019) sRNA profiling combined with gene function analysis reveals a lack of evidence for cross-kingdom RNAi in the wheat–*Zymoseptoria tritici* pathosystem. *Front Plant Sci*. https://consensus.app/papers/details/w:e1b0959a-f0d2-4127-9276-b78b82b204eb/
32. Qin S. *et al.* (2022) Molecular characterization reveals no functional evidence for naturally occurring cross-kingdom RNA interference in the early stages of *Botrytis cinerea*–tomato interaction. *Mol Plant Pathol*. https://consensus.app/papers/details/w:acbf9b72-fc3f-49db-91a2-89dc2e841073/ — and Sammut S. *et al.* (2025) miRBench: novel benchmark datasets for microRNA binding site prediction that mitigate against prevalent microRNA frequency class bias. *Bioinformatics*. https://consensus.app/papers/details/w:ddea3411-267a-4d8d-8da7-f74f36422c70/
33. Pinzón N. *et al.* (2016) microRNA target prediction programs predict many false positives. *Genome Research*. https://consensus.app/papers/details/w:f23725cf-9465-4e67-ae33-3bc690ce7f2d/
34. Helwak A., Kudla G., Dudnakova T., Tollervey D. (2013) Mapping the human miRNA interactome by CLASH reveals frequent noncanonical binding. *Cell*. https://core.ac.uk/download/28971533.pdf
35. Arvey A. *et al.* (2010) Target mRNA abundance dilutes microRNA and siRNA activity. *Mol Syst Biol*. https://consensus.app/papers/details/w:344ac66d-fa9f-4395-b938-d252843e283b/
36. Khan A.A. *et al.* (2009) Transfection of small RNAs globally perturbs gene regulation by endogenous microRNAs. *Nature Biotechnology*. https://consensus.app/papers/details/w:... ; Denzler R. *et al.* (2014) *Molecular Cell*
37. Dai X., Zhao P.X. (2011) psRNATarget: a plant small RNA target analysis server. *NAR*. https://consensus.app/papers/details/w:... (NAR 2011)
38. Dai X., Zhuang Z., Zhao P.X. (2018) psRNATarget: a plant small RNA target analysis server (2017 release). *NAR*. https://consensus.app/papers/details/w:85465c20-fcc0-4d90-a284-b4c9728a0155/
39. Kusch S. *et al.* (2023) Site-specific analysis reveals candidate cross-kingdom small RNAs, tRNA and rRNA fragments, and signs of fungal RNA phasing in the barley–powdery mildew interaction. *Mol Plant Pathol*. https://www.biorxiv.org/content/10.1101/2022.07.26.501657v2.full.pdf
40. Bonnet E. *et al.* (2010) TAPIR, a web server for the prediction of plant microRNA targets, including target mimics. *Bioinformatics*. https://consensus.app/papers/details/w:566119a8-a9d6-45ab-bfed-f2934b94ae7a/
41. Srivastava P.K. *et al.* (2014) A comparison of performance of plant miRNA target prediction tools and the characterization of features for genome-wide target prediction. *BMC Genomics*. https://consensus.app/papers/details/w:0e5fa73b-6a41-48d7-9c17-18c54edb4006/
42. Ajila V. *et al.* (2023) P-TarPmiR accurately predicts plant-specific miRNA targets. *Scientific Reports*. https://consensus.app/papers/details/w:23312f3a-6650-4ab9-b03a-1801f58687dd/
43. Zanini S.F. *et al.* (2018) A bioinformatics pipeline for the analysis and target prediction of RNA effectors in bidirectional communication during plant–microbe interactions. *Front Plant Sci*. https://consensus.app/papers/details/w:13ea3671-00db-4a34-aabd-be165a029eb6/
44. Weiberg A. *et al.* (2013) Fungal small RNAs suppress plant immunity by hijacking host RNA interference pathways. *Science*.
45. Cai Q. *et al.* (2018) Plants send small RNAs in extracellular vesicles to fungal pathogen to silence virulence genes. *Science*.
46. Wong-Bajracharya J. *et al.* (2022) The ectomycorrhizal fungus *Pisolithus microcarpus* encodes a microRNA involved in cross-kingdom gene silencing during symbiosis. *PNAS*.
47. Silvestri A. *et al.* (2024) A fungal sRNA silences a host plant transcription factor to promote arbuscular mycorrhizal symbiosis. *New Phytologist*. — and Silvestri A. *et al.* (2019) *BMC Genomics*, https://consensus.app/papers/details/w:1032ac44-0488-4ccb-ba80-607c65ef46d0/
48. Johnson N.R. *et al.* (2019) Compensatory sequence variation between trans-species small RNAs and their target sites. *eLife*. https://consensus.app/papers/details/w:9b51a466-f999-41ac-85b3-ed66ad7cf5b1/
49. Mayoral J.G. *et al.* (2014) *Wolbachia* small noncoding RNAs and their role in cross-kingdom communications. *PNAS*. https://pmc.ncbi.nlm.nih.gov/articles/PMC4284532
50. RNAhybrid / RNAcalibrate manual pages. https://manpages.debian.org/wheezy/rnahybrid/RNAhybrid.1.en.html ; https://manpages.debian.org/wheezy/rnahybrid/RNAcalibrate.1.en.html
51. Rehmsmeier M. *et al.* (2004) Fast and effective prediction of microRNA/target duplexes. *RNA*; Krüger J., Rehmsmeier M. (2006) RNAhybrid: microRNA target prediction easy, fast and flexible. *NAR*. https://consensus.app/papers/details/w:... 
52. Melamed S. *et al.* (2016) Global mapping of small RNA-target interactions in bacteria (RIL-seq). *Molecular Cell*. https://consensus.app/papers/details/w:118fe378-ee19-40ef-8c65-f7f9f9d4bc9b/
53. King A. *et al.* (2019) sRNA Target Prediction Organizing Tool (SPOT). *mSphere*. https://consensus.app/papers/details/w:... 
54. Min S. *et al.* (2021) TargetNet: functional microRNA target prediction with deep neural networks. *Bioinformatics*. https://consensus.app/papers/details/w:2470bb60-15fc-4a79-97b4-1fda28d581dd/

---

## 9. Recommendations for SmithHunter2 module B

A concrete, opinionated design. Every value has a stated reason.

### 9.1 Architecture: three independent evidence channels, never one score

Produce, per (sRNA, transcript) pair, three scores that are **reported separately and combined only at the reporting layer**:

- **Channel S — site class.** Best seed-class call among 8mer / 7mer-m8 / 7mer-A1 / 6mer / offset-6mer at guide positions 2–8 (already built).
- **Channel E — accessibility-aware energy.** IntaRNA as primary; RNAhybrid as an independent structure-naive second opinion.
- **Channel F — folding annotation.** MFE, MFEI, dinucleotide z-score of the sRNA locus. **Annotation only, never a filter.**

Rationale: §4. Seed scoring and thermodynamic scoring fail in opposite directions, and when the mechanism is unknown the agreement between them *is* the information. Collapsing them into a weighted sum hides exactly what the reader needs to judge.

### 9.2 Default (seed) mode — concrete parameters

**Search region.** 3'UTR + full CDS for the host; **do not restrict to 3'UTR**. CLASH put ~50% of miRNA sites in CDS [34], and the one well-characterised trans-species system targets conserved coding regions almost exclusively [48]. Restricting to 3'UTR would make a true cross-kingdom effector invisible by construction.

**Seed pinning: mandatory.** Guide positions 2–8 only. Scanning all sRNA windows multiplies the chance-hit base rate by (L−6) — ~94× for a 100-nt sRNA (§6.1). Expose an `--all-windows` escape hatch, off by default, with a loud warning in the log stating the expected chance-hit count.

**Drop the A1 requirement to a flag, don't enforce it.** The A1 anchor is a property of the AGO MID pocket, not of base pairing. Emit 7mer-A1 and 8mer calls, but also emit the A1-agnostic 7mer-m8 / 7mer equivalent, and carry a boolean `a1_anchor_present`. Let the downstream analysis decide whether AGO loading is assumed.

**3' supplementary pairing.** Compute pairing at guide 13–16 and store as `supp_bp_13_16` (count of WC pairs). This is the one context feature with a generic base-pairing rationale, it was the most conserved and most repressive CLASH class [34], and it costs almost nothing. **Do not weight it into the score** — just report it.

**Explicitly do NOT implement:** 3'UTR AU content, distance from stop codon, position within UTR, target-site abundance / AGO competition. These are quantitative AGO/UTR-biology constants fitted in human, mouse and fly cells [23,24,25]. There is no basis for their coefficients in a bacterium-to-insect setting, and importing them would manufacture false confidence. State this explicitly in the methods.

### 9.3 Bacterial mode — concrete parameters

**Search region:** −200 / +100 nt around each start codon. Justification: this is the field-standard window and Pain et al. showed replacing it with RNA-seq-derived real UTRs gives only 4.9–17.5% AUC change, none significant [1]. Do not spend engineering effort on exact UTR boundaries.

**No seed rule.** Correct — TargetRNA3's SHAP analysis found seed length 8–10 contributed little in the largest bacterial feature study to date [7].

**IntaRNA parameters** (from [5], the only systematic sweep; it recovered ~30% more verified interactions at ~20% of baseline runtime):

```
--seedBP 7
--seedNoGU
--seedMaxE -2          # no-GU variant of the tuned optimum
--seedMinPu 0.001
--intLenMax 60
--intLoopMax 8
--outMaxE -10
--qAcc C --qAccW 0 --qAccL 0        # fold the whole sRNA: its structure is real
--tAcc C --tAccW 150 --tAccL 100    # local window on the target; global folding of long
                                    # mRNA is slow and unreliable [20]
--mode H
--outMode C --outNumber 1 --outOverlap N   # one best interaction per transcript, for FDR
--threads <n>
```

Add an optional second run with the **Andronescu-2007** energy model (+5 verified pairs over Turner-2004 in [5]) and report hits found by both parameterisations.

**Add Beisel's refinement as an opt-in:** restrict `--qRegion` to predicted single-stranded regions of the sRNA. Using only Spot 42's three unstructured regions substantially improved recovery of known and novel targets [17]. This is cheap and is the closest thing the bacterial field has to a free lunch after parameter tuning.

**RNAhybrid as second opinion:** `-f 2,8` in seed mode (this is how you pin a seed in RNAhybrid), `-n` raised above the sRNA length (**default is 30 and long sRNAs silently fail**), `-e` with a generous `-b` so the hit cap does not truncate, and **`-d ξ,θ` from `RNAcalibrate` run on your own target set** — never `-s 3utr_human/fly/worm` [50].

### 9.4 Empirical FDR — concrete design

**Two nulls. Report both q-values. Disagreement is a demotion signal.**

**Null A (primary):** 500 Altschul–Erikson **dinucleotide-preserving** shuffles per sRNA, each scored against the **unshuffled** target set through the identical pipeline.
- 500, not 100, because at target q = 0.01 you need enough decoy mass in the tail; 500 decoys × ~4,000 transcripts = 2×10⁶ null scores, which resolves q to ~5×10⁻⁷.
- Dinucleotide, not mononucleotide: stacking free energies are dinucleotide quantities, and the folding literature explicitly says mononucleotide shuffles must not be used for z-scores [28, 29-folding]. Mononucleotide shuffling makes real sequences look spuriously significant.
- **Shuffle the whole sRNA including positions 2–8.** Holding the seed fixed collapses the null for the feature under test.

**Null B (orthogonal):** 200+ real sRNAs of matched length from a phylogenetically and ecologically distant bacterium with no shared host. This controls for everything a shuffle destroys: real folding, real k-mer vocabulary, real 5'-end nucleotide bias. Document the donor organism and acknowledge the conservative bias if the donor could plausibly target the same host.

**Never use shuffled targets as the primary null** — it alters the search space rather than the query and destroys real UTR/CDS composition structure, giving an anti-conservative FDR.

**Permuted sRNA–target assignment** is an evaluation-only control, with the miRBench frequency-class-bias caveat [32] stated in the docs. Do not use it to generate ML training negatives.

**Stratify.** Compute q separately per (site class × transcript-length decile × GC quartile). A pooled q-value hides that the 6mer stratum has ~293 expected chance hits per sRNA in the bacterial window versus ~18 for 8mer (§6.1), and that long GC-rich transcripts win on energy for reasons that have nothing to do with biology.

**Report q, not p.** Storey–Tibshirani q-values from the empirical null. Surface the tool-reported p-value alongside, clearly labelled as uncalibrated.

### 9.5 What to report and how to frame it

Per candidate: site class, A1 flag, `supp_bp_13_16`, IntaRNA ΔG_total with its ΔG_hybrid / ΔG_unfold decomposition, RNAhybrid normalised energy with recalibrated p, q_A, q_B, mode flag, and a `both_channels` boolean (seed class ≥ 7mer **and** q_E < 0.1).

**Headline the `both_channels` set.** It is small, it is the only set with two independent lines of evidence, and it is the honest deliverable.

**State the base rate in the output header**, computed for the actual run: "search space L nt; expected chance 7mer matches ≈ L/16384 = N." A reader who sees N = 73 next to 91 reported hits will calibrate correctly without being told.

**State plainly in the methods** that no validated method exists for bacterial-sRNA → animal-host-mRNA prediction, that module B is a calibrated ranked hypothesis generator, and that absolute precision is unknown and on the evidence from better-resourced fields is likely below 50% (§6.3). This is not hedging; it is the difference between this tool and the 3,185-predicted-targets-no-phenotype literature [32].

### 9.6 Implementation priority order

1. **Tune the IntaRNA parameters** (§9.3). Biggest single performance gain available, costs a config file [5].
2. **Seed pinning + the base-rate line in the output header.** Two orders of magnitude of false positives, removed by a default [§6.1].
3. **Null A (dinucleotide shuffle q-values), stratified.** The core scientific contribution over existing cross-kingdom practice.
4. **RNAcalibrate wiring for RNAhybrid.** Small, and it makes the second opinion mean something.
5. **Null B (unrelated-organism decoys).** The orthogonal check.
6. **CDS search region + `supp_bp_13_16`.** Cheap, and required if the mechanism turns out to be the *Cuscuta* type [48].
7. **Beisel single-stranded `--qRegion` restriction.** Opt-in [17].
8. **Folding annotation fields.** Last, and explicitly non-filtering.

**Deliberately not on this list:** TargetScan context++ features, a trained ML model, conservation/comparative scoring. The first has no mechanistic warrant here; the second has no training data for this problem (the whole bacterial field has ~4,386 pairs, none cross-kingdom [7]); the third is structurally unavailable given one bacterial genome per library — and it is precisely the signal that made CopraRNA the best tool in family 1. Module B should say openly that it is working without the field's strongest feature, rather than hoping no one notices.
