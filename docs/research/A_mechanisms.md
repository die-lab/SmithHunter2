# Strand A — MECHANISM: how cross-kingdom RNAi actually works, at the resolution a target-prediction algorithm needs

Literature review for **SmithHunter2** (detection of endosymbiont-derived small RNAs that may target host mRNAs).
Compiled 2026-11-13. Sources retrieved via Consensus (≈600 records screened across 10 broad searches) plus targeted web retrieval of primary full text where numbers were needed.

**Evidence grading used throughout:**
- **[D] Demonstrated** — causal, replicated or multi-assay (AGO-IP + mutant + reporter + phenotype), ideally by more than one group.
- **[R] Reported once / unreplicated** — a single paper, usually correlative plus one functional assay.
- **[S] Speculation / extrapolation** — plausible but not tested in the stated system.
- **[C] Contested** — a direct negative or failed-replication result exists.

A deliberate warning up front: the mechanistic detail that a target predictor would most like to have — a quantitative pairing rule learned from a set of validated cross-kingdom pairs — **does not exist in the literature in usable form for bacteria**, and exists only in a weak, partly contested form for fungi. Most of this report is therefore about constraining the search space and bounding the false-positive rate, not about optimising a scoring function.

---

## 1. Pairing rules

### 1.1 Plant hosts: plant-like, near-complementary, cleavage-competent pairing — not animal seed pairing

The best-characterised system is *Botrytis cinerea* → *Arabidopsis*/tomato. Weiberg et al. 2013 (Science) filtered 832 candidate Bc-sRNAs (present in both infected *Arabidopsis* and tomato libraries, enriched over axenic culture, exact genome matches) down to **73** predicted to target host genes in both hosts under stringent criteria; **63 of those 73 were 20–22 nt with a 5′-terminal U**, the configuration favoured for AtAGO1 loading (Weiberg et al. 2013, Science). The three functionally characterised effectors — Bc-siR3.1, Bc-siR3.2, Bc-siR5 — are **21 nt**; Bc-siR3.1/3.2 derive from the same LTR locus offset by a 4-nt shift. Targets: Bc-siR3.2 → At-MPK1/MPK2 (and predicted Sl-MAPKKK4), Bc-siR3.1 → peroxiredoxin PRXIIF, Bc-siR5 → a cell-wall-associated kinase WAK. **[D]**

Mechanistically this is **plant-type siRNA pairing**, i.e. extensive near-complementarity over most of the guide, with cleavage, and the prediction was done with plant target-prediction logic (penalty-score schemes of the psRNATarget/TargetFinder family), **not** with a 7-mer/8-mer seed rule. Cross-kingdom degradomics in *Brachypodium distachyon*–*Fusarium graminearum* recovered **258 ckRNAi-mediating sRNAs** by joint sRNA-seq + transcriptome + degradome sequencing and then used simulation to sort them into "functionality classes" by sequence characteristics; the same study reports that plant ck-sRNAs are organised in long non-coding RNAs and that ck-sRNAs show evolution toward invariant mRNA target regions (Werner et al. 2025, bioRxiv). **[R]** — this is, to my knowledge, the only genome-wide attempt to learn pairing-feature classes from validated cross-kingdom pairs, and it is a preprint.

The oomycete case is the same shape: *Hyaloperonospora arabidopsidis* sRNAs act through AtAGO1/RISC, demonstrated with a purpose-built CRISPR-Csy4/GUS *in situ* reporter, with *atago1* mutants and short-tandem-target-mimic (STTM) plants both showing enhanced resistance, and knockouts of the targeted plant genes showing quantitatively enhanced susceptibility (Dunker et al. 2020, eLife). **[D]**

A single AMF example gives a cleavage-level answer: *Rhizophagus irregularis* Rir2216 is loaded into *Medicago truncatula* AGO1 and **cleaves** MtWRKY69, shown by 5′ RACE plus AGO1-IP plus heterologous co-expression (Silvestri et al. 2024, New Phytologist). **[R]**

**Practical read for a predictor targeting plant hosts:** near-complementary, centrally paired, cleavage-competent sites; plant-style penalty scoring; no categorical requirement for a canonical 2–8 seed but strong pairing in the 5′ half is what the plant scoring schemes effectively enforce.

### 1.2 Animal hosts: the few validated pairs look miRNA-like, with a conserved seed, and most are mimicry

- *Fasciola hepatica* fhe-miR-125b is **homologous to mammalian hsa-miR-125b** and loads onto host AGO2 in macrophages during infection; the authors frame it explicitly as mimicry of a host miRNA (Tran et al. 2021, Scientific Reports). **[R]**
- Nematode miR-5352 is conserved across GI nematodes and shows **seed-sequence conservation with the mammalian miR-92a family**, which the authors attribute to convergent evolution onto an existing host regulatory network; it represses Klf-4 and the IL-22 receptor (Pérez et al. 2025, Frontiers in Immunology). **[R]**
- Helminth EV let-7 is highly conserved and **highly homologous in the host**; its function in crosstalk is argued to reflect active selection for compatibility with host miRNA machinery (Zhong et al. 2024, Frontiers in Immunology). **[S/R]**
- Conversely, Hpo-miR-71-5p is **absent in mammals** and still validates against host *Irf4*, and Asu-miR-791-3p against *Mtor*; the authors are candid that "the absolute repression is subtle" and propose synergy across the dozens of miRNAs released (Soichot et al. 2022, Frontiers in Molecular Biosciences). **[R]**
- A systematic comparison found that parasite miRNA **seed sequences are more often identical to those of *C. elegans* than to those of their hosts**, but that the *secreted* fraction shares seeds with host miRNAs more often than the non-secreted fraction (Duguet et al. 2020, Frontiers in Cellular and Infection Microbiology). **[R]** — this is the single most useful quantitative-ish statement in the parasite literature for a predictor: selection for host-seed-compatibility is detectable but weak.

The one true cross-kingdom AGO-RIP in a *plant* host for a *metazoan* parasite: AGO1-RIP from tomato galls showed selective *in planta* loading of **10 *Meloidogyne incognita* miRNAs** into host AGO1; degradome + target prediction + dual-luciferase validated **9 tomato transcripts silenced by 4 nematode miRNAs**; AGO1-RIP in *Arabidopsis* recovered a conserved subset including miR-2b and **miR-100** (Dussutour et al. 2026, bioRxiv). **[R]**

### 1.3 Bacteria → eukaryote: pairing is essentially unconstrained in the literature

This is the weakest layer and the one most relevant to SmithHunter2. Across the bacterial examples I found, **every reported pair was obtained by running a eukaryotic miRNA target-prediction tool (miRanda / RNAhybrid / TargetScan-like) on a bacterial sRNA and then testing one or a few hits in a reporter assay.** No paper derives a pairing rule from data, and no paper reports a seed-class distribution, mismatch distribution, or free-energy distribution across a set of validated bacterial→host pairs.

Concrete instances:
- *Pseudomonas aeruginosa* OMV sRNA52320, a **fragment of a methionine tRNA**, reduces LPS- and OMV-induced IL-8 in human airway epithelium and KC/neutrophil infiltration in mouse lung; selected "based on its stable secondary structure and predicted mRNA targets" (Koeppen et al. 2016, PLoS Pathogens). **[D]** for the phenotype, **[S]** for the pairing.
- *Legionella pneumophila* RsmY binds the **UTR of *ddx58*** and *cRel*; tRNA-Phe binds *ddx58* and *irak1*; the authors call them "miRNA-like" trans-kingdom regulators (Sahr et al. 2022, Nature Communications). **[R]**
- *Porphyromonas gingivalis* OMV msRNA sRNA45033 targets the **3′ UTR of CBX5**, validated by dual-luciferase plus overexpression/knockdown (Montes-García et al. 2026, Archives of Microbiology, reviewing the primary work). **[R]** — this is repeatedly described as the most robust bacterial→host mechanistic example.
- *Lactobacillus plantarum* EV sRNA71 suppresses Tp53 by binding its **3′ UTR**, validated by luciferase (Yu et al. 2022, Frontiers in Microbiology). **[R]**
- *Escherichia coli* Ile-tRF-5X is released in OMVs, transferred into HCT116 cells, and **promotes** MAP3K4 expression (Diallo et al. 2022, PLoS Pathogens). **[R]** — note the direction: up, not down. A silencing-only predictor would never have proposed this.
- *Brucella melitensis* OMV sRNA t0000835 is **predicted** to target GPR45 "by binding to AGO2"; the only functional datum is that transfection increased intracellular survival (Chen et al. 2026, Kafkas Univ Vet Fak Derg). **[R]**, weak.
- *Xanthomonas oryzae* pv. *oryzicola* Xosr001 in OMVs suppresses rice *OsJMT1*, with a ΔXosr001 mutant phenotype and synthetic-sRNA application on a tagged overexpressor (Wu et al. 2024, Plant Communications). **[R]**, among the better-controlled bacterial examples.
- *Vibrio fischeri* SsrA (**tmRNA**) is loaded into OMVs and delivered into squid light-organ epithelium; a ΔssrA mutant produces an abnormal host transcriptional response (10 genes differentially up-regulated) and reduced host robustness (Moriano-Gutierrez et al. 2020, PLoS Biology). **[D]** for the symbiotic phenotype; **critically, the mechanism here is not demonstrated to be base-pairing at all** — it is equally compatible with innate-immune RNA sensing.

The closest thing to a bacterial→eukaryote pairing rule with real machinery behind it is **rhizobial tRFs hijacking soybean AGO1**: three tRF families target nodulation genes through host AGO1, with tRF silencing and target overexpression both reducing nodulation, and target repression by amiRNA/CRISPR promoting it (Ren et al. 2019, Science). A follow-up found other rhizobial tRFs co-immunoprecipitated with *Phaseolus vulgaris* AGO5 (reviewed in Ledford et al. 2023, New Phytologist). **[D]** — and note it is a *plant* host, so it inherits plant pairing rules, not animal seed rules.

### 1.4 What the numbers actually are

| Quantity | Best available value | Source | Grade |
|---|---|---|---|
| Fraction of candidate fungal ck-sRNAs that are 20–22 nt with 5′U | 63/73 (86%) | Weiberg 2013 | D |
| Candidate pool → predicted host-targeting, *B. cinerea* | 832 → 73 (8.8%) | Weiberg 2013 | D |
| Predicted *B. cinerea* sRNA → tomato mRNA pairs, unfiltered | 7,042 sRNAs → 3,185 mRNAs; **163** down-regulated in planta (5.1%) | Qin 2022 | C |
| Predicted plant → fungus coverage, unfiltered | ~1 M unique tomato sRNAs could target **97% of all fungal genes** | Qin 2023 | C |
| ckRNAi-mediating sRNAs recovered by degradome, *Fg*–*Bd* | 258 | Werner 2025 | R |
| AGO1-RIP-confirmed nematode miRNAs in tomato | 10 loaded; 9 targets validated by 4 miRNAs | Dussutour 2026 | R |
| HITS-CLIP miRNA–target map, *S. japonicum* (intra-parasite) | 769 miRNAs, 11,854 target sites, ~15 sites/miRNA, **50% verification rate** on a small test set | Zhao 2015 | R |
| Off-target silencing threshold in plants | as few as **11 contiguous nt**, or **15/19 bp**, suffices | reviewed in Ledford 2023 / SIGS design guidance | D |

That last row is the one that should frighten anyone building a predictor. **The 97% figure from Qin et al. 2023 is the key negative statistic in the whole field**: with ~10⁶ sRNA species and a few-mismatch-tolerant pairing rule, essentially every transcript in the other organism acquires a predicted target site. Sequence complementarity alone carries almost no information.

---

## 2. Which silencing machinery

### 2.1 Plant hosts

**AGO1 is the near-universal answer.** Demonstrated by AGO1 co-IP for *B. cinerea* in *Arabidopsis* (Weiberg 2013) and in tomato (He 2023, bioRxiv); by *atago1* phenotype + STTM for *H. arabidopsidis* (Dunker 2020); by AGO1-IP for *Serendipita indica* (Nasfi 2024, J Exp Bot); by AGO1-IP for *Verticillium dahliae* VdrsR-1 (Zhang 2022, Front Plant Sci); by AGO1-IP + 5′RACE for *R. irregularis* Rir2216 (Silvestri 2024); by OsAGO1 binding for *Rhizoctonia solani* (Ni 2026, Nat Commun); by AGO1-RIP for *M. incognita* (Dussutour 2026); by AGO1 for rhizobial tRFs (Ren 2019). In the *B. cinerea* case the effectors were present in the AGO1 fraction and **not** in the AGO2- or AGO4-associated fractions (Weiberg 2013). **[D]**

Exceptions worth noting:
- *Fusarium oxysporum* f. sp. *lycopersici* Fol-milR1 binds **SlyAGO4a**, and *SlyAGO4a* VIGS knockdown reduces susceptibility (Ji et al. 2021, New Phytologist). **[R]** — a genuine deviation from the AGO1 rule.
- A systematic *Arabidopsis* DCL/AGO mutant panel across many microbes found AGO1 acts as a **positive or negative** regulator of immunity depending on the pathogen, AGO4 as a positive regulator, and AGO2/AGO10 as modulators of colonisation by a mutualist — with several interactions showing wild-type-like phenotypes in all mutants (Ruf et al. 2024, Plant Direct). **[D]** — i.e. the AGO1-hijack model does not generalise across pathosystems.
- The reverse direction uses a different AGO: AGO10 is required for pathogen-induced trans-species silencing in *Arabidopsis* (Wang et al. 2026, bioRxiv). **[R]**

**Slicing vs repression:** in plants the demonstrated mode is **slicer-dependent cleavage** — 5′ RACE cleavage products for Rir2216–MtWRKY69 (Silvestri 2024) and infection-specific degradation products by RLM-RACE for three barley transcripts targeted by *F. graminearum* tRFs (Werner et al. 2021, PLoS ONE). **[D]**

**Amplification matters and is often ignored.** Endophyte-induced silencing in *N. benthamiana* required **RDR6**, produced secondary sRNAs mainly 22–24 nt with a conspicuous absence of 21-nt species, and failed in an intron-containing target (Kellari et al. 2025, Plant Biology). **[R]** Conversely, HIGS into *V. dahliae* showed **no transitivity** — siRNAs occurred exclusively within the dsRNA region (Zhang et al. 2022, IJMS). **[D]** So secondary amplification is host- and system-specific and cannot be assumed.

A fungal-side twist: BcAGO1 binds both fungal and plant sRNAs and acts bidirectionally, and BcAGO2 is required for efficient **delivery** of pathogen sRNAs into host cells (Weiberg et al. 2025, PNAS). **[R]**

### 2.2 Animal hosts

Much thinner. The positive evidence:
- fhe-miR-125b loads onto host **Ago-2** in macrophages (Tran 2021). **[R]**
- Hpo-miR-71-5p "can incorporate into mouse Argonaute protein complexes" and transfected miRNAs persist for days in recipient cells (Soichot 2022). **[R]**
- *Beauveria bassiana* bba-milR1 binds mosquito **AGO1** and silences *Spz4* (Cui et al. 2019, Nature Communications). **[D]** — this is probably the single best-controlled fungus→animal cross-kingdom RNAi paper: expression timing, target silencing, and a biphasic virulence rationale.
- Helminth EVs carry a nematode **Argonaute** of their own: exWAGO is secreted at multiple copies per EV and associates with **secondary siRNAs produced by RdRPs from transposons and newly evolved repeats** — not miRNAs, as had been assumed (Chow et al. 2019, Nucleic Acids Research). **[D]** This is important: it means a parasite can in principle deliver a pre-loaded RISC rather than relying on host AGO loading.

No bacterial sRNA has been shown by AGO-IP or CLIP to be loaded into a host Argonaute. The Brucella paper merely *predicts* AGO2 binding (Chen 2026). **[gap]**

**Plant vs animal, summarised:** plant hosts load 21–22 nt 5′U guides into AGO1 and slice; animal hosts would have to load a guide into AGO1/AGO2 and act through seed-driven, mostly non-slicing repression. The piRNA pathway is perturbed by *Wolbachia* (Mayoral et al. 2014, PLoS ONE; Dou & Li 2026, Front Microbiol) but there is **no evidence that bacterial sRNAs enter piRNA pathways**. **[gap]**

---

## 3. Length and 5′ nucleotide

**Fungi/oomycetes → plants.** Functional effectors are **20–22 nt, overwhelmingly 21 nt, with 5′U** (Weiberg 2013). 5′U enrichment is described as "a common feature of RNAi-competent fungi" and is the basis of the (hypothesised) ability to hijack plant AGO1 (Ledford et al. 2023, New Phytologist). **[D for Botrytis, S as a general rule.]**

Important counter-data on length:
- AMF (*R. irregularis*, *Gigaspora margarita*) sRNAs have a **unimodal length distribution centred at 24, 25 or 26 nt**, not 21 (Ledford 2023) — yet Rir2216 still loads into AGO1 (Silvestri 2024). The 21-nt rule is therefore not necessary.
- *Colletotrichum higginsianum* sRNAs are dominated by **29 nt** species in appressorial/biotrophic/mycelial stages and shift to **18 nt** in necrotrophy, which the authors attribute to RNA degradation during host cell death (Armijos et al. 2026, Plant Direct). **[R]** This is a direct warning: length distributions in infected tissue are confounded by degradation state.
- Oomycete sRNAs cluster at **21, 25 and 26 nt with species-specific 5′-base preferences** (noted in the Consensus corpus).
- *V. dahliae* VdrsR-1 is **24 nt and rRNA-derived**, and *increases* miR157d accumulation by targeting the MIR157d precursor (Zhang 2022). **[R]** — length 24, non-canonical mechanism, net effect opposite in sign.

**Bacteria.** There is no coherent size class.
- "msRNAs" in periodontal pathogens are defined as **15–25 nt** with predicted hairpin precursors (Choi et al. 2017, J Dent Res). **[R]**
- *E. coli* and five other species carry abundant **vsRNAs of 8–15 nt** that are *more* abundant than 16–30 nt species and account for about a third of reads in the 8–30 nt window; tRFs (Ile-tRF, Ala-tRF) dominate and are **selectively loaded into OMVs** (Diallo et al. 2022, Front Mol Biosci). **[R]** An 8–15 nt "guide" cannot support canonical RISC loading; whatever these do, it is not classical RNAi.
- Classical bacterial regulatory sRNAs are **50–500 nt** (Mayoral et al. 2014, PNAS) — far too long for AGO.
- *M. tuberculosis* sncRNAs produced in infected macrophages are **<30 nt** with stable predicted stem-loops and evidence of hairpin-dependent processing — but the demonstrated function is on *mycobacterial* transcripts (oleic-acid biosynthesis), not host transcripts (Coskun et al. 2020, Front Microbiol). **[R]**

**Does the 5′ nt match host AGO preference?** For fungi→plants, yes and this is the central argument (5′U → AtAGO1). For bacteria, **nobody has reported a 5′-nucleotide distribution for the sRNAs they claim target host mRNAs**. This is a conspicuous and testable gap, and one SmithHunter2 could fill. **[gap]**

---

## 4. Delivery

**Fungi → plants: EVs, demonstrated.** *B. cinerea* secretes Bc-sRNAs in EVs marked by the tetraspanin BcPLS1; these enter plant cells by **clathrin-mediated endocytosis**; CME mutants show increased resistance and **attenuated Bc-sRNA loading into AtAGO1 and attenuated target suppression** (He et al. 2023, Nature Communications). **[D]** The reverse direction uses exosome-like EVs with AGO1, RNA helicases RH11/RH37 and annexins ANN1/ANN2 as selective loading factors (Cai et al. 2018, Science; He et al. 2021, Nature Plants). **[D]** Plants also move whole **mRNAs** into fungi in EVs, and those mRNAs are translated in the recipient (Wang et al. 2023, Cell Host & Microbe) — the channel is not sRNA-specific. **[R]**

Non-vesicular routes exist too: in the *Arabidopsis* apoplast, three populations of active extracellular sRNAs were recovered — non-vesicular protein-associated, EV-contained, and **free dsRNA ("efsRNAs") unbound to protein**; both transgene- and genome-derived efsRNAs were found inside *Pseudomonas syringae* cells (Ravet et al. 2025, Nature Communications). **[D]** This is the strongest evidence that plant sRNAs reach bacterial cytoplasm.

**Bacteria → host: OMVs.** OMV delivery of RNA into host cells is demonstrated by imaging and by functional readouts across many systems (Koeppen 2016; Choi 2017; Ha 2020; Sahr 2022; Dauros-Singorenko 2018). Proof-of-principle for functional biomolecule transfer in vivo comes from *E. coli*^Cre into Rosa26.tdTomato mice, where OMV-delivered Cre induced reporter expression in intestinal epithelium, stem cells, macrophages and distal organs (Bittel et al. 2021, J Extracell Vesicles). **[D]** New structural work shows OMV RNA is not only luminal but can be **inserted into the vesicle membrane in an Hfq-dependent manner**, i.e. surface-presented (Mosca et al. 2026, J Extracell Vesicles). **[R]**

Caveats that a predictor must respect:
- Uptake efficiency is species- and cell-type-specific. dsRNA uptake is efficient in *B. cinerea*, *S. sclerotiorum*, *R. solani*, *A. niger*, *V. dahliae*; absent in *Colletotrichum gloeosporioides*; weak in *Trichoderma virens*; limited and stage-variable in *Phytophthora infestans*; and *Zymoseptoria tritici* is **incapable of dsRNA uptake** (Qiao et al. 2021, Plant Biotechnol J; Kettles et al. 2019, Front Plant Sci). **[D]**
- The confound of co-purified LPS and other non-RNA cargo is explicitly flagged as the central unresolved methodological problem in bacterial-EV RNA work (Dauros-Singorenko et al. 2018, FEMS Microbiol Lett). **[D]**
- RNA inside OMVs can act as a **TLR ligand** rather than a guide. *S. aureus* MVs stimulate IFN-β through endosomal TLRs (Rodriguez & Kuehn 2020, Sci Rep); *B. fragilis* OMVs activate TLR2/3/7 (Sheikh et al. 2025, J Extracell Biol); *A. actinomycetemcomitans* exRNA activates IL-6/NF-κB in microglia (Ha et al. 2020, Front Mol Biosci); tRF-rich circulating sncRNAs in TB patients activate macrophages via **TLR7** (Gumas et al. 2024, Mol Ther Nucleic Acids). **[D]** For every reported "bacterial sRNA silences host gene X", innate RNA sensing is a live alternative explanation.

**Intracellular symbionts.** An intracellular, vertically transmitted bacterium is not behind a vesicle barrier in the same way — but it *is* usually inside a host-derived vacuolar membrane. The literature does not resolve how *Wolbachia* RNA would cross that membrane. **[gap]** The only directional evidence in an obligate endosymbiosis runs the *other* way: aphid miR-1, miR-10, miR-29 and miR-927 localise **inside *Buchnera*** cells by miRNA-FISH, with aphid Ago1 also translocated (Ding et al. 2026, bioRxiv). **[R]**

**Does the delivery route constrain candidacy?** Yes, strongly, and this is one of the few hard filters available: a candidate must be demonstrably *exported*, i.e. enriched in a vesicle or extracellular fraction relative to the whole-cell transcriptome. Selective enrichment in OMVs relative to the cytosol is reported repeatedly (Koeppen 2016; Chapagain 2024; Diallo 2022; Naaz 2026) and OMV sncRNA abundance is normally **~1/10 of cytosolic** (Naaz et al. 2026, Non-Coding RNA), so export is selective and measurable.

---

## 5. Target region: 3′UTR vs CDS vs 5′UTR

**Animal hosts.** Where a region is reported, it is **3′UTR**, consistent with canonical animal miRNA biology: P. gingivalis sRNA45033 → CBX5 3′UTR; *L. plantarum* sRNA71 → Tp53 3′UTR; Legionella RsmY → *ddx58* UTR; *Wolbachia*-induced host miRNAs → *psq* 3′UTR (Zheng et al. 2019, J Insect Physiol) and → *crebA* 3′UTR (Bi et al. 2018, IBMB). Reviews of parasite miRNAs state the binding is "most often in the 3′UTR" (Britton et al. 2015, Vet Parasitol; Rojas-Pirela et al. 2022, Open Biology). **[R, consistently]**

But the honest statement is: **this is almost certainly reporting bias, not a measured enrichment.** Every one of these studies ran a 3′UTR-restricted prediction tool (miRanda/TargetScan defaults) and then tested a hit in a 3′UTR luciferase construct. **No cross-kingdom study has tested CDS and 5′UTR sites on equal footing and reported the relative enrichment.** A 3′UTR prior for cross-kingdom pairs is therefore circular. **[gap]**

One concrete warning against the 3′UTR prior as a scoring feature: a species-matched camel-forage analysis found that predicted targeting multiplicity **scales with 3′UTR length** (Pearson r = 0.63, Spearman ρ = 0.58), so apparent "hub" targets are largely long-3′UTR genes (Zoziuk et al. 2026, Genes). **[D]** Any predictor scoring 3′UTRs must normalise for UTR length or it will simply rediscover long UTRs.

**Plant hosts.** The question barely arises: plant targeting is by near-complementary sites that can lie anywhere in the transcript, most validated sites are in the CDS, and the demonstrated outcome is cleavage. The *V. dahliae* case targets a **miRNA precursor** (MIR157d), not an mRNA at all (Zhang 2022) — a reminder that non-mRNA targets are possible. Fungal sRNAs also target host **transposon- and pseudogene-derived** loci and phased loci (Kusch et al. 2023, Mol Plant Pathol). **[R]**

---

## 6. Processing: Dicer products vs degradation fragments

### 6.1 Fungi: genuinely Dicer-dependent, and genuinely transposon-derived

*B. cinerea* effector sRNAs are **DCL1/DCL2-dependent**: the *dcl1 dcl2* double mutant loses them and shows reduced pathogenicity (Weiberg 2013; Wang et al. 2016, Nature Plants). **[D]** They are overwhelmingly **retrotransposon-derived**: of the 73 host-targeting Bc-sRNAs, **52 came from six LTR loci**, 13 from intergenic regions, 8 from protein-coding genes (Weiberg 2013). Less aggressive *B. cinerea* strains have lower LTR copy numbers and have lost retrotransposon-derived BcsRNA production; reconstituting retrotransposon expression restores aggressiveness in a BcsRNA-dependent way (Porquier et al. 2021, Genome Biology). **[D]** Production also requires **BcRDR1** (Cheng et al. 2023, PLoS Pathogens). **[R]**

### 6.2 But a large and growing fraction of reported cross-kingdom sRNAs are tRFs and rRFs

- *F. graminearum*: the three functionally implicated sRNAs targeting HvEOL1, HvBAK1, HvSERK2 and BdSERK2 are **tRNA-derived fragments**, and their abundance is DCL-dependent (Werner et al. 2021, PLoS ONE). **[R]**
- *V. dahliae* VdrsR-1 is **rRNA-derived**, 24 nt (Zhang 2022). **[R]**
- Barley–*Blumeria hordei*: 31–33-base **5′ fragments of barley 5.8S rRNA** enriched in EVs and infected epidermis, and *B. hordei* **tRFs** in haustoria (Kusch et al. 2023). **[R]**
- Bacteria: *P. aeruginosa* sRNA52320 is a **Met-tRNA fragment** (Koeppen 2016); *E. coli* Ile-tRF-5X is a tRF (Diallo 2022); *Legionella* uses **tRNA-Phe** (Sahr 2022); *V. fischeri* uses **tmRNA/SsrA** (Moriano-Gutierrez 2020). A dedicated review concludes tRFs are the underappreciated regulatory class in microbial pathogenesis (Li & Stanton 2021, Front Microbiol).
- Salmonella OMV RNA is dominated by **rRNA** in most conditions, with some transcripts exported as full-length and others in processed/degraded form (Malabirade et al. 2018, Front Microbiol). **[D]**

**Rough fraction, bacterial side:** of the ~8 bacterial→host "sRNA effector" claims with a named molecule that I could find, **at least 5 are tRF-, rRF- or tmRNA-derived** (P. aeruginosa, E. coli, Legionella ×2, V. fischeri), plus the rhizobial tRFs in plants. That is ≥60%, and the OMV RNA pool as a whole is rRNA-dominated.

**Does this affect credibility?** Two ways, pulling in opposite directions:
1. **Against:** tRNA and rRNA fragments are exactly what you get from RNA degradation, they are the most abundant species in any library, they are the species most likely to be present as environmental/reagent contamination, and their sequences are the most conserved across the tree of life — so they are the most likely to align spuriously to host transcripts. Witwer 2018 (F1000Research) showed that the one putative plant miRNA that mapped consistently above background in human plasma, peu-MIR2910, is **100% identical to a human rRNA sequence**. Degradation-fragment identity is, empirically, the signature of an artefact.
2. **For:** tRFs and rRFs are *selectively* loaded into OMVs relative to the cytosol (Diallo 2022; Chapagain 2024), which is hard to explain by degradation alone, and the rhizobial tRF case has full genetic support (Ren 2019).

**Operational conclusion for SmithHunter2: tRF/rRF origin should be a flag, not an exclusion — but a candidate whose only support is "it's a tRF and it aligns to a host 3′UTR" should be scored near zero.** The discriminating evidence is selective export, not sequence.

---

## 7. The skeptical case

This section deserves disproportionate weight for a tool whose job is to avoid false positives.

### 7.1 Direct failed replications *inside* cross-kingdom RNAi (not just the dietary controversy)

- **Qin et al. 2022 (Molecular Plant Pathology)** sequenced *B. cinerea* and tomato sRNAs and mRNAs in early interaction, predicted 7,042 fungal sRNAs targeting 3,185 tomato mRNAs, found only 163 targets actually down-regulated, generated mutants deleting a transposon region producing ~10% of fungal sRNAs, and generated *Bcdcl1 Bcdcl2* double mutants with **>99% reduction of transposon-derived sRNA production** — and **neither mutant was significantly reduced in virulence on any plant species tested**. Title and conclusion: no functional evidence for naturally occurring ckRNAi in early *B. cinerea*–tomato interaction. **[C — a direct, well-powered negative against the field's founding system.]**
- **Qin et al. 2023 (Frontiers in Plant Science)**, reverse direction: made a fungal mutant with a **5-nt synonymous substitution** abolishing the predicted tomato-sRNA/Bcspl1 interaction without changing protein sequence. The mutant transcript showed the *same* transient down-regulation as wild type, and virulence was identical. **[C — the cleanest kind of negative control available, and it came out negative.]**
- **He et al. 2023 (bioRxiv)** rebut Qin 2022 with reanalysis plus new tomato AGO1-IP showing Bc-sRNA association and target silencing. **[R]** The dispute is unresolved in the published record.
- **Kettles et al. 2019 (Frontiers in Plant Science)**: *Zymoseptoria tritici* — 389 in-planta-induced sRNA loci, predicted wheat targets, but molecular validation of targeting **failed**; DCL and AGO deletion mutants retained full virulence; HIGS was ineffective; and the fungus is **incapable of dsRNA uptake**. **[C]**
- **Ma et al. 2020 (Frontiers in Plant Science)**, same pathosystem, independent group: degradome analysis found **no clear evidence of cross-kingdom RNAi acting by mRNA cleavage**, and fungal sRNA enrichment was *lower* in planta than in vitro. **[C, replicating Kettles.]**
- **Porquier et al. 2026 (bioRxiv)**: *B. cinerea* on grapevine — retrotransposon-derived sRNAs had predicted host targets, but **predicted cleavage sites could not be experimentally validated**, and Dicer mutants unable to produce them **remained fully pathogenic**. **[C]**
- **Ruf et al. 2024 (Plant Direct)**: across a panel of microbes, several interactions showed wild-type-like infection in *all* DCL/AGO mutants. **[C, scope-limiting.]**

So: the fungal ckRNAi literature contains at least **five independent negative or failed-validation results**, two of them in the flagship *Botrytis* system. Any Bayesian prior for a predicted cross-kingdom pair must be set with these in view.

### 7.2 The dietary plant-miRNA controversy and its technical lessons

- Zhang et al. 2011 (Cell Research) reported MIR168a from rice in human/mouse serum repressing LDLRAP1. **[R, heavily contested.]**
- Zhang et al. 2012 (BMC Genomics): plant miRNAs present in animal datasets, **miR168 extremely over-represented**, >96% of miR168 sequences monocot-derived even in insects reared on dicots; feeding studies showed no accumulation. Conclusion: **the signal can originate in the sequencing process itself**. **[D]**
- Tosar et al. 2014 (RNA): designed an experiment showing **cross-contamination during library preparation** generates exogenous reads, and that these contamination-derived sequences **resist periodate oxidation** (i.e. the standard plant-miRNA authentication control fails). They further showed the plant-miRNA spectrum in the originating group's human sera correlated strongly with that group's earlier amphioxus libraries. **[D]** This is the most damaging single methodological result in the field.
- Lusk 2014 (PLoS ONE): contamination explains all observations behind a high-profile "genes pass from food to human blood" claim; contamination can **replicate across independent experiments and vary across conditions**, so replication is not proof. **[D]**
- Kang et al. 2017 (RNA): meta-study of **824 datasets** — xenomiRs present in 17% of tissues and 69% of body fluids but at **0.001% of host miRNA counts**; **no enrichment** in diet-exposed tissues such as liver; **no depletion** in tissues separated from the bloodstream such as brain and CSF; **81% of body-fluid xenomiRs come from rodents** (rare in the human diet, common in labs); samples cluster by study, indicating batch effects; and controlled feeding experiments detected no transfer. **[D — effectively definitive for the dietary claim.]**
- Witwer 2018 (F1000Research): the main "plant miRNA" in human plasma, peu-MIR2910, is 100% identical to a human rRNA. **[D]**
- Huang et al. 2017/2018 (PLoS ONE / Nutrients): with properly validated periodate controls and spike-ins, false positives in liver disappear and **no corn miRNA is detected** in mouse plasma or liver; <1% of dietary miRNA survives simulated digestion. **[D]**
- Heintz-Buschart et al. 2018 (BMC Biology): **sRNA contaminants are present in widely used miRNA extraction columns**, traced by qPCR, and the same artefactual sequences appear in published datasets; they explicitly call for re-evaluation of reports of exogenous microbial and dietary RNA in plasma. **[D]**
- Beatty et al. 2014 (BMC Genomics) illustrates the trap positively: abundant non-human sRNAs in plasma, mostly rRNA-derived, with >90% of fungal reads from the order **Hypocreales** — a classic reagent/environment signature presented as a microbiome.

On the other side, Chen et al. 2021 (PLoS ONE, 421 datasets), Zhao et al. 2020 (PLoS ONE, 388 datasets), Luo et al. 2017 (Sci Rep, pigs) and Olmi et al. 2023 report positive evidence; the most recent reviews (Matsagani et al. 2026, IJMS; Tambaro et al. 2025, Advances in Nutrition) conclude the hypothesis remains unproven, with Tambaro stating plainly that **no direct molecular evidence shows plant xenomiRs binding mammalian silencing machinery**. **[C]**

### 7.3 Index hopping and cross-talk — the quantitative numbers

- Index swapping on patterned flow cells (HiSeqX, HiSeq 3000/4000, NovaSeq) occurs at **0.2–6%** of reads; non-redundant dual indexing eliminates it (Costello et al. 2018, BMC Genomics). **[D]**
- Combinatorial adapters give cross-talk up to **0.29%**, i.e. ~110,000–1,121,000 misassigned reads per lane; unique dual-matched indices reduce this to **≤1 misassigned read per lane** (MacConaill et al. 2018, BMC Genomics). **[D]**
- Model-based estimates of per-read hop probability of **0.003–0.009** nonetheless yield >8% phantom molecules in >25% of samples, and up to **85% in low-complexity samples** (Farouni et al. 2020, Nature Communications). **[D]** — *a Wolbachia sRNA library is exactly a low-complexity sample.*
- **miRTrace** found evidence that **>7% of 700+ public miRNA datasets are cross-contaminated** (Kang et al. 2018, Genome Biology). **[D]**
- Typical RNA-seq carries **10³–10⁵ contaminating microbial reads per million host reads** (Park et al. 2019, BMC Biology). **[D]**
- Database/classifier artefacts: in a controlled plasma study, **per-base-shuffled reads** that destroy all biological sequence were still abundantly classified, still clustered by specimen type, and still correlated (ρ≈0.7) with a published cancer-microbiome catalogue; genus counts scaled with k-mer representation in the database (r²=0.74 real, 0.85 shuffled) (Brumit et al. 2026, bioRxiv). **[R, but methodologically important]** The recommended controls — specimen-matched negative controls, a second conservative classifier, per-base shuffling as an artefact floor, GC-aware analysis — are directly transferable.
- Salzberg et al. 2026 (Nature Cancer) set out required QC and validation for any claim of microbial sequence in an atypical site, including the risk that **contaminated reference genomes cause misclassification of host reads**.

### 7.4 Synthesis of the skeptical case

The field's positive results cluster in a few well-resourced labs and a few pathosystems; the negative results are independent, span multiple systems, and in two cases (Qin 2022, Qin 2023) use the strongest available controls. Simultaneously, the technical literature establishes that (i) sRNA libraries routinely contain exogenous sequence at the 0.001–1% level from reagents, index hopping and library cross-talk; (ii) those artefacts replicate across experiments and resist the field's standard authentication control; and (iii) sequence complementarity between two transcriptomes is essentially uninformative at the scale of modern sRNA catalogues. A tool that reports a ranked list of predicted cross-kingdom targets without addressing (i)–(iii) will produce a list dominated by artefacts.

---

## 8. Specifically for endosymbionts: what changes, and what the literature does not say

### 8.1 What is actually known about *Wolbachia* sRNA

- **Mayoral et al. 2014 (PNAS)** identified two abundant *Wolbachia* snRNAs, **WsnRNA-46 and WsnRNA-49**, expressed from non-coding transcripts whose precursors contain **stem-loop structures**, detected in *Ae. aegypti* (wMelPop-CLA), *D. melanogaster* (wMelPop) and *D. simulans* (wAu) — i.e. conserved across host species and *Wolbachia* strains. The paper states they "may potentially regulate host genes and *Wolbachia* genes". **[R — and note the hedged language: the targeting is predicted, not demonstrated.]** The processed length of the mature WsnRNAs is not stated in material I could verify here; the paper's framing is of prokaryotic snRNAs of 50–500 nt being processed from hairpins.
- **Woolfit et al. 2015 (PLoS ONE)** identified candidate *Wolbachia* sRNAs by two independent routes (comparative genomics and RNA-seq), experimentally characterised one across four strains, and showed it is **differentially regulated by host tissue and sex**. They hypothesise a role in host interaction. **[R]** No host target is proposed.
- **Zhang et al. 2025 (iScience)** identified **12 WsnRNAs** in *Wolbachia*-infected *Tetranychus truncatus* with predicted precursors and hairpin structures; predicted targets of five highly expressed WsnRNAs are enriched for reproductive development; **inhibition of WsnRNA-744 and WsnRNA-3640 reduced fecundity, while inhibition of WsnRNA-6108 increased it**. **[R]** The functional assay is inhibition-of-WsnRNA → host phenotype, which is a real result, but the intervening step — WsnRNA binding a specific host transcript in a host AGO — is **predicted only**. Opposing effects of different WsnRNAs on the same phenotype also mean the system is not a simple one-way silencing channel.

That is the entire direct *Wolbachia*-sRNA-targets-host literature: **three papers, no AGO-IP, no degradome, no 5′ RACE, no validated pair.**

By contrast, the **host-miRNA-side** *Wolbachia* literature is large and mechanistically solid: aae-miR-2940 is induced by *Wolbachia* and its manipulation of a metalloprotease gene is required for endosymbiont maintenance (Hussain et al. 2011, PNAS) **[D]**; aae-miR-12 regulates MCT1 and MCM6 and affects *Wolbachia* density (Osei-Amo et al. 2012, PLoS ONE) **[R]**; nov-miR-12 → *psq* 3′UTR links *Wolbachia* to CI phenotypes (Zheng et al. 2019) **[R]**; dme-miR-92b → *crebA* 3′UTR links *Wolbachia* to learning and memory (Bi et al. 2018) **[R]**; aae-miR-34-3p → MyD88/Toll mediates dengue blocking (She et al. 2024, PLoS Pathogens) **[D]**; Wolbachia-driven host miRNAs affect fecundity and symbiont density in *T. truncatus* (Zhang et al. 2026, Communications Biology) **[R]**. And one strong negative: in a *Drosophila* cell/SFV model, *Wolbachia* antiviral protection occurred **without any transcriptional response by the bacterium and without changes in host miRNA expression or involvement of the antiviral siRNA pathway** (Rainey et al. 2016, PLoS Pathogens). **[C]**

**The field's own centre of gravity is therefore "Wolbachia changes host miRNAs", not "Wolbachia sends sRNAs that act as miRNAs".** A tool predicting the latter is working in a space where the former is the better-supported alternative explanation for almost any observed host phenotype — and the two are experimentally confounded unless the candidate sRNA is shown to be of bacterial sequence origin *and* in a host AGO.

### 8.2 What genuinely changes for an intracellular, vertically transmitted, reduced-genome symbiont

**Plausibly favourable (all [S] unless noted):**
1. **No extracellular barrier.** Secretion into the host cytoplasm does not require OMV budding, survival in an extracellular milieu, and endocytic uptake. The RNase-protection argument that motivates vesicle packaging is weaker. *Wolbachia* does sit inside a host-derived vacuolar membrane, so one membrane still has to be crossed — unaddressed in the literature. **[gap]**
2. **Long coevolution and strict vertical transmission** align symbiont and host fitness and make a slow, quantitative regulatory channel evolutionarily plausible in a way that is harder to argue for a transient pathogen. The *Buchnera* result (host miRNAs + host Ago1 translocated into the symbiont; all translocated miRNAs predicted to target symbiosis-function genes; Ding et al. 2026) establishes that *something* crosses in obligate endosymbiosis. **[R]**
3. **Tissue specificity and dosage.** *Wolbachia* sRNA candidates are differentially regulated by host tissue and sex (Woolfit 2015) **[R]**, which is the kind of conditional expression a regulatory channel would need, and which also gives a usable experimental contrast.
4. **Pre-existing host small-RNA machinery is demonstrably perturbed** by *Wolbachia*: an overall increase in 18–28 nt small RNAs in both nucleus and cytoplasm of infected cells, altered miRNA nucleo-cytoplasmic shuttling, and altered piRNA abundance (Mayoral et al. 2014, PLoS ONE). **[R]** A symbiont that is already modulating RISC is better placed to exploit it.

**Plausibly unfavourable:**
5. **Reduced genome, reduced regulatory apparatus.** *Wolbachia* is a genome-reduced obligate intracellular bacterium. Whether it retains Hfq (the chaperone that mediates OMV RNA loading in free-living Gram-negatives; Mosca 2026, El Hamidi/Hfq literature) is not addressed in the sRNA papers I retrieved. **[gap]** Without Hfq-type chaperoning, selective RNA export is harder to rationalise.
6. **No DCL, no RdRP, no AGO.** Bacteria have no RNAi machinery, so there is no biogenesis pathway that produces a defined 21–22 nt, 5′-U, 2-nt-3′-overhang duplex. Any bacterial guide reaching host AGO must be produced by **nuclease processing of a hairpin or a tRNA/rRNA**, with essentially no control over 5′ end identity — the exact feature that governs plant AGO sorting and matters for animal AGO loading too. This is, in my judgement, the single strongest a-priori argument against bacterium→host canonical RNAi, and it is **not** confronted anywhere in the bacterial sRNA-effector literature. **[S, mine]**
7. **No published AGO-IP from any symbiont-infected animal tissue.** The entire argument for host-AGO loading of a bacterial sRNA rests on prediction. **[gap]**

**Where we would be extrapolating, explicitly:**
- From plant AGO1 5′U/21-nt loading rules → to insect/arachnid AGO1/AGO2. Not tested for any bacterial sRNA.
- From OMV-mediated delivery in free-living Gram-negatives → to cytoplasmic release by an intracellular symbiont. Not tested.
- From helminth EV-miRNA mimicry (where the guides are *bona fide* Dicer products of a eukaryote) → to bacterial hairpin/tRNA fragments. The biogenesis is not comparable.
- From *Wolbachia*-responsive host miRNA phenotypes → to *Wolbachia*-encoded sRNA function. These are routinely conflated and are not the same claim.

---

## 9. Implications for a target predictor

### 9.1 Pairing rules to implement

1. **Do not implement a single rule. Implement two modes, chosen by host kingdom.**
   - *Plant-host mode* (relevant if SmithHunter2 is ever applied to plant-associated bacteria, e.g. rhizobia, *Xanthomonas*): plant-style penalty scoring over the full guide (psRNATarget/TargetFinder-class), allowing ≤4–5 mismatch-equivalents, penalising mismatches at positions 2–13 more heavily, and **requiring pairing across the predicted cleavage site (positions 10–11)**. Justification: every validated plant-host cross-kingdom pair is cleavage-type (Weiberg 2013; Dunker 2020; Silvestri 2024; Ren 2019).
   - *Animal-host mode* (the Wolbachia case): require a **contiguous Watson–Crick 7-mer in guide positions 2–8**, optionally extended (7mer-m8, 8mer with a t1A), because every validated animal-host cross-kingdom pair for which a mode is known is seed-driven miRNA mimicry (Tran 2021; Pérez 2025; Zhong 2024). Treat 6-mers as unscored.
2. **Do not advertise free-energy thresholds as evidence-based** — the literature reports none for cross-kingdom pairs. Use ΔG (RNAhybrid/RNAduplex) only as a tiebreaker within an already-filtered set, and report it, never threshold on it alone.
3. **Allow central bulges in plant mode; forbid seed bulges in animal mode.** No cross-kingdom study quantifies bulge tolerance, so this inherits from endogenous-miRNA biology and should be labelled as such in the tool's documentation.
4. **Score a bidirectional possibility.** Ile-tRF-5X *increased* MAP3K4 (Diallo 2022) and VdrsR-1 *increased* miR157d by targeting a precursor (Zhang 2022). A silencing-only output model is known to be incomplete.
5. **Normalise against a shuffled/dinucleotide-preserving null.** Given the "97% of fungal genes targetable" result (Qin 2023), any hit must be reported with an empirical p-value against sequence-composition-matched shuffles of the same sRNA, not just a raw score. Per-base shuffling as an artefact floor is now an explicit recommendation in the low-biomass literature (Brumit 2026).

### 9.2 Length and 5′-nucleotide filters that are defensible

| Filter | Defensible setting | Justification | Risk |
|---|---|---|---|
| Length, primary window | **19–24 nt** | AGO-loadable range in both plants and animals; 63/73 validated-class fungal effectors were 20–22 nt (Weiberg 2013) | Excludes AMF 24–26 nt class and any 8–15 nt vsRNA activity |
| Length, reported-but-downweighted | 25–30 nt | AMF sRNAs centre at 24–26 nt and Rir2216 still works (Silvestri 2024); *Colletotrichum* 29-nt class | — |
| Length, exclude from RNAi scoring | **<16 nt and >30 nt** | 8–15 nt vsRNAs dominate bacterial libraries (Diallo 2022) but cannot support canonical RISC; classical bacterial sRNAs are 50–500 nt (Mayoral 2014) | May discard a real non-canonical mechanism; report separately |
| 5′ nucleotide | **Report it; bonus for 5′U, do not require it** | 5′U is what makes the plant-AGO1 story work; for animal AGO the preference is weaker; and **no bacterial study has reported a 5′-nt distribution at all** | Requiring 5′U would be unjustified extrapolation for an animal host |
| 5′-end definition | Require a **well-defined, reproducible 5′ end** across replicates | Bacteria lack Dicer, so a sharply defined 5′ end is itself evidence of regulated processing rather than degradation | — |
| Hairpin precursor | Bonus, not requirement | Mayoral 2014, Woolfit 2015, Zhang 2025 and Coskun 2020 all report stem-loop precursors; Choi 2017 defines msRNAs by predicted hairpins | Hairpin prediction is cheap and over-permissive on AT-rich genomes |

### 9.3 Regions to search

- **Search the entire mature transcript**, and report the region. Do not restrict to 3′UTR.
- **Normalise for 3′UTR length.** Predicted targeting multiplicity scales with 3′UTR length at r = 0.63 (Zoziuk 2026); without normalisation the output will be a list of long-UTR genes.
- Include **miRNA precursors and other non-coding RNAs** as candidate targets — the VdrsR-1/MIR157d precedent (Zhang 2022) shows mRNA-only target spaces miss real biology.
- For plant hosts, weight CDS sites equally with UTR sites.

### 9.4 Realistic prior probability that a predicted pair is real

Assembling the available numbers:
- In the best-characterised system, **73/832 = 8.8%** of exported-and-enriched fungal sRNAs passed stringent prediction (Weiberg 2013), of which **3** were functionally validated at the time → roughly **0.4%** of the exported pool has a demonstrated host target, and ~4% of the stringently-predicted set.
- Where a transcriptome-wide correlation test was run, **163/3,185 = 5.1%** of predicted targets were actually down-regulated in planta (Qin 2022) — and that 5.1% did not survive genetic testing.
- A HITS-CLIP-anchored prediction within a single organism achieved a **50% verification rate** (Zhao 2015) — i.e. *even with direct AGO-binding evidence in hand*, half the predicted sites were wrong.
- Unanchored complementarity prediction is uninformative: 97% of fungal genes become "targetable" (Qin 2023).
- For the specific case of a bacterial endosymbiont sRNA in an animal host, there is **not one validated pair in the literature**.

**Defensible prior, bacteria→animal host, sequence prediction alone: ≲1%, and arguably <0.1%.** With selective-export evidence (enrichment in symbiont-derived vesicles or in a host-cytoplasm fraction over whole-symbiont RNA), a defined 5′ end, a reproducible hairpin precursor, and an expression-anticorrelation with the candidate target, a prior of a **few percent** is arguable. Only AGO-IP or CLIP evidence moves it into the tens of percent, and by the *S. japonicum* benchmark even then only to ~50%. The tool should state this explicitly in its output, and should be designed to **rank and triage for experimental follow-up**, not to assert targets.

### 9.5 Hard filters I would implement before any pairing score is computed

1. **Export/enrichment filter.** Candidate must be enriched in the extracellular/vesicular or host-cytoplasmic fraction relative to whole-symbiont RNA (baseline expectation: OMV sncRNA is ~1/10 cytosolic abundance and selectively composed — Naaz 2026; Diallo 2022; Chapagain 2024).
2. **Contamination filter.** (a) Reject sequences with exact or near-exact identity to host rRNA, tRNA or other abundant host ncRNA — the peu-MIR2910 lesson (Witwer 2018). (b) Require unique dual-index libraries or explicitly model index hopping at 0.003–0.009 per read (Farouni 2020; Costello 2018; MacConaill 2018). (c) Require specimen-matched negative/reagent controls; sRNA contaminants are present in standard extraction columns (Heintz-Buschart 2018). (d) Flag candidates whose abundance correlates with batch rather than with biology (Kang 2017). (e) Run a miRTrace-style taxonomic provenance check (Kang 2018).
3. **Non-RNAi alternative-hypothesis flag.** Any candidate that is a tRF, rRF, tmRNA fragment or Y-RNA fragment should carry a flag that its host effect may be **TLR7/TLR3-mediated innate sensing** rather than base-pairing (Gumas 2024; Rodriguez 2020; Sheikh 2025; Ha 2020).
4. **Confounder flag.** For symbiont systems, flag that the observed host phenotype may be mediated by **symbiont-induced changes in host miRNAs** rather than by symbiont-encoded sRNAs — the better-supported mechanism in *Wolbachia* (Hussain 2011; She 2024) and the explicit negative control that is almost always missing.

### 9.6 What the field does not know (and SmithHunter2 should not pretend to)

1. **No measured pairing statistics for any set of validated cross-kingdom pairs.** No published mismatch distribution, bulge-tolerance measurement, ΔG distribution or seed-class distribution.
2. **No 5′-nucleotide distribution reported for any bacterial sRNA claimed to target a host mRNA.**
3. **No AGO-IP or CLIP evidence for any bacterial sRNA in any host Argonaute.** Not one.
4. **No degradome or 5′-RACE evidence for any bacterial→animal pair.**
5. **No test of 3′UTR vs CDS vs 5′UTR enrichment for cross-kingdom targets on equal footing.** The 3′UTR consensus is a methodological artefact of the tools used.
6. **No mechanism for RNA crossing the symbiosome/vacuolar membrane** in an intracellular endosymbiont.
7. **No resolution of whether *Wolbachia* retains the Hfq-type machinery** that mediates selective RNA export in free-living Gram-negatives.
8. **No validated *Wolbachia* sRNA–host transcript pair.** The three existing papers report prediction plus, in one case, a host phenotype after sRNA inhibition.
9. **Unresolved in the flagship system:** whether *B. cinerea* ckRNAi is required for virulence at all (Weiberg 2013/Wang 2016/He 2023 vs Qin 2022/2023/Porquier 2026).
10. **Unknown whether stoichiometry is ever sufficient.** The one honest statement on effect size in the parasite literature is that repression is "subtle" and must rely on synergy across dozens of species (Soichot 2022). Nobody has measured how many copies of a symbiont sRNA per host cell would be needed to occupy enough RISC to matter.

---

## References

1. Weiberg A. et al. (2013). *Fungal Small RNAs Suppress Plant Immunity by Hijacking Host RNA Interference Pathways.* Science. https://consensus.app/papers/details/w:999aaa2e-ddf2-4d63-8f8f-8981293722b3/
2. Weiberg A. et al. (2025). *Fungal Argonaute proteins act in bidirectional cross-kingdom RNA interference during plant infection.* PNAS. https://consensus.app/papers/details/w:a4162bfa-3fe1-4d7d-aee6-0ad87909add6/
3. He B. et al. (2023). *Fungal small RNAs ride in extracellular vesicles to enter plant cells through clathrin-mediated endocytosis.* Nature Communications. https://consensus.app/papers/details/w:1c72a505-fcbe-4c55-b5a4-92165443d150/
4. Wang M. et al. (2016). *Bidirectional cross-kingdom RNAi and fungal uptake of external RNAs confer plant protection.* Nature Plants. https://consensus.app/papers/details/w:cab5fcb3-578e-4f12-b617-4db7a52fec13/
5. He B. et al. (2021). *RNA-binding proteins contribute to small RNA loading in plant extracellular vesicles.* Nature Plants. https://consensus.app/papers/details/w:7c95139b-dedc-47ba-80a1-06c71af2b497/
6. Wang M. et al. (2017). *Botrytis small RNA Bc-siR37 suppresses plant defense genes by cross-kingdom RNAi.* RNA Biology. https://consensus.app/papers/details/w:51becb29-69f3-4335-847b-2ee7cd6ed8e2/
7. Cai Q. et al. (2018). *Plants send small RNAs in extracellular vesicles to fungal pathogen to silence virulence genes.* Science. https://consensus.app/papers/details/w:69cc5ea5-2f44-4411-8fd8-a36dc75ddfc7/
8. Cheng A.-P. et al. (2023). *A fungal RNA-dependent RNA polymerase is a novel player in plant infection and cross-kingdom RNA interference.* PLoS Pathogens. https://consensus.app/papers/details/w:28c0baf6-d6df-4dde-8827-02fe3577c648/
9. Qin S. et al. (2022). *Molecular characterization reveals no functional evidence for naturally occurring cross-kingdom RNA interference in the early stages of Botrytis cinerea–tomato interaction.* Molecular Plant Pathology. https://consensus.app/papers/details/w:acbf9b72-fc3f-49db-91a2-89dc2e841073/
10. Qin S. et al. (2023). *Molecular characterization of cross-kingdom RNA interference in Botrytis cinerea by tomato small RNAs.* Frontiers in Plant Science. https://consensus.app/papers/details/w:18a2e25a-34e5-4301-8eae-df4963cd3f36/
11. He B. et al. (2023). *Botrytis cinerea small RNAs are associated with tomato AGO1 and silence tomato defense-related target genes supporting cross-kingdom RNAi.* bioRxiv. https://consensus.app/papers/details/w:eb4ac8a0-4e36-4090-b62c-2a6b5cbbec97/
12. Dunker F. et al. (2020). *Oomycete small RNAs bind to the plant RNA-induced silencing complex for virulence.* eLife. https://consensus.app/papers/details/w:efc7a66d-2adc-46ef-9414-1dbf7f3d21a6/
13. Dunker F. et al. (2021). *Plant ARGONAUTE Protein Immunopurification for Pathogen Cross Kingdom Small RNA Analysis.* Bio-protocol. https://consensus.app/papers/details/w:e21bebf5-8ff7-48d4-80a0-6e3d190a5eca/
14. Kusch S. et al. (2023). *Site-specific analysis reveals candidate cross-kingdom small RNAs, tRNA and rRNA fragments, and signs of fungal RNA phasing in the barley–powdery mildew interaction.* Molecular Plant Pathology. https://consensus.app/papers/details/w:23dd2295-e084-474e-8378-e9cf734b9dbf/
15. Ni J. et al. (2026). *A rice ceRNA module suppresses Rhizoctonia solani–induced cross-kingdom RNAi to reduce fungal pathogenicity.* Nature Communications. https://consensus.app/papers/details/w:855551ab-724e-4a97-9f58-78f4c2be84f0/
16. Silvestri A. et al. (2024). *A fungal sRNA silences a host plant transcription factor to promote arbuscular mycorrhizal symbiosis.* New Phytologist. https://consensus.app/papers/details/w:24295e90-f619-4374-b6b6-dd0140bba739/
17. Silvestri A. et al. (2019). *In silico analysis of fungal small RNA accumulation reveals putative plant mRNA targets in the symbiosis between an arbuscular mycorrhizal fungus and its host plant.* BMC Genomics. https://consensus.app/papers/details/w:1032ac44-0488-4ccb-ba80-607c65ef46d0/
18. Usländer A. et al. (2026). *Cross-kingdom RNA interference promotes arbuscular mycorrhiza development.* Nature Plants. https://consensus.app/papers/details/w:86c99291-4e80-432a-938c-e91a8ee41227/
19. Kellari L. M. et al. (2025). *Cross-kingdom RNAi: a universal mechanism of inter-organismal communication with many unknowns.* Journal of Experimental Botany. https://consensus.app/papers/details/w:ea49a3f5-11f1-421e-bd86-b12852d9e574/
20. Kellari L. M. et al. (2025). *Cross-kingdom RNAi induced by a beneficial endophytic fungus to its host requires transitivity and amplification of silencing signals.* Plant Biology. https://consensus.app/papers/details/w:71cc9703-fff8-4fca-a09a-9e1b8e82d5b5/
21. Nasfi S. et al. (2024). *A pipeline for validation of Serendipita indica effector-like sRNA suggests cross-kingdom communication in the symbiosis with Arabidopsis.* Journal of Experimental Botany. https://consensus.app/papers/details/w:b599665d-2a43-45ad-8de3-9ff7d5486a42/
22. Ruf A. et al. (2024). *Broad-scale phenotyping in Arabidopsis reveals varied involvement of RNA interference across diverse plant-microbe interactions.* Plant Direct. https://consensus.app/papers/details/w:f5479280-c5bc-48e1-b6fe-133e2ae6034b/
23. Werner B. T. et al. (2025). *Cross-kingdom degradomics identifies natural plant small RNAs with disease protecting activity.* bioRxiv. https://consensus.app/papers/details/w:0d13771f-4dfe-49d7-8425-dacd506a7456/
24. Werner B. T. et al. (2021). *Fusarium graminearum DICER-like-dependent sRNAs are required for the suppression of host immune genes and full virulence.* PLoS ONE. https://consensus.app/papers/details/w:be72c50c-23fb-4e44-bced-0769edb02707/
25. Zhang B. et al. (2022). *Verticillium dahliae Secretes Small RNA to Target Host MIR157d and Retard Plant Floral Transition During Infection.* Frontiers in Plant Science. https://consensus.app/papers/details/w:2593e663-9a60-4dcb-b2fd-95d9a8140e0d/
26. Zhu C. et al. (2022). *A fungal effector suppresses the nuclear export of AGO1–miRNA complex to promote infection in plants.* PNAS. https://consensus.app/papers/details/w:2b538812-0f00-4139-9eba-cb9335535190/
27. Ji H.-M. et al. (2021). *Fol-milR1, a pathogenicity factor of Fusarium oxysporum, confers tomato wilt disease resistance by impairing host immune responses.* New Phytologist. https://consensus.app/papers/details/w:240712ac-61c9-484f-a75c-24b83881659b/
28. Porquier A. et al. (2021). *Retrotransposons as pathogenicity factors of the plant pathogenic fungus Botrytis cinerea.* Genome Biology. https://consensus.app/papers/details/w:40922767-369d-4c05-b74b-4f0076801b9d/
29. Porquier A. et al. (2026). *Role of transposons in the specialization of Botrytis cinerea to grapevine: Insights into small RNAs and a new Starship.* bioRxiv. https://consensus.app/papers/details/w:4b42a3cd-9c5f-4615-8df4-7960b5b0afe8/
30. Kettles G. J. et al. (2019). *sRNA Profiling Combined With Gene Function Analysis Reveals a Lack of Evidence for Cross-Kingdom RNAi in the Wheat–Zymoseptoria tritici Pathosystem.* Frontiers in Plant Science. https://consensus.app/papers/details/w:e1b0959a-f0d2-4127-9276-b78b82b204eb/
31. Ma X. et al. (2020). *Small RNA Bidirectional Crosstalk During the Interaction Between Wheat and Zymoseptoria tritici.* Frontiers in Plant Science. https://consensus.app/papers/details/w:1f4e755c-28fa-4490-97ea-7c8b7e8a0249/
32. Derbyshire M. C. et al. (2019). *Small RNAs from the plant pathogenic fungus Sclerotinia sclerotiorum highlight host candidate genes associated with quantitative disease resistance.* Molecular Plant Pathology. https://consensus.app/papers/details/w:4aae2f1f-9dc5-473b-8f85-612e3d03a177/
33. Armijos C. E. et al. (2026). *Stage-Specific RNA Turnover Drives Small RNA Dynamics in Arabidopsis–Colletotrichum Interactions.* Plant Direct. https://consensus.app/papers/details/w:709fa091-cdfd-4bfb-9616-59a39d0855b1/
34. Qiao L. et al. (2021). *Spray-induced gene silencing for disease control is dependent on the efficiency of pathogen RNA uptake.* Plant Biotechnology Journal. https://consensus.app/papers/details/w:f086cb46-8445-4fd8-a11a-71624773cf84/
35. Zhang T. et al. (2022). *Exploring the Effectiveness and Durability of Trans-Kingdom Silencing of Fungal Genes in the Vascular Pathogen Verticillium dahliae.* IJMS. https://consensus.app/papers/details/w:e22265ef-da8c-4acd-b592-3894540721f6/
36. Wang S. et al. (2023). *Plant mRNAs move into a fungal pathogen via extracellular vesicles to reduce infection.* Cell Host & Microbe. https://consensus.app/papers/details/w:a2498be8-4f9e-45b4-b3a6-a82d71b5541c/
37. Wang M. et al. (2026). *A specialized ARGONAUTE enables trans-species RNA interference in plant immunity.* bioRxiv. https://consensus.app/papers/details/w:1d1098b8-dd28-4c73-a832-2eb91ed36c42/
38. Cui C. et al. (2019). *A fungal pathogen deploys a small silencing RNA that attenuates mosquito immunity and facilitates infection.* Nature Communications. https://consensus.app/papers/details/w:cddb578d-6ded-4821-9246-40d37d212468/
39. Ravet A. et al. (2025). *Vesicular and non-vesicular extracellular small RNAs direct gene silencing in a plant-interacting bacterium.* Nature Communications. https://consensus.app/papers/details/w:726cbc2f-b3ca-44ac-8339-fe118b04a0ab/
40. Chen X. & Rechavi O. (2021). *Plant and animal small RNA communications between cells and organisms.* Nature Reviews Molecular Cell Biology. https://consensus.app/papers/details/w:cab5edba-b76a-4fce-9ad0-5ad9d7d49d0b/
41. Ledford W. C. et al. (2023). *A journey into the world of small RNAs in the arbuscular mycorrhizal symbiosis.* New Phytologist. https://iris.unito.it/retrieve/handle/2318/1945777/1214325/
42. Koeppen K. et al. (2016). *A Novel Mechanism of Host-Pathogen Interaction through sRNA in Bacterial Outer Membrane Vesicles.* PLoS Pathogens. https://consensus.app/papers/details/w:cb737203-e34d-4dc4-a88a-b44eeeeea042/
43. Sahr T. et al. (2022). *Translocated Legionella pneumophila small RNAs mimic eukaryotic microRNAs targeting the host immune response.* Nature Communications. https://consensus.app/papers/details/w:89560f15-8b21-4861-9ff7-bbee6846d2c2/
44. Choi J.-W. et al. (2017). *Secretable Small RNAs via Outer Membrane Vesicles in Periodontal Pathogens.* Journal of Dental Research. https://consensus.app/papers/details/w:5b9bb4cf-6e6a-4ba8-857f-3a4083a12faf/
45. Diallo I. et al. (2022). *A tRNA-derived fragment present in E. coli OMVs regulates host cell gene expression and proliferation.* PLoS Pathogens. https://consensus.app/papers/details/w:8c70e42f-6a88-4b73-bb0b-9a9e01b4a752/
46. Diallo I. et al. (2022). *RNA Sequencing Unveils Very Small RNAs With Potential Regulatory Functions in Bacteria.* Frontiers in Molecular Biosciences. https://consensus.app/papers/details/w:c23d46c9-adbd-4299-b906-368a91c11acd/
47. Li Z. & Stanton B. A. (2021). *Transfer RNA-Derived Fragments, the Underappreciated Regulatory Small RNAs in Microbial Pathogenesis.* Frontiers in Microbiology. https://consensus.app/papers/details/w:0903e0f9-7e6b-488c-9616-d7440aa74729/
48. Xie Z. et al. (2024). *Pseudomonas aeruginosa outer membrane vesicle-packed sRNAs can enter host cells and regulate innate immune responses.* Microbial Pathogenesis. https://consensus.app/papers/details/w:9265f445-3e72-4f84-835f-02cf4d2b9dfa/
49. Wu Y. et al. (2024). *Suppression of host plant defense by bacterial small RNAs packaged in outer membrane vesicles.* Plant Communications. https://consensus.app/papers/details/w:b0a24d4f-09c8-479d-a7df-420d6223ea05/
50. Yu S. et al. (2022). *Biological Functions and Cross-Kingdom Host Gene Regulation of Small RNAs in Lactobacillus plantarum-Derived Extracellular Vesicles.* Frontiers in Microbiology. https://consensus.app/papers/details/w:965207fd-73ef-4f84-b328-d21c7c73f481/
51. Zhang H. et al. (2019). *sncRNAs packaged by Helicobacter pylori outer membrane vesicles attenuate IL-8 secretion in human cells.* Int. J. Medical Microbiology. https://consensus.app/papers/details/w:8c90d146-60c4-47fb-aea7-d4b1c1419879/
52. Moriano-Gutierrez S. et al. (2020). *The noncoding small RNA SsrA is released by Vibrio fischeri and modulates critical host responses.* PLoS Biology. https://consensus.app/papers/details/w:34246a46-5582-4868-809a-67d9eb7ae0b9/
53. Montes-García J. F. et al. (2026). *Bacterial small RNAs and extracellular vesicles in inflammatory and osteoimmune diseases.* Archives of Microbiology. https://consensus.app/papers/details/w:cb64478f-141d-47e5-aaf5-42239a1e3408/
54. (P. gingivalis OMV msRNA sRNA45033 / CBX5 primary study). *Porphyromonas gingivalis Outer Membrane Vesicles Promote Apoptosis via msRNA-Regulated DNA Methylation in Periodontitis.* https://consensus.app/papers/details/w:daa7f657-3651-43c2-b305-82cbae5b9be5/
55. Ha J. Y. et al. (2020). *Delivery of Periodontopathogenic Extracellular Vesicles to Brain Monocytes and Microglial IL-6 Promotion by RNA Cargo.* Frontiers in Molecular Biosciences. https://consensus.app/papers/details/w:6bb2446a-65e9-4573-bc9a-df5d08bc9cfa/
56. Chen P. et al. (2026). *A sRNA Subset Present in Brucella melitensis 16M OMVs Regulates Intracellular Survival.* Kafkas Univ. Vet. Fak. Derg. https://consensus.app/papers/details/w:64e5265c-ef6c-4ac0-86dd-11b60c2b502a/
57. Malabirade A. et al. (2018). *The RNA Complement of Outer Membrane Vesicles From Salmonella enterica Serovar Typhimurium Under Distinct Culture Conditions.* Frontiers in Microbiology. https://consensus.app/papers/details/w:c4a57566-aa74-4def-856a-408b944c2b49/
58. Naaz S. et al. (2026). *P22 Small Noncoding RNAs Are Actively Secreted in Salmonella Outer Membrane Vesicles During Bacteriophage Infection.* Non-Coding RNA. https://consensus.app/papers/details/w:02635ec7-fcdd-4b14-bace-a64ee1ff2049/
59. Chapagain P. et al. (2024). *Characterisation of sRNAs enriched in outer membrane vesicles of pathogenic Flavobacterium psychrophilum.* Journal of Extracellular Biology. https://consensus.app/papers/details/w:cb28123e-dc1c-4a88-b8f9-2c2a9bb0eb53/
60. Mosca K. et al. (2026). *RNAs Associated With Bacterial Outer Membrane Vesicles: Structural Insights Into Surface Composition.* Journal of Extracellular Vesicles. https://consensus.app/papers/details/w:a36bd6ff-0663-4011-9b8e-946cc55267f0/
61. Bittel M. et al. (2021). *Visualizing transfer of microbial biomolecules by outer membrane vesicles in microbe-host-communication in vivo.* Journal of Extracellular Vesicles. https://consensus.app/papers/details/w:2520b789-5563-44d8-a49d-ab4fa8d5214b/
62. Dauros-Singorenko P. et al. (2018). *The functional RNA cargo of bacterial membrane vesicles.* FEMS Microbiology Letters. https://consensus.app/papers/details/w:5887b838-8738-4139-b6e3-550e7435665b/
63. Rodriguez B. V. & Kuehn M. J. (2020). *Staphylococcus aureus secretes immunomodulatory RNA and DNA via membrane vesicles.* Scientific Reports. https://consensus.app/papers/details/w:5ab832cd-d500-4a58-afd4-a37a22e54450/
64. Sheikh A. et al. (2025). *Outer Membrane Vesicles From Bacteroides fragilis Contain Coding and Non-Coding Small RNA Species That Modulate Inflammatory Signalling in Intestinal Epithelial Cells.* Journal of Extracellular Biology. https://consensus.app/papers/details/w:423d60d6-978c-44df-9f0f-0b2b01a3408b/
65. Gumas J. et al. (2024). *Immunostimulatory short non-coding RNAs in the circulation of patients with tuberculosis infection.* Molecular Therapy — Nucleic Acids. https://consensus.app/papers/details/w:965dd946-3bbb-4aa4-a6e4-9bdc7ddeec99/
66. Coskun F. S. et al. (2020). *sncRNA-1 Is a Small Noncoding RNA Produced by Mycobacterium tuberculosis in Infected Cells.* Frontiers in Microbiology. https://consensus.app/papers/details/w:41ffee15-9d94-4e82-b359-ac18355480b9/
67. Ren B. et al. (2019). *Rhizobial tRNA-derived small RNAs are signal molecules regulating plant nodulation.* Science. https://www.purdue.edu/newsroom/releases/2019/Q3/purdue-scientists-show-how-soybeans-may-get-more-nitrogen-from-atmosphere.html
68. Buck A. H. et al. (2015). *Exosomes secreted by nematode parasites transfer small RNAs to mammalian cells and modulate innate immunity.* Nature Communications. https://consensus.app/papers/details/w:0ce27e71-f15f-4485-9fb8-9a9c65ddc4da/
69. Chow F. W.-N. et al. (2019). *Secretion of an Argonaute protein by a parasitic nematode and the evolution of its siRNA guides.* Nucleic Acids Research. https://consensus.app/papers/details/w:cdf56b92-9b20-4135-80e8-423a42393438/
70. Soichot J. et al. (2022). *Nematode microRNAs can Individually Regulate Interferon Regulatory Factor 4 and mTOR in Differentiating T Helper 2 Lymphocytes.* Frontiers in Molecular Biosciences. https://consensus.app/papers/details/w:a01a6129-ba1d-431e-be49-115d5f7f4e0d/
71. Pérez M. G. et al. (2025). *A secreted helminth microRNA suppresses gastrointestinal cell differentiation required for innate immunity.* Frontiers in Immunology. https://consensus.app/papers/details/w:56a4de0a-9884-4d6e-8016-846ebe734dd3/
72. Zhong H. et al. (2024). *Roles of helminth extracellular vesicle-derived let-7 in host–parasite crosstalk.* Frontiers in Immunology. https://consensus.app/papers/details/w:4f37b5d6-9e63-4214-a26b-40c8c2e34c31/
73. Duguet T. B. et al. (2020). *Extracellular Vesicle-Contained microRNA of C. elegans as a Tool to Decipher the Molecular Basis of Nematode Parasitism.* Frontiers in Cellular and Infection Microbiology. https://consensus.app/papers/details/w:b49fd9f7-ba06-4f97-ac54-ac82de3902d4/
74. Tran N. et al. (2021). *Fasciola hepatica hijacks host macrophage miRNA machinery to modulate early innate immune responses.* Scientific Reports. https://consensus.app/papers/details/w:f1795b07-13cc-4d03-aad1-9f811b4bb076/
75. Liu J. et al. (2019). *Schistosoma japonicum extracellular vesicle miRNA cargo regulates host macrophage functions facilitating parasitism.* PLoS Pathogens. https://consensus.app/papers/details/w:507afc0a-29ce-422f-85ed-2d813ec5d0cf/
76. Zhao J. et al. (2015). *HITS-CLIP reveals Argonaute-associated microRNAs and targets in Schistosoma japonicum.* Parasites & Vectors. https://consensus.app/papers/details/w:72a83de0-f843-4acc-96b9-95f841c6e5cb/
77. Dussutour A. et al. (2026). *Plant-parasitic nematode microRNAs hijack plant AGO1 to induce host-cell reprogramming.* bioRxiv. https://consensus.app/papers/details/w:7e2e82eb-958e-481c-a3ec-6e4e24a83991/
78. Britton C. et al. (2015). *Application of small RNA technology for improved control of parasitic helminths.* Veterinary Parasitology. https://consensus.app/papers/details/w:276b03aa-2f33-42eb-8f42-8c38c579e3ee/
79. Sotillo J. et al. (2020). *The protein and microRNA cargo of extracellular vesicles from parasitic helminths.* International Journal for Parasitology. https://consensus.app/papers/details/w:858abc84-9d55-4d86-bb8f-f94f3cec9396/
80. Whitehead B. J. et al. (2026). *Selective Packaging of miR-22 and miR-29 Orthologues into Extracellular Vesicles of the Porcine Whipworm Trichuris suis.* Journal of Extracellular Vesicles. https://consensus.app/papers/details/w:688d706d-c550-4183-99eb-147f2cf6551c/
81. Mayoral J. G. et al. (2014). *Wolbachia small noncoding RNAs and their role in cross-kingdom communications.* PNAS. https://consensus.app/papers/details/w:37c60c81-99f1-41da-b4f7-2c1ab5abdc92/
82. Mayoral J. G. et al. (2014). *Wolbachia Infection Modifies the Profile, Shuttling and Structure of MicroRNAs in a Mosquito Cell Line.* PLoS ONE. https://consensus.app/papers/details/w:abaa016c-65e9-4f68-ab71-8d68efa3f49d/
83. Woolfit M. et al. (2015). *Discovery of Putative Small Non-Coding RNAs from the Obligate Intracellular Bacterium Wolbachia pipientis.* PLoS ONE. https://consensus.app/papers/details/w:31030765-68f2-4d4a-b7b2-011c2788532f/
84. Zhang Y. et al. (2025). *Wolbachia-derived small non-coding RNAs exhibit cross-kingdom regulatory effects on host reproduction of a polyphagous mite.* iScience. https://consensus.app/papers/details/w:740f8c2d-4c99-4d82-a568-c4fbc7619538/
85. Zhang Y. M. et al. (2026). *Wolbachia-driven host miRNAs mediate arthropod reproduction in a Wolbachia density-dependent manner.* Communications Biology. https://consensus.app/papers/details/w:8d240e07-a511-48b7-aa79-aeaae40e35d5/
86. Hussain M. et al. (2011). *Wolbachia uses host microRNAs to manipulate host gene expression and facilitate colonization of the dengue vector Aedes aegypti.* PNAS. https://consensus.app/papers/details/w:a64dc3df-3a89-4664-a352-36e7c6a60cee/
87. Osei-Amo S. et al. (2012). *Wolbachia-Induced aae-miR-12 miRNA Negatively Regulates the Expression of MCT1 and MCM6 Genes.* PLoS ONE. https://consensus.app/papers/details/w:d6f96844-34ad-4436-bf83-27da9906aa1f/
88. Zheng Y. et al. (2019). *Small RNA analysis provides new insights into cytoplasmic incompatibility in Drosophila melanogaster induced by Wolbachia.* Journal of Insect Physiology. https://consensus.app/papers/details/w:cf6f859f-9d87-4c4f-9a68-6292e6120dc5/
89. Bi J. et al. (2018). *Wolbachia infection may improve learning and memory capacity of Drosophila by altering host gene expression through microRNA.* Insect Biochemistry and Molecular Biology. https://consensus.app/papers/details/w:efcb2fd0-8821-48a0-9dff-f30b90e0140f/
90. She L. et al. (2024). *Wolbachia mediates crosstalk between miRNA and Toll pathways to enhance resistance to dengue virus in Aedes aegypti.* PLoS Pathogens. https://consensus.app/papers/details/w:f624f5d3-0949-4b8b-b2b9-7f28cf0bf0c2/
91. Rainey S. M. et al. (2016). *Wolbachia Blocks Viral Genome Replication Early in Infection without a Transcriptional Response by the Endosymbiont or Host Small RNA Pathways.* PLoS Pathogens. https://consensus.app/papers/details/w:b80cffb1-3e13-49cd-8339-ff02dfcd2c6f/
92. Dell'Aglio E. et al. (2025). *Epigenetics and non-coding RNAs in host–endosymbiont interactions: insights from Wolbachia and beyond.* Current Opinion in Insect Science. https://consensus.app/papers/details/w:31c974e8-a130-482e-9c9d-f395fd1dd9fc/
93. Ding Y. M. et al. (2026). *Host insect microRNAs are translocated to an obligate endosymbiont.* bioRxiv. https://consensus.app/papers/details/w:1dd7624d-dfe9-4f07-8ce2-b3db5e8e1d79/
94. Dou W. & Li (2026). *Wolbachia-induced cytoplasmic incompatibility triggers intergenerational dysregulation of the small RNA regulatory network in offspring.* Frontiers in Microbiology. https://consensus.app/papers/details/w:fa4f96b9-b3c4-4f53-a99c-cab92a1eb331/
95. Zhang L. et al. (2011). *Exogenous plant MIR168a specifically targets mammalian LDLRAP1.* Cell Research. https://consensus.app/papers/details/w:591f89a4-7713-4426-a497-0a5c065fb4b6/
96. Zhang Y. et al. (2012). *Analysis of plant-derived miRNAs in animal small RNA datasets.* BMC Genomics. https://consensus.app/papers/details/w:28b27f16-ac62-4177-a9fc-038c303e8cb3/
97. Tosar J. P. et al. (2014). *Mining of public sequencing databases supports a non-dietary origin for putative foreign miRNAs: underestimated effects of contamination in NGS.* RNA. https://consensus.app/papers/details/w:463e1b6f-e19b-4576-b241-f4a56b9130b2/
98. Lusk R. W. (2014). *Diverse and Widespread Contamination Evident in the Unmapped Depths of High Throughput Sequencing Data.* PLoS ONE. https://consensus.app/papers/details/w:9c925296-65e8-4c50-acd4-588f2f091071/
99. Kang W. et al. (2017). *Survey of 800+ data sets from human tissue and body fluid reveals xenomiRs are likely artifacts.* RNA. https://consensus.app/papers/details/w:0fe53ac7-b922-4602-b0b5-2f01dd92dbbb/
100. Kang W. et al. (2018). *miRTrace reveals the organismal origins of microRNA sequencing data.* Genome Biology. https://consensus.app/papers/details/w:e9947b04-bdd2-4628-87c6-c2dc50ab5e37/
101. Witwer K. W. (2018). *Alternative miRNAs? Human sequences misidentified as plant miRNAs in plant studies and in human plasma.* F1000Research. https://consensus.app/papers/details/w:23bbe746-80bc-4114-bd75-592cf396eda4/
102. Witwer K. W. & Hirschi K. D. (2014). *Transfer and functional consequences of dietary microRNAs in vertebrates: Concepts in search of corroboration.* BioEssays. https://consensus.app/papers/details/w:20d298b1-4777-4f99-aa58-d61f0c5c054c/
103. Heintz-Buschart A. et al. (2018). *Small RNA profiling of low biomass samples: identification and removal of contaminants.* BMC Biology. https://consensus.app/papers/details/w:ac605041-222f-47d0-8c8a-57f931997ecf/
104. Huang H. et al. (2018). *Extensive Degradation and Low Bioavailability of Orally Consumed Corn miRNAs in Mice.* Nutrients. https://consensus.app/papers/details/w:0f4f8612-0d93-4df5-820b-b3b8fa2647ba/
105. Huang H. et al. (2017). *An improved method to quantitate mature plant microRNA in biological matrices using modified periodate treatment and inclusion of internal controls.* PLoS ONE. https://consensus.app/papers/details/w:871690fd-8d89-42c6-8a53-f71fbc1fdc35/
106. Mar-Aguilar F. et al. (2020). *Evidence of transfer of miRNAs from the diet to the blood still inconclusive.* PeerJ. https://consensus.app/papers/details/w:30d6b0d4-091f-4713-b2fb-bd061086d470/
107. Tambaro F. et al. (2025). *The Conundrum of XenomiRs and Human Health.* Advances in Nutrition. https://consensus.app/papers/details/w:85e5b3ca-2f93-4345-9c19-ac9033ab631d/
108. Matsagani R. et al. (2026). *An Apple a Day? The Hypothesis of Cross-Kingdom Gene Regulation by Plant miRNAs in Mammals and Its Controversies.* IJMS. https://consensus.app/papers/details/w:f6c8ac27-f699-43a0-8d07-d0a030672df8/
109. Chen X. et al. (2021). *Large-scale identification of extracellular plant miRNAs in mammals implicates their dietary intake.* PLoS ONE. https://consensus.app/papers/details/w:e6386aa7-43f8-4572-afee-6554267b8588/
110. Zhao Q. et al. (2020). *Evidence for plant-derived xenomiRs based on a large-scale analysis of public small RNA sequencing data from human samples.* PLoS ONE. https://consensus.app/papers/details/w:53106b2a-4d5d-414a-ac06-d20fd948c3bb/
111. Luo Y. et al. (2017). *Detection of dietetically absorbed maize-derived microRNAs in pigs.* Scientific Reports. https://consensus.app/papers/details/w:e75ef8fd-66c3-4c36-97c8-e1c77f20e922/
112. Zoziuk M. O. et al. (2026). *Five Conserved microRNAs Dominate the Small-RNA Pool of Three Arid-Zone Camel-Forage Plants.* Genes. https://consensus.app/papers/details/w:7cdd555b-95c7-4fb7-bbb6-e24bbebb461e/
113. Costello M. et al. (2018). *Characterization and remediation of sample index swaps by non-redundant dual indexing on massively parallel sequencing platforms.* BMC Genomics. https://consensus.app/papers/details/w:e690897c-4afa-4ee2-a77e-03fdf096daaa/
114. MacConaill L. E. et al. (2018). *Unique, dual-indexed sequencing adapters with UMIs effectively eliminate index cross-talk.* BMC Genomics. https://consensus.app/papers/details/w:7a0de842-eb28-49b5-bb40-6b1bfede2996/
115. Farouni R. et al. (2020). *Model-based analysis of sample index hopping reveals its widespread artifacts in multiplexed single-cell RNA-sequencing.* Nature Communications. https://consensus.app/papers/details/w:a0c478ac-4447-4f1b-a5ae-3e36c3dd15a1/
116. Park S.-J. et al. (2019). *A systematic sequencing-based approach for microbial contaminant detection and functional inference.* BMC Biology. https://consensus.app/papers/details/w:84b44e67-8697-4989-9699-8deed0bf455d/
117. Strong M. J. et al. (2014). *Microbial Contamination in Next Generation Sequencing: Implications for Sequence-Based Analysis of Clinical Samples.* PLoS Pathogens. https://consensus.app/papers/details/w:b179cd52-c8d0-41c9-bf92-00b30d7d6515/
118. Salzberg S. L. et al. (2026). *Setting higher standards for reports of microbial species in human cancers.* Nature Cancer. https://consensus.app/papers/details/w:7b35b3bf-60e6-4d7e-b427-be66b764cad8/
119. Brumit D. F. et al. (2026). *Circulating Microbial DNA as a Potential Cancer Biomarker: Technical Challenges and Controlled Evaluation.* bioRxiv. https://consensus.app/papers/details/w:f5529099-d248-4699-82ca-1ccd274a3e03/
120. Beatty M. et al. (2014). *Small RNAs from plants, bacteria and fungi within the order Hypocreales are ubiquitous in human plasma.* BMC Genomics. https://consensus.app/papers/details/w:32cb7c88-7a47-40b4-8f83-15e0f681a293/
121. Wang K. et al. (2012). *The Complex Exogenous RNA Spectra in Human Plasma: An Interface with Human Gut Biota?* PLoS ONE. https://consensus.app/papers/details/w:3512c861-8c10-460f-9379-c79782163596/
122. Rozowsky J. S. et al. (2019). *exceRpt: A Comprehensive Analytic Platform for Extracellular RNA Profiling.* Cell Systems. https://consensus.app/papers/details/w:8bb127dc-5d6d-47fe-924a-3c176f160a72/

*Primary full-text used for the Weiberg 2013 numbers in §1.1, §3 and §6.1 (lengths, 5′U fraction, LTR locus breakdown, AGO1/AGO2/AGO4 co-IP): author manuscript at https://escholarship.org/content/qt3w66b649/qt3w66b649.pdf*
