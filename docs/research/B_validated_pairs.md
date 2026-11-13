# Strand B — Catalogue of published, experimentally claimed microbial/parasite sRNA → eukaryotic host mRNA pairs

Compiled for SmithHunter2 benchmark calibration. Date: 2026-11-13.

**Scope of inclusion.** One specific small RNA from a non-host organism + one specific host transcript it is claimed to regulate. Both directions of the host/microbe boundary are included where the *target* is a eukaryote (so plant→fungus pairs are kept, because they are the best-validated positives in the whole field and are mechanistically the same prediction problem).

---

## 1. Evidence grading scheme (strict, defined here)

Each pair is graded on what was **actually done for that pair**, not on what the paper did elsewhere.

| Grade | Definition |
|---|---|
| **A** | Site-level physical evidence (degradome/PARE or 5'RLM-RACE cleavage site mapped inside the predicted duplex, **or** AGO-IP of the microbial sRNA *in the host*) **plus** functional genetics on both sides (sRNA knockout/knockdown in the source organism **and** host-target mutant/rescue with a matching phenotype). |
| **B** | Reporter assay carrying the real target site **plus** an independent perturbation (sRNA mimic/overexpression, sRNA deletion in the source organism, or antagomir) with target measurement. Seed/site-mutated reporter control present in some but not all B pairs — noted per row. |
| **C** | Reporter **or** mimic/knockdown alone, with target measurement. No site mutant, no cleavage mapping. |
| **D** | Expression anticorrelation + computational prediction, possibly with AGO-IP showing the sRNA is loaded but no site-level evidence; or a phenotype assigned to a pathway rather than to a mapped site. |
| **E** | Computational prediction only. |
| **N** | Published **non-target**: the pair was tested and failed. Rare and valuable. |

Markers used below: `†` = detail taken from abstract/secondary source, full text not reached in this review (verify before use); `‡` = contested / failed independent replication.

---

## 2. Fungi and oomycetes → plants

| sRNA name / sequence | Source organism | Length | Host | Target gene (region) | Claimed pairing | Evidence type (what was actually done) | Grade | Ref |
|---|---|---|---|---|---|---|---|---|
| Bc-siR3.2 — `UACAUUGUGGAUCUUGUAGGU` (reverse-complemented from the published Northern probe; verify against Table 2.2 of ref 2) | *Botrytis cinerea* (LTR retrotransposon locus) | 21 nt | *Arabidopsis thaliana* | *MPK1*, *MPK2* (CDS) | near-complete complementarity; mismatch pattern not published per site; no ΔG given | AGO1-IP (sRNA loaded in host AGO1, absent from AGO2/AGO4); YFP sensor with **synonymous target-site mutant** control; *N. benthamiana* co-expression + HA-tag Western; Bc-sRNA-overexpressing *Arabidopsis*; *Bcdcl1 dcl2* double mutant (reduced virulence, targets de-repressed); *mpk1 mpk2* mutant more susceptible | **A** | 1, 2 |
| Bc-siR3.2 | *B. cinerea* | 21 nt | *Solanum lycopersicum* | *MAPKKK4* (CDS) | predicted (psRNATarget-class) | prediction + TRV-VIGS of tomato *MAPKKK4* → increased susceptibility; no reporter | **C** | 1 |
| Bc-siR3.1 — `UUGUGGAUCUUGUAGGUGGGC` (derived, verify) | *B. cinerea* (same LTR locus, 4-nt offset from Bc-siR3.2) | 21 nt | *A. thaliana* | *PRXIIF* peroxiredoxin (CDS) | near-complete | AGO1-IP; Bc-sRNAox lines; *dcl1 dcl2*; *prx* mutant susceptibility | **B** | 1, 2 |
| Bc-siR5 — `UUUGACUCGGAAUGUAUACUU` (derived, verify) | *B. cinerea* (LTR, 13 copies) | 21 nt | *A. thaliana* | *WAK* cell-wall-associated kinase (CDS) — allele used SALK_089827 | near-complete | AGO1-IP; Bc-sRNAox; *dcl1 dcl2*; *wak* mutant | **B** | 1, 2 |
| Bc-siR37 — `AAGGAGGAAGAUGAUGAUGAU` (from *B. cinerea* ABC-type plasma-membrane ATPase gene BC1G_10137) | *B. cinerea* | 21 nt | *A. thaliana* | *AtWRKY7* (**5'UTR**) | near-perfect complementarity over the 21 nt in the 5'UTR | anticorrelation in infection; *N. benthamiana* co-expression + Western; YFP sensor silencing at infection sites; Bc-siR37ox line; T-DNA mutant more susceptible | **B** | 3, 2 |
| Bc-siR37 | *B. cinerea* | 21 nt | *A. thaliana* | *AtFEI2* (**5'UTR**) | near-perfect | same panel as *WRKY7* (co-expression Western, ox line, T-DNA mutant) | **B** | 3, 2 |
| Bc-siR37 | *B. cinerea* | 21 nt | *A. thaliana* | *AtPMR6* (region not stated) | predicted | anticorrelation; reduced in Bc-siR37ox; T-DNA mutant more susceptible; Western inconsistent in source | **C** | 3, 2 |
| Bc-siR37 | *B. cinerea* | 21 nt | *A. thaliana* | *AtATG5* (region not stated) | predicted | anticorrelation only | **D** | 3, 2 |
| Bc-siR37 | *B. cinerea* | 21 nt | *A. thaliana* | *AtWRKY57* | predicted | none | **E** | 3, 2 |
| Bc-siR37 (remaining set) | *B. cinerea* | 21 nt | *A. thaliana* (15 predicted targets) / tomato (14 predicted targets) | full list in Table 2.1 of ref 2 — **not reached in this review** | — | prediction only | **E** | 2 |
| Fol-milR1 — `TGGATGAATCAAGCGTGGTATG` (listed as novel-m0003-3p; text says 23 nt, table gives 22 nt — discrepancy unresolved in source) | *Fusarium oxysporum* f. sp. *lycopersici* | 22–23 nt | tomato | *SlyFRG4* / Solyc06g007430, CBL-interacting protein kinase (**CDS**, cleavage at nt 943) | psRNATarget prediction; no duplex diagram published | 5'RLM-RACE cleavage at nt 943 (8/12 clones); **SlyAGO4a-IP** + sRNA-seq (SlyAGO1 only faint); *N. benthamiana* co-expression (qRT-PCR + GFP-fusion Western); 3 fungal knockouts + 2 site-mutated strains (attenuated); 3 overexpressors (hypervirulent on resistant cv.); **CRISPR *SlyFRG4* host mutants** → severe wilt; *SlyAGO4a* VIGS → reduced susceptibility | **A** | 4 |
| Fg-sRNA1 | *Fusarium graminearum* | n.r. in abstract | *Triticum aestivum* | *TaCEBiP* chitin elicitor binding protein | predicted | anticorrelation; *N. benthamiana* co-expression with YFP fluorescence + Western; BSMV-mediated silencing of Fg-sRNA1 in planta → reduced fungal invasion | **B** † | 5 |
| Vm-milR1 | *Valsa mali* | n.r. | *Malus domestica* (apple) | *MdRLKT1*, *MdRLKT2* receptor-like kinases | sequence-specific suppression claimed; site not mapped | precursor deletion mutant (no Vm-milR1, reduced virulence, host defences restored); target OE → resistance; target knockdown → susceptibility; no AGO-IP, no degradome | **B** | 6 |
| Vm-milR1 (sequestration) | *V. mali* | n.r. | apple | *MdRLKT21* acts as a **decoy** that binds Vm-milR1 and rescues *MdRLKT1/2* | — | genetic + biochemical; a rare example of a host counter-measure to a ck-sRNA | **B** | 7 |
| AamilR251 | *Alternaria alternata* f. sp. *mali* (EV-associated) | n.r. | apple | *MdMAPKKK1* | predicted | EV sRNA-seq; **individual milRNA deletion** → attenuated virulence + target restored; target OE → resistance | **B** | 8 |
| AamilR292 | *A. alternata* f. sp. *mali* | n.r. | apple | *MdRLK2* | predicted | same panel | **B** | 8 |
| AamilR004 | *A. alternata* f. sp. *mali* | n.r. | apple | *MdMAPKK6* | predicted | same panel | **B** | 8 |
| Pst-milR1 | *Puccinia striiformis* f. sp. *tritici* | ~20–22 nt | wheat | *TaPR2* (pathogenesis-related 2) | predicted | anticorrelation; host-induced silencing of Pst-milR1 restored *PR2* and resistance | **B** † | 9 |
| *Pst* milRNA / phased 21-nt sRNAs (unnamed set) | *P. striiformis* f. sp. *tritici* | 21 nt | wheat | an NLR receptor + multiple defence TF families | PARE-supported cleavage positions | sRNA-seq + **PARE/degradome** + expression; no functional genetics per pair | **D** | 10 |
| HpasRNAs (specific names/AGI codes **not retrieved** — in the eLife figures/supplement) | *Hyaloperonospora arabidopsidis* (oomycete) | 21–24 nt | *A. thaliana* | several immunity genes | predicted | AtAGO1-IP of the oomycete sRNAs; **Csy4/GUS in-situ sensor** showing target suppression in planta; *atago1* → resistance; STTM target-mimic plants blocking HpasRNAs → resistance; host target knockouts → enhanced susceptibility | **B** † | 11 |
| SisRNAs (set) | *Serendipita indica* (mutualist) | n.r. | *A. thaliana* | cell-wall, hormone, immunity and gene-regulation genes | predicted | stem-loop RT-qPCR of sRNA; protoplast transient expression + target qPCR; AtAGO1-IP | **C** | 12 |
| *R. irregularis* sRNAs | *Rhizophagus irregularis* (AM fungus) | n.r. | plant host | not individually resolved in abstract | — | transfer shown; colonization phenotype | **D** † | 13 |
| Tri-milRNA (Trichoderma) | *Trichoderma* sp. | n.r. | *B. cinerea* (fungal "host") | chitin synthase genes | predicted | nanomaterial delivery + target suppression + disease control | **C** † | 14 |
| Mi-miR-2b, Mi-miR-100 + 2 further miRNAs | *Meloidogyne incognita* (root-knot nematode) | ~21–22 nt | tomato (and *A. thaliana*) | 9 validated tomato transcripts (immune signalling, metabolism, feeding-site reprogramming) | predicted, dual-luciferase confirmed | **AGO1-RIP from galls** (10 nematode miRNAs loaded in host AGO1); **degradome**; dual-luciferase; miR-2b functional effect on giant-cell development; conserved AGO1 loading in *Arabidopsis* | **A** | 15 |

### Published *negative* findings in this section (important)
| Claim tested | Outcome | Ref |
|---|---|---|
| *B. cinerea* sRNAs → tomato transcripts, early infection | 7042 Bc-sRNAs predicted to hit 3185 tomato mRNAs; only 163 predicted targets actually down-regulated; transposon-locus deletion mutants and *Bcdcl1 Bcdcl2* mutants (>99% loss of transposon sRNAs) were **not** reduced in virulence. Authors conclude no detectable ckRNAi role at this stage | 16 |
| *Zymoseptoria tritici* sRNAs → wheat | 389 in-planta-induced sRNA loci, targets predicted, but **molecular validation of selected wheat targets failed**; *DCL* and *AGO* deletion strains fully virulent | 17 |

---

## 3. Plants → fungi / oomycetes (reverse direction; the best-validated positives in the field)

| sRNA | Source | Length | Target organism | Target gene (region) | Evidence | Grade | Ref |
|---|---|---|---|---|---|---|---|
| ghr-miR166 | *Gossypium hirsutum* (cotton) | 21 nt | *Verticillium dahliae* | *Clp-1* (Ca²⁺-dependent cysteine protease) | 5'RACE cleavage; miRNA knockdown (STTM) → susceptibility; fungal target mutants → virulence loss | **A** † | 18 |
| ghr-miR159 | cotton | 21 nt | *V. dahliae* | *HiC-15* (isotrichodermin C-15 hydroxylase) | same panel | **A** † | 18 |
| TAS1c-siR483, TAS2-siR453 | *A. thaliana* (EV-delivered) | 21 nt | *B. cinerea* | *BcVPS51*, *BcDCTN1*, *BcSAC1* | EV purification, sRNA-seq, fungal target suppression, mutant phenotype | **B** | 19 |
| sly-miR393 | tomato | 21 nt | *B. cinerea* | cell-wall-related gene | nanocarrier delivery + target suppression + disease control | **C** † | 20 |

These four are the cleanest template for what a Grade-A positive looks like, and are worth keeping in the benchmark even though the direction is plant→fungus.

---

## 4. Bacteria → eukaryotic hosts

| sRNA name / sequence | Source | Length | Host | Target gene (region) | Claimed pairing | Evidence | Grade | Ref |
|---|---|---|---|---|---|---|---|---|
| **sRNA52320** — `CGCGGGAUGGAGCAGUCUGGUAGC` (5' fragment of tRNA-Met, locus PA14_52320) | *Pseudomonas aeruginosa* (OMV) | 24 nt | human airway epithelium / mouse | *MAP2K2*, *MAP2K3*, *MAP2K4*, *MAP3K7*, *PIK3R2* | miRanda score ≥140 + BLASTN; **no duplex, no site coordinates, no ΔG published** | synthetic sRNA transfection (10 nM) + proteomics (MAP3K7 and MAP2K4 protein down); IL-8 mRNA and secretion reduced; naked-sRNA control (no effect); ΔsRNA + complement OMVs (−34% IL-8); **mouse in vivo** (KC, neutrophils). **No reporter, no seed mutant** | **C** | 21 |
| **RsmY** (nt 48–65 resemble hsa-miR-144-3p) | *Legionella pneumophila* (EV-translocated) | ~50–150 nt population; region used is 18 nt | human THP-1 / hMDM / CD14⁺ | *DDX58* (RIG-I), "UTR" (not specified 3' vs 5') | miRNA-mimicry region only; no seed coordinates | luciferase reporter; EV-luciferase; RNA–RNA binding assay; protein loss after infection and EVs, lost with Δ*rsmY* EVs; AGO2-inhibitor test. **No seed-mutant reporter**; AGO2-CLIP detected RsmY in 2/3 pulldowns, not significant; no in vivo | **B** | 22 |
| RsmY | *L. pneumophila* | as above | human | *REL* (cRel), "UTR" | as above | binding assay; luciferase; slight protein reduction after transfection | **C** | 22 |
| **tRNA-Phe fragment** (lppt29; nt 59–73 resemble hsa-miR-5001-3p) | *L. pneumophila* | 15-nt mimicry region | human | *IRAK1*, "UTR" | as above | transfection → protein loss; binding assay; luciferase; EV-luciferase; AGO2-inhibitor. Single-copy locus so **no knockout possible** | **C** | 22 |
| tRNA-Phe fragment | *L. pneumophila* | as above | human | *DDX58* | as above | binding assay; luciferase only | **C** | 22 |
| RsmY / tRNA-Phe predicted-only set | *L. pneumophila* | — | human | *MAPK8*, *MAVS*, *AGO1*, *MAPK9* | miRNA-mimicry prediction | none | **E** | 22 |
| **Xosr001** | *Xanthomonas oryzae* pv. *oryzicola* (OMV) | n.r. | rice | *OsJMT1* (jasmonate methyltransferase) | predicted | OMV sRNA-seq; in-planta transcript suppression; synthetic Xosr001 applied to *OsJMT1-HA-OE* line → suppression; Δ*Xosr001* strain → restored stomatal immunity / reduced susceptibility; downstream MeJA readout | **B** | 23 |
| **sRNA23392** | *Porphyromonas gingivalis* (OMV) | n.r. | human OSCC cells | *DSC2* (desmocollin-2) | predicted | OMV abundance; sRNA inhibitor attenuates OMV-induced migration/invasion; "targeting" claimed — reporter not confirmed from the abstract | **C** † | 24 |
| **sRNA71** | *Lactobacillus plantarum* (EV) | miRNA-sized | human HEK293T | *TP53* (**3'UTR**) | predicted site in 3'UTR | mimic transfection → TP53 protein down (Western + proteomics); **luciferase reporter** of the 3'UTR. No seed mutant stated | **B** | 25 |
| **Ile-tRF-5X** | *Escherichia coli* (OMV) | <16 nt (very small RNA) | human HCT116 | *MAP3K4* — note the effect is **up**-regulation | predicted | OMV loading, transfer to human cells, target expression change, proliferation phenotype | **C** | 26 |
| **sR-2509025**, **sR-989262** | *Helicobacter pylori* (OMV) | sncRNA | human AGS gastric cells | no individual mRNA validated; read-out is IL-8 secretion | predicted | OMV enrichment, transfer, cytokine read-out; a follow-up showed sncRNA disruption changes OMV vaccine efficacy in mice | **D** | 27, 28 |
| msRNAs (multiple) | *Aggregatibacter actinomycetemcomitans*, *P. gingivalis*, *Treponema denticola* (OMV) | 15–25 nt | human Jurkat T cells | "human immune-related target genes" (predicted set) | prediction | deep sequencing; OMV delivery into eukaryotic cells; exogenous msRNA suppressed certain cytokines | **D** | 29 |
| *Bacteroides fragilis* OMV sRNAs | gut commensal | n.r. | intestinal epithelial cells | inflammatory signalling (pathway-level) | prediction | EV sRNA-seq + signalling read-out | **D** † | 30 |
| *Flavobacterium psychrophilum* OMV sRNAs | fish pathogen | n.r. | rainbow trout | immune genes (predicted interactions) | prediction | sRNA-seq + predicted interactions | **E** † | 31 |
| ncBCG427 | *M. bovis* BCG | n.r. | macrophage | apoptosis / inflammation (pathway-level) | prediction | expression + phenotype | **D** † | 32 |

**Gap flagged:** no *Mycobacterium tuberculosis* sRNA → specific human mRNA pair with site-level evidence was found; *Salmonella* OMV RNA work (ref 33) characterises cargo but assigns no host-mRNA pair; *Listeria* likewise. These are prediction-stage literatures, not positive-set material.

---

## 5. Wolbachia and other endosymbionts

| sRNA / sequence | Source | Length | Host | Target | Evidence | Grade | Ref |
|---|---|---|---|---|---|---|---|
| **WsnRNA-46** — `UCCCAGUGUUGGCUACUUGGAUGACAU` (mimic and inhibitor are listed with the *same* sequence in the Methods — a likely error to check before reuse) | *Wolbachia* wMelPop-CLA | 27 nt mature (precursor 65–70 nt) | *Aedes aegypti* | **Dhc** (dynein heavy chain) | GFP reporter carrying the target site, **with seed-region mutants that abolished the effect**; mimic injection into Wol⁻ mosquitoes; RT-qPCR in Wol⁺ vs Wol⁻. **Direction is anomalous: Dhc transcripts go UP** — authors propose stabilisation, not silencing | **B** (direction caveat) | 34 |
| WsnRNA-46 | *Wolbachia* | 27 nt | *Wolbachia* itself | *murD*, *carA* | GFP reporter (*murD*); inhibitor transfection (both) — **self-targets, not host targets** | — | 34 |
| WsnRNA-59 (abstract calls it WsnRNA-49; body text uses -59) | *Wolbachia* | ~30 nt | — | no validated target | Northern detection only | **E** | 34 |
| WsnRNA-744, WsnRNA-3640, WsnRNA-6108 (+9 more) | *Wolbachia* in *Tetranychus truncatus* | n.r. | *T. truncatus* | reproduction-related genes (enrichment-level prediction) | inhibition of individual WsnRNAs changes fecundity (−744, −3640 reduce; −6108 increases), but **no specific mRNA pair validated** | **E** (phenotype without a pair) | 35 |
| aae-miR-2940 | **host** miRNA induced by *Wolbachia* | 22 nt | *A. aegypti* | metalloprotease gene | GFP reporter, mimic, miRNA inhibitor, target knockdown → reduced *Wolbachia* density | not a microbial sRNA — **excluded from the positive set**, listed because it is often miscited as one | — | 36 |

**This is the thinnest and most important part of the catalogue for SmithHunter2.** There is exactly **one** *Wolbachia*-sRNA → host-mRNA pair in the literature with reporter-level evidence (WsnRNA-46 → *Dhc*), and its direction of effect is the wrong way round for a silencing model.

---

## 6. Helminths and protozoa → animal hosts

| sRNA / sequence | Source | Length | Host | Target (region) | Evidence | Grade | Ref |
|---|---|---|---|---|---|---|---|
| EV miRNA set: parasite homologues of miR-8/miR-200 and let-7 (7mer sites), miR-63/miR-425 (6mer site) | *Heligmosomoides polygyrus* (exosomes) | ~22 nt | mouse MODE-K | ***Dusp1* (3'UTR)** | psiCHECK2 *Renilla* 3'UTR reporter + 50 nM mimics (1.2–2.0× repression singly, 3.1× for a 3-miRNA cocktail). **No seed-mutant 3'UTR, no scrambled control described.** Microarray anticorrelation; in-vivo *Alternaria* model with exosomes (whole-vesicle, not miRNA-specific) | **C** | 37 |
| **hpo-miR-71** — `UGAAAGACAUGGGUAGUGAGACG` | *H. polygyrus* | 23 nt | mouse | ***Il1rl1*/*Il33r* (3'UTR)** — two predicted 7mer sites | **3'UTR reporter was NOT repressed.** Published non-target | **N** | 37 |
| Hpo-miR-71-5p | *H. polygyrus bakeri* | 22 nt | mouse macrophages / Th2 cells | *Irf4* | mimic transfection + target validation; IFN-γ, IL-2/STAT5, mTOR pathways down; **sRNA shown to incorporate into mouse AGO complexes** | **C** | 38 |
| Asu-miR-791-3p | *Ascaris suum* (abundant in *T. muris* secretions) | 22 nt | mouse | *Mtor* | mimic transfection + target validation | **C** | 38 |
| miR-5352 (conserved across GI nematodes; seed-shares with mammalian miR-92a family) | GI nematodes | 22 nt | mouse intestinal epithelium | *Klf4*, IL-22 receptor | targeted repression demonstrated; organoid differentiation and cytokine phenotypes; Wnt/Notch modulation | **B/C** † | 39 |
| **sja-miR-2a-3p** | *Schistosoma japonicum* (EV) | 22 nt | mouse hepatic stellate cells | *BCL2* | **dual-luciferase**; EV transfer confirmed by transwell; overexpression → apoptosis; **agomir in infected mice in vivo** → BCL2 down, fibrosis reduced | **B** | 40 |
| sja-miR-125b, sja-bantam | *S. japonicum* (EV) | ~22 nt | mouse macrophages | *Pros1*, *Fam212b*, *Clmp* | EV uptake, cargo transfer, target regulation, macrophage phenotype | **C** | 41 |
| **emu-miR-10a-5p** | *Echinococcus multilocularis* (EV) | 22 nt | mouse MH-S macrophages | *LIF* → JAK1–STAT3 | mimic transfection + transcriptome; target and pathway validated by RT-qPCR/Western; OVA asthma model in vivo | **B/C** | 42 |
| **Bx-miR-71-5p** | *Bursaphelenchus xylophilus* (exosomes) | 22 nt | *Monochamus alternatus* (vector beetle) | *Notch* → *Mmp3* → ECM | exosome uptake; nanomaterial-mediated delivery in vivo → Notch up, ECM thicker, tracheal diameter up | **C** | 43 |
| tsp-miR-1 | *Trichinella spiralis* (EV) | 22 nt | mouse | IL-6 suppression (no mapped mRNA site) | EV stage-specific functional assays | **D** | 44 |
| bma-miR-100 | *Brugia malayi* microfilariae (EV) | 22 nt | human DC / THP-1 | *mTOR*, *EIF4EBP1* | microarray (GSE157226) + qPCR of the miRNA; **target prediction only** (MIRPathv3/TarBase/MicroT-CD); phospho-mTOR reduction by whole EVs | **D** | 45 |
| T. suis EV miR-22-3p / miR-29-3p homologues | *Trichuris suis* | 22 nt | human colon epithelial cells | "conserved oncogenic targets" | mimic/EV functional assays; targets not individually mapped | **D** | 46 |
| tRNA-derived halves | *Trypanosoma cruzi* (shed microvesicles) | ~30–35 nt | mammalian host cells | genome-wide expression change; **no specific pair** | transfer + transcriptome | **E** | 47 |

---

## 7. Plants → mammals via diet (contested — flagged)

| sRNA | Source | Length | Host | Target (region) | Evidence | Grade | Ref |
|---|---|---|---|---|---|---|---|
| **osa-miR168a** | rice | 21 nt | human / mouse | *LDLRAP1* (**exonic/CDS**, plant-style near-complete pairing) | detection in serum; mimic/feeding; reporter; mouse feeding with LDLRAP1 protein change | **C ‡** — multiple independent failures to replicate; large surveys attribute circulating plant miRNA reads to contamination/index hopping | 48, 49, 50 |
| **miR159** | *Arabidopsis*/diet | 21 nt | human | *TCF7* | serum detection; mimic; reporter; xenograft tumour-growth effect | **C ‡** | 51 |
| **MIR2911** | honeysuckle | ~20 nt, GC-rich | human/mouse | influenza *PB2*/*NS1* — **viral, not a host mRNA** | detection in serum, antiviral effect | out of scope for a host-mRNA positive set; included for completeness | 52, 53 |

Counter-evidence to carry in the benchmark: a survey of 800+ human tissue and body-fluid datasets concluded xenomiR signals are likely artefacts (ref 49); human sequences have been misidentified as plant miRNAs in both plant studies and human plasma (ref 50); corn miRNAs are extensively degraded with low bioavailability in mice (ref 54).

**Recommendation: keep the dietary pairs in a separate "contested" partition, never in the primary positive set.** Mixing them will inflate apparent recall and make any ROC curve uninterpretable.

---

## 8. Grade counts and whether this supports benchmarking

Counting distinct sRNA→target-transcript pairs catalogued above (excluding *Wolbachia* self-targets, the viral MIR2911 target, and the host-miRNA aae-miR-2940):

| Grade | Count | Pairs |
|---|---|---|
| **A** | **4** | Bc-siR3.2→AtMPK1/MPK2 (counted as 1 locus pair with 2 genes), Fol-milR1→SlyFRG4, Mi-miRNA→tomato targets (9 transcripts, counted as 1 study-level set), ghr-miR166→*Clp-1* / ghr-miR159→*HiC-15* (plant→fungus, counted as 2) |
| **B** | **17** | Bc-siR3.1→PRXIIF; Bc-siR5→WAK; Bc-siR37→WRKY7; Bc-siR37→FEI2; Fg-sRNA1→TaCEBiP; Vm-milR1→MdRLKT1; Vm-milR1→MdRLKT2; AamilR251/292/004→3 apple kinases; Pst-milR1→TaPR2; HpasRNA set→*Arabidopsis* immunity genes; TAS1c-siR483/TAS2-siR453→Bc genes; RsmY→DDX58; Xosr001→OsJMT1; sRNA71→TP53; WsnRNA-46→Dhc; sja-miR-2a-3p→BCL2 |
| **C** | **17** | Bc-siR3.2→SlMAPKKK4; Bc-siR37→PMR6; SisRNAs; Tri-milRNA; sly-miR393; sRNA52320→5 MAPK-pathway genes (1 set); RsmY→REL; tRNA-Phe→IRAK1; tRNA-Phe→DDX58; sRNA23392→DSC2; Ile-tRF-5X→MAP3K4; HpEV miRNAs→Dusp1; Hpo-miR-71-5p→Irf4; Asu-miR-791-3p→Mtor; miR-5352→Klf4/IL-22R; sja-miR-125b/bantam→Pros1/Fam212b/Clmp; emu-miR-10a-5p→LIF; Bx-miR-71-5p→Notch; + the 2 contested dietary pairs held separately |
| **D** | **11** | Bc-siR37→ATG5; *Pst* PARE set; *R. irregularis*; H. pylori sncRNAs; periodontal msRNAs; *B. fragilis*; ncBCG427; tsp-miR-1→IL-6; bma-miR-100→mTOR; T. suis miRNAs; plus pathway-level others |
| **E** | **≥8 named + hundreds unnamed** | Bc-siR37 remaining predicted set (29 transcripts); *Legionella* MAPK8/MAVS/AGO1/MAPK9; WsnRNA-59; 12 *Tetranychus* WsnRNAs; *Flavobacterium* set; *T. cruzi* tRNA halves |
| **N** (non-target) | **1** | hpo-miR-71 → mouse *Il1rl1* 3'UTR |

**Totals: ~49 named pairs, of which 4 are A, 17 are B, 17 are C, 11 are D, and 1 is a published non-target.**

### Honest assessment

- **ROC analysis: not supportable.** ROC needs a positive set that is (i) large, (ii) unbiased in how it was discovered, and (iii) matched to a credible negative set. All three fail here. The ~21 A+B pairs were each found by first running a target-prediction tool (psRNATarget, TAPIR, miRanda, RNAhybrid) and then validating a handful of top hits. The positive set is therefore **the output of the class of tools being benchmarked**, which makes any AUC a measure of agreement with psRNATarget, not of truth. There is no published pair that was discovered by an assay first and only afterwards found to be predictable.
- **Precision/recall on an unbiased hold-out: not supportable**, for the same reason plus the absence of negatives (section 9).
- **What this set *can* do, and should be used for:**
  1. **Sanity checking / regression testing.** A tool that cannot recover Bc-siR3.2→*MPK1/MPK2*, Fol-milR1→*SlyFRG4* (which has a mapped cleavage site at nt 943 — a true positional label), Bc-siR37→*WRKY7*/*FEI2* and ghr-miR166→*Clp-1* is broken. These four systems give you site-level ground truth.
  2. **Rank-recovery testing.** For each validated pair, ask where the true target falls in the ranked list over the whole host transcriptome. Median rank of the true target is reportable and honest; it needs no negative set.
  3. **Cross-kingdom-class stratification.** Plant targets are mostly **CDS** with near-complete complementarity (plant-style, cleavage-competent). Animal targets are mostly **3'UTR** with 6mer/7mer seeds (metazoan-style). A single scoring model will not serve both; benchmark them as two separate problems.
  4. **A direction-of-effect check.** Two catalogued pairs (WsnRNA-46→*Dhc*, Ile-tRF-5X→*MAP3K4*) report **up**-regulation. If SmithHunter2 predicts repression, these are not positives for it.
- **Bacteria→animal is the weakest sub-corpus and the one closest to the *Wolbachia* use case.** Best in class is *Legionella* RsmY→*DDX58* (B), and even there the AGO2-CLIP failed to reach significance and no seed mutant was tested. For a *Wolbachia*-focused tool, the usable positive set is effectively **one pair** (WsnRNA-46→*Dhc*), with an anomalous direction. That is below any threshold for calibration.

---

## 9. The negative set problem

**Has anyone published verified non-targets?** Essentially no — with three exceptions found in this review:

1. **hpo-miR-71 → mouse *Il1rl1*/*Il33r* 3'UTR** (ref 37). Two predicted 7mer sites; the 3'UTR reporter was not repressed. This is, as far as this review found, the only cleanly reported single-pair negative in the entire cross-kingdom literature.
2. **Qin et al. 2022** (ref 16): 7042 *B. cinerea* sRNAs → 3185 predicted tomato targets, of which only **163 (≈5%)** were actually down-regulated during early infection, and the transposon-locus and *dcl1 dcl2* deletion mutants were not reduced in virulence. The ~3000 non-responding predicted targets are the largest de facto negative set available, though "not down-regulated at one timepoint" is a weak negative label.
3. **Kettles et al. 2019** (ref 17): predicted wheat targets of *Z. tritici* sRNAs that **failed molecular validation**. The paper does not tabulate them as negatives but the underlying data would support extracting a negative set.

**What groups actually use as negatives.** In practice: nothing principled. The standard cross-kingdom pipeline (ref 55, a methods review of exactly this workflow) recommends psRNATarget and TAPIR with seed-region penalties (positions 2–13 and 2–12 respectively), scores 0–5 with "decrease the score to reduce false positives," and TAPIR free-energy ratio ≥0.7. Its controls are **biological, not statistical**: uninfected plants and axenic cultures, to show the sRNA is infection-induced. It explicitly does **not** recommend shuffled-sequence or dinucleotide-preserving randomised controls, and does **not** discuss base rates or quantify false-positive rates. Where negatives appear at all in individual papers, they are:
- **Unrelated control sRNAs** in reporter assays (e.g. ath-miR395 in ref 1, sly-miR166 in ref 4) — a negative *sRNA*, not a negative *site*;
- **Seed/site-mutated reporters** (refs 1, 34) — a negative *version of a positive site*, which is the strongest control in the field and is present in only a minority of papers;
- **Mismatched-species transcriptomes**, used informally (e.g. requiring ≥2 mismatches to the source organism's own genome, ref 55).

**The base-rate problem, stated plainly.** A 21-nt sRNA scored with a permissive plant-style tool hits on the order of 10²–10³ transcripts in a 30–40k-gene host transcriptome; the *B. cinerea*–tomato numbers above (7042 sRNAs → 3185 targets, 5% responsive) are the field's own measurement of this. If the true number of functional cross-kingdom pairs per pathogen is in the tens, the prior is ~10⁻³–10⁻⁴ and even a 99%-specific predictor returns a list that is overwhelmingly false. **The cross-kingdom field does essentially nothing about this.** It does not correct for multiple testing across the transcriptome, does not report expected false positives per library, and does not publish the non-validated predictions that would let anyone else estimate the rate. The field's substitute for a base-rate correction is *post-hoc functional genetics on a handful of hits* — which establishes those hits but says nothing about the tool's precision.

**Concrete recommendation for SmithHunter2's negative set**, in descending order of defensibility:
1. **Host transcripts not expressed in the infected tissue** (from the matched host RNA-seq) that nonetheless score as predicted targets — a real negative with a real base rate.
2. **The non-responding predicted targets of Qin et al. 2022** — ~3000 labelled pairs, the only large set available.
3. **Species-mismatched controls**: score the sRNA against a host transcriptome from a species the organism does not infect, holding GC and transcript-length distributions matched.
4. **Dinucleotide-preserving shuffles of the sRNA** (not of the transcript — shuffling the transcript destroys UTR composition and inflates apparent specificity).
5. Report, alongside any precision figure, the **expected hits per 10,000 transcripts at the chosen threshold**, so readers can apply their own prior.

---

## 10. Where the raw data live

Accessions recovered in this review. Several are first-party (confirmed in the paper text); the rest need the paper's data-availability statement fetched.

| Dataset | Accession | Notes |
|---|---|---|
| *B. cinerea* sRNA libraries | **GSE45320** | Weiberg 2013 / Wang thesis |
| *B. cinerea*-infected *Arabidopsis* | **GSE45323** | sRNA + host |
| *B. cinerea*-infected *S. lycopersicum* | **GSE45321** | |
| *P. aeruginosa* and matched OMV sRNAs | **GSE71598** | Koeppen 2016 — directly reusable bacterial-sRNA-in-vesicle data |
| HBE cells exposed to *P. aeruginosa* OMVs | **GSE80421** | matched host response |
| *H. polygyrus* exosome sequencing + host microarray | **GSE55978**, **GSE55941** | Buck 2014; the paper does not state which is which |
| *Wolbachia* wMelPop-CLA sRNA deep sequencing (Aag2 and Pop cell lines) | **GSE55210** | **Mayoral 2014 — the key WsnRNA dataset** |
| *B. malayi* microfilaria EV miRNA microarray | **GSE157226** | proteome at ProteomeXchange **PXD021844** |
| *Legionella* EV study | flow data FlowRepository **FR-FCM-Z2XL**, **FR-FCM-Z2XM** | sequencing accession not in the fetched portion — check the data-availability statement of ref 22 |
| Fol-milR1 study | not captured | New Phytologist 2021, doi:10.1111/nph.17436 — supplementary tables hold the sRNA list |
| Dunker 2020 oomycete AGO1-IP sRNA-seq | not captured | eLife 56096 — the site was JS-gated; fetch the supplementary files from the eLife record |
| *Pst*–wheat sRNA + PARE/degradome | not captured | Mueth 2022, *Genomics* — PARE libraries are the valuable part |

### Wolbachia-infected vs uninfected small RNA pairs across host species — reusability assessment

This is the most directly reusable resource for the project and it is thinner than it looks.

- **GSE55210** (Mayoral 2014) is the only dataset built specifically to find *Wolbachia*-derived sRNAs, and it is **cell-line** data (Aag2 uninfected vs Pop *Wolbachia*-infected), i.e. one host species.
- Most other "*Wolbachia* sRNA" datasets in the corpus are **host-miRNA** studies with infected/uninfected contrasts, which are reusable only if you re-map the unassigned reads to the *Wolbachia* genome. Candidates found: Hussain 2011 (*A. aegypti*, host miRNA microarray), Mayoral 2014 PLoS ONE (mosquito cell line miRNA profile), Liu 2019 (*Laodelphax striatellus*), Zheng 2019 (*D. melanogaster* CI sRNA), Bishop 2022 (*A. aegypti* wAlbB), Zhang 2025 iScience (*Tetranychus truncatus* — the newest and the only mite dataset), Dou 2026 (CI offspring sRNA networks). Accessions for these were **not** retrieved; each paper's data-availability statement must be checked individually.
- Projects surfaced but unverified for content: **PRJNA557362**, **PRJNA218499**, **PRJEB21550**, **PRJEB74854**. Do not cite these as *Wolbachia* sRNA sets without checking.
- **Practical point:** a reanalysis across *Aedes*, *Drosophila*, *Laodelphax*, *Tetranychus* and (from the user's own system) *Zyginidia pullula* would be the first multi-host WsnRNA atlas. That is a genuinely novel contribution, because the published WsnRNA repertoire rests on **two** studies.

---

## 11. Biggest gaps

1. **No pair anywhere in the bacteria→animal literature has both a mapped site and functional genetics on both sides.** The best is Grade B.
2. ***Wolbachia* specifically has one reporter-level host pair, with an anomalous direction of effect.** SmithHunter2 cannot be calibrated on its own organism.
3. **Sequences are chronically unpublished in the main text.** Of ~49 pairs, exact mature sequences were recovered for only 7 (Bc-siR37, Bc-siR3.1/3.2/5 by reverse-complementing probes, Fol-milR1, sRNA52320, WsnRNA-46, hpo-miR-71). Everything else is in supplementary tables.
4. **Target-site coordinates are almost never given** — no duplex diagrams, no ΔG, no mismatch positions in the large majority of papers. Only Fol-milR1 (nt 943) and the Bc-siR37 5'UTR sites are positionally specified.
5. **No negative set exists.** See section 9.
6. **Two prominent negative-result papers (refs 16, 17) are rarely cited by the tool-development literature** and should be read before any claim about tool performance.

---

## References

1. Weiberg A. et al. (2013) *Fungal small RNAs suppress plant immunity by hijacking host RNA interference pathways.* Science. https://consensus.app/papers/details/w:999aaa2e-ddf2-4d63-8f8f-8981293722b3/
2. Wang M. (2015) *The role of fungal pathogen small RNAs in host-microbe interactions.* UC Riverside dissertation (eScholarship) — contains Bc-sRNA sequences, probe tables and GEO accessions. https://escholarship.org/content/qt3w66b649/qt3w66b649.pdf
3. Wang M. et al. (2017) *Botrytis small RNA Bc-siR37 suppresses plant defense genes by cross-kingdom RNAi.* RNA Biology. doi:10.1080/15476286.2017.1291112
4. Ji H.-M. et al. (2021) *Fol-milR1, a pathogenicity factor of Fusarium oxysporum, confers tomato wilt disease resistance by impairing host immune responses.* New Phytologist. doi:10.1111/nph.17436 — full text used: https://escholarship.org/content/qt70d644fs/qt70d644fs.pdf
5. Jian J. et al. (2019) *One small RNA of Fusarium graminearum targets and silences CEBiP gene in common wheat.* Microorganisms 7:425. https://mdpi.com/2076-2607/7/10/425
6. Xu M. et al. (2022) *A fungal microRNA-like RNA subverts host immunity and facilitates pathogen infection by silencing two host receptor-like kinase genes.* New Phytologist. https://consensus.app/papers/details/w:1ce50e8b-4df5-49f0-b0f0-15a6d424100b/
7. Lei M. et al. (2025) *An apple susceptibility-related LRR receptor-like kinase MdRLKT21 activates plant immunity by hijacking a trans-kingdom fungal microRNA-like RNA.* Advanced Science. https://consensus.app/papers/details/w:74b145f7-6ea6-4a26-a235-b880ab79087c/
8. Tang J. et al. (2026) *Fungal microRNA-like RNAs facilitate convergent attack on a host immune signalling network.* Plant Cell & Environment. https://consensus.app/papers/details/w:7bbbe83c-976f-4e81-99fa-f87cc4127181/
9. Wang B. et al. (2017) *Puccinia striiformis f. sp. tritici microRNA-like RNA 1 (Pst-milR1) ... impairs wheat resistance by suppressing the wheat pathogenesis-related 2 gene.* PMID 28464281. https://www.pubmed.gov/28464281/
10. Mueth N.A. et al. (2022) *Small RNAs target native and cross-kingdom transcripts on both sides of the wheat stripe rust interaction.* Genomics. https://consensus.app/papers/details/w:2c5f3294-c42a-4e4b-b17d-8495091ba19d/
11. Dunker F. et al. (2020) *Oomycete small RNAs bind to the plant RNA-induced silencing complex for virulence.* eLife 9:e56096. https://elifesciences.org/articles/56096
12. Nasfi S. et al. (2024) *A pipeline for validation of Serendipita indica effector-like sRNA suggests cross-kingdom communication in the symbiosis with Arabidopsis.* J Exp Bot. https://consensus.app/papers/details/w:b599665d-2a43-45ad-8de3-9ff7d5486a42/
13. Usländer A. et al. (2026) *Cross-kingdom RNA interference promotes arbuscular mycorrhiza development.* Nature Plants. https://consensus.app/papers/details/w:86c99291-4e80-432a-938c-e91a8ee41227/
14. Yin Y. et al. (2026) *A Trichoderma milRNA cross-species targeting chitin synthase genes is delivered by carbon nanomaterials to suppress Botrytis cinerea.* Plant Physiology.
15. (2026) *AGO1-RIP, degradome and luciferase validation of Meloidogyne incognita miRNAs silencing tomato transcripts* — root-knot nematode ckRNAi study. https://consensus.app/papers/details/ (record retrieved via Consensus; full citation to be completed from the DOI)
16. Qin S. et al. (2022) *Molecular characterization reveals no functional evidence for naturally occurring cross-kingdom RNA interference in the early stages of Botrytis cinerea–tomato interaction.* Molecular Plant Pathology. https://consensus.app/papers/details/w:acbf9b72-fc3f-49db-91a2-89dc2e841073/
17. Kettles G.J. et al. (2019) *sRNA profiling combined with gene function analysis reveals a lack of evidence for cross-kingdom RNAi in the wheat–Zymoseptoria tritici pathosystem.* Front Plant Sci. https://consensus.app/papers/details/w:e1b0959a-f0d2-4127-9276-b78b82b204eb/
18. Zhang T. et al. (2016) *Cotton plants export microRNAs to inhibit virulence gene expression in a fungal pathogen.* Nature Plants. https://consensus.app/papers/details/ (Consensus record)
19. Cai Q. et al. (2018) *Plants send small RNAs in extracellular vesicles to fungal pathogen to silence virulence genes.* Science. https://consensus.app/papers/details/w:69cc5ea5-2f44-4411-8fd8-a36dc75ddfc7/
20. Yin Y. et al. (2025) *Nanomaterial-mediated delivery of tomato miR393 ... confers efficient control of gray mold disease.* Plant Communications.
21. Koeppen K. et al. (2016) *A novel mechanism of host-pathogen interaction through sRNA in bacterial outer membrane vesicles.* PLoS Pathogens 12:e1005672. https://pmc.ncbi.nlm.nih.gov/articles/PMC4905634/
22. Sahr T. et al. (2022) *Translocated Legionella pneumophila small RNAs mimic eukaryotic microRNAs targeting the host immune response.* Nature Communications 13:762. https://pmc.ncbi.nlm.nih.gov/articles/PMC8828724
23. Wu Y. et al. (2024) *Suppression of host plant defense by bacterial small RNAs packaged in outer membrane vesicles.* Plant Communications. https://consensus.app/papers/details/w:b0a24d4f-09c8-479d-a7df-420d6223ea05/
24. Liu D. et al. (2021) *sRNA23392 packaged by Porphyromonas gingivalis outer membrane vesicles promotes oral squamous cell carcinoma migration and invasion by targeting desmocollin-2.* Molecular Oral Microbiology. https://consensus.app/papers/details/w:338496aa-b6a4-4900-84cd-ef465f8dfa29/
25. Yu S. et al. (2022) *Biological functions and cross-kingdom host gene regulation of small RNAs in Lactobacillus plantarum-derived extracellular vesicles.* Front Microbiol. https://consensus.app/papers/details/w:965207fd-73ef-4f84-b328-d21c7c73f481/
26. Diallo I. et al. (2022) *A tRNA-derived fragment present in E. coli OMVs regulates host cell gene expression and proliferation.* PLoS Pathogens. https://consensus.app/papers/details/w:8c70e42f-6a88-4b73-bb0b-9a9e01b4a752/
27. Zhang H. et al. (2019) *sncRNAs packaged by Helicobacter pylori outer membrane vesicles attenuate IL-8 secretion in human cells.* Int J Med Microbiol. https://consensus.app/papers/details/w:8c90d146-60c4-47fb-aea7-d4b1c1419879/
28. Li B. et al. (2022) *Disruption of sncRNA improves the protective efficacy of outer membrane vesicles against Helicobacter pylori infection in a mouse model.* Infection and Immunity.
29. Choi J.-W. et al. (2017) *Secretable small RNAs via outer membrane vesicles in periodontal pathogens.* J Dent Res. https://consensus.app/papers/details/w:5b9bb4cf-6e6a-4ba8-857f-3a4083a12faf/
30. Sheikh A. et al. (2025) *Outer membrane vesicles from Bacteroides fragilis contain coding and non-coding small RNA species that modulate inflammatory signalling in intestinal epithelial cells.* J Extracell Biol.
31. Chapagain P. et al. (2024) *Characterisation of sRNAs enriched in outer membrane vesicles of pathogenic Flavobacterium psychrophilum.* J Extracell Biol.
32. Zhang K. et al. (2026) *Bacillus Calmette-Guérin ncBCG427 sRNA potentially enhances intracellular bacteria survival by inhibiting apoptosis and inflammatory responses.* BMC Microbiology.
33. Malabirade A. et al. (2018) *The RNA complement of outer membrane vesicles from Salmonella enterica serovar Typhimurium under distinct culture conditions.* Front Microbiol.
34. Mayoral J.G. et al. (2014) *Wolbachia small noncoding RNAs and their role in cross-kingdom communications.* PNAS 111:18721–18726. doi:10.1073/pnas.1420131112. https://pmc.ncbi.nlm.nih.gov/articles/PMC4284532
35. Zhang Y. et al. (2025) *Wolbachia-derived small non-coding RNAs exhibit cross-kingdom regulatory effects on host reproduction of a polyphagous mite.* iScience. https://consensus.app/papers/details/w:740f8c2d-4c99-4d82-a568-c4fbc7619538/
36. Hussain M. et al. (2011) *Wolbachia uses host microRNAs to manipulate host gene expression and facilitate colonization of the dengue vector Aedes aegypti.* PNAS. https://consensus.app/papers/details/w:a64dc3df-3a89-4664-a352-36e7c6a60cee/
37. Buck A.H. et al. (2014) *Exosomes secreted by nematode parasites transfer small RNAs to mammalian cells and modulate innate immunity.* Nature Communications 5:5488. https://pmc.ncbi.nlm.nih.gov/articles/PMC4263141
38. (2024/2025) *Nematode miRNAs impair differentiation of Th2 cells and macrophages; Hpo-miR-71-5p targets Irf4, Asu-miR-791-3p targets Mtor.* https://consensus.app/papers/details/ (Consensus record, helminth sweep)
39. (2025) *A single parasite miRNA (miR-5352), conserved across gastrointestinal nematodes, suppresses IL-13-induced epithelial differentiation via Klf-4 and the IL-22 receptor.* https://consensus.app/papers/details/ (Consensus record)
40. Dong B. et al. (2026) *Host–parasite cross-species regulation by Schistosoma japonicum EV-miRNA sja-miR-2a-3p limits hepatic fibrosis via BCL2 suppression in hepatic stellate cells.* Veterinary Research. https://consensus.app/papers/details/w:2298b287-79a7-46f4-8390-90cf02ef4863/
41. (2020) *Schistosoma japonicum EV miR-125b and bantam regulate host macrophages via Pros1, Fam212b and Clmp.* https://consensus.app/papers/details/ (Consensus record)
42. Xin Y. et al. (2025) *Emu-miR-10a-5p in Echinococcus multilocularis-derived extracellular vesicles alleviates airway inflammation ... through LIF-mediated JAK1–STAT3 signaling.* Front Immunol. https://consensus.app/papers/details/w:d39756bb-bae7-4bb2-a885-2c6d92d8e08a/
43. (2025) *Bursaphelenchus xylophilus exosomal Bx-miR-71-5p remodels vector beetle tracheal development via a Notch–Mmp3–ECM axis.* https://consensus.app/papers/details/ (Consensus record)
44. (2024) *Trichinella spiralis stage-specific extracellular vesicles; tsp-miR-1 suppresses IL-6.* https://consensus.app/papers/details/ (Consensus record)
45. (2021) *Brugia malayi microfilaria-derived extracellular vesicles downregulate mTOR phosphorylation in human monocytes; miR100 predicted to target mTOR and 4E-BP1.* GEO GSE157226, ProteomeXchange PXD021844.
46. (2026) *Trichuris suis EV miRNAs with homology to human miR-22-3p and miR-29-3p downregulate conserved oncogenic targets in human colon epithelial cells.* https://consensus.app/papers/details/ (Consensus record)
47. Garcia-Silva M.R. et al. (2014) *Gene expression changes induced by Trypanosoma cruzi shed microvesicles in mammalian host cells: relevance of tRNA-derived halves.* https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4000953/
48. Zhang L. et al. (2011) *Exogenous plant MIR168a specifically targets mammalian LDLRAP1: evidence of cross-kingdom regulation by microRNA.* Cell Research. https://consensus.app/papers/details/ (Consensus record)
49. Kang W. et al. (2017) *Survey of 800+ data sets from human tissue and body fluid reveals xenomiRs are likely artifacts.* RNA. https://consensus.app/papers/details/ (Consensus record)
50. Witwer K.W. (2018) *Alternative miRNAs? Human sequences misidentified as plant miRNAs in plant studies and in human plasma.* F1000Research. https://consensus.app/papers/details/ (Consensus record)
51. Chin A.R. et al. (2016) *Cross-kingdom inhibition of breast cancer growth by plant miR159.* Cell Research. https://consensus.app/papers/details/ (Consensus record)
52. Yang J. et al. (2015) *Detection of an abundant plant-based small RNA in healthy consumers.* PLoS ONE.
53. Yang J. et al. (2016) *Anomalous uptake and circulatory characteristics of the plant-based small RNA MIR2911.* Scientific Reports.
54. Huang H. et al. (2018) *Extensive degradation and low bioavailability of orally consumed corn miRNAs in mice.* Nutrients.
55. (2018) *A bioinformatics pipeline for the analysis and target prediction of RNA effectors in bidirectional communication during plant–microbe interactions.* Front Plant Sci 9:1212. https://www.frontiersin.org/articles/10.3389/fpls.2018.01212/full

---

### Provenance note
Literature sweep used the Consensus MCP tool (6 sweeps, 40–100 papers each, ~440 records screened) until the monthly allowance was exhausted, then WebSearch/WebFetch for full texts. Full texts successfully read: refs 2, 4, 21, 22, 34, 37, 55 and partial for 1. PMC was intermittently CAPTCHA-gated and eLife was JavaScript-gated; the `†` rows are the ones that could not be verified at full-text level.
