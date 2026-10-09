# Module B – target prediction and precursor folding: design

Status: agreed design, not implemented. The decisions of section 11 were accepted as
proposed (2026-10-09). The general seed-based logic (`ago`) is the default for every target
set; `bacterial` mode runs only on target sets the user declares as bacterial.

## 1. What module B must answer

For each candidate small RNA from module A:

1. Which transcripts could it regulate, through which sites, and how much should we trust
   each prediction?
2. Could it come from a hairpin precursor, as a miRNA would?

The second answer is independent of the first and is reported side by side, not used as
a filter.

The design follows four principles from the plan:

- **Every score is kept.** v0 kept a target only when PITA and RNAhybrid both passed
  fixed cut-offs. Here each predictor's score goes in the table, and consensus is a
  column, not a filter.
- **Trust is measured.** Each score gets an empirical false discovery rate from decoy
  small RNAs that go through the identical pipeline.
- **The site comes first.** In `ago` mode every prediction is anchored to a canonical
  seed site class. Energy predictors score sites; they do not search freely across the
  transcriptome. That keeps the run time bounded and the output interpretable.
- **One general logic, one opt-in exception.** Every target set is checked with the
  general seed-based logic (`ago`: Argonaute-guided small RNAs, seed site classes, UTRs).
  Only target sets the user explicitly declares `mode: bacterial` are searched the
  bacterial way: no seed rule, windows around the start codon, for sRNAs that pair with
  bacterial mRNAs (often through Hfq). Nothing is switched to bacterial mode
  automatically, not even a circular genome.

## 2. What v0 does, and what changes

| step | v0 (smithHunterB.sh) | module B |
|---|---|---|
| target regions | transcripts FASTA + UTR BED, `bedtools getfasta` without `-s` | transcript regions built from a genome + GFF3 (spliced, strand-aware), or transcripts + regions table; 3'UTR by default |
| seed | positions 4–10, reverse complement searched with `agrep`, 0–2 mismatches | positions 2–8 by default; site classes 8mer / 7mer-m8 / 7mer-A1 / 6mer / offset-6mer; 4–10 available as `legacy` |
| extra pairing | `blastn-short`, ≥ 11 matching nt | dropped: IntaRNA and RNAhybrid already score the full duplex |
| accessibility | PITA (column 7 read as ddG, but it is dGduplex) | IntaRNA (ED term, accessibility of both partners) |
| duplex | RNAhybrid `-f 3,10 -e -15 -s 3utr_fly`; seed options not passed | RNAhybrid with the same seed as the scan; p-values calibrated on the actual target regions with RNAcalibrate |
| decision | PITA ≤ −9 and RNAhybrid ≤ −15 | no fixed cut-off: q-value per predictor from decoys, `support` = number of predictors passing |
| precursor | two windows (−15/+50, −50/+15), RNAfold −T 25, MFE only; 1-based / 0-based mix-up | same two windows, configurable; MFE, AMFE, MFEI, shuffle p-value, hairpin check; correct coordinates, circular genomes |

## 3. Inputs and configuration

Module B starts from `discovery/candidates.fasta` and `discovery/loci.tsv`. As in v0,
the user may review candidates between the two modules. The `candidates` setting points
to any FASTA or TSV of candidates, so an edited list can replace the default.

Targets are declared as **target sets**. A set is a group of transcripts with one
interaction mode; `mode` is optional and defaults to `ago`. Each focus genome lists the
sets its small RNAs are tested against. This replaces the single per-genome
`targets: ago | bacterial | none` field of module A: the mechanism belongs to the target
cell, not to the genome that makes the small RNA. The old field stays accepted and is read
as "all target sets of that mode".

```yaml
genomes:
  mito:
    fasta: data/mt.fasta
    role: focus
    topology: circular
    target_sets: [host]           # small RNAs of mito are tested against host transcripts
  symbiont:
    fasta: data/symbiont.fasta
    role: focus
    target_sets: [symbiont_mrna]

target_sets:
  host:                                      # no mode: the general logic (ago)
    # either a genome + GFF3 ...
    genome: data/nuclear.fasta
    annotation: data/nuclear.gff3
    # ... or transcripts + an optional table of regions in transcript coordinates
    # transcripts: data/transcripts.fasta
    # region_table: data/regions.tsv         # transcript, region, start, end (1-based)
    regions: [three_prime_UTR]               # add five_prime_UTR, CDS if wanted
    whole_transcript_if_no_utr: false        # many non-model transcriptomes lack UTRs
  symbiont_mrna:
    mode: bacterial                          # only where the user asks for it
    genome: data/symbiont.fasta
    annotation: data/symbiont.gff3
    window: [-200, 100]                      # around each start codon, as in CopraRNA

targets:
  candidates: results/discovery/candidates.fasta
  seed: [2, 8]                 # guide positions; legacy: [4, 10]
  site_classes: [8mer, 7mer-m8, 7mer-A1, 6mer, offset-6mer]
  predictors: [intarna, rnahybrid]          # + miranda
  decoys: 20                   # shuffled decoys per candidate
  fdr: 0.05

precursors:
  flanks: [15, 50]             # windows: [start-15, end+50] and [start-50, end+15]
  temperature: 25              # required, no default (v0 used 25 for bivalves)
  shuffles: 1000
```

### Building target regions

- **Genome + GFF3** (preferred). Exons of each mRNA are joined in transcript order and
  reverse-complemented on the minus strand. 3'UTR and 5'UTR come from `three_prime_UTR` /
  `five_prime_UTR` features. When these are missing they are derived as exon parts after
  or before the CDS. Each region keeps a map back to genomic coordinates, so sites can be
  reported in both systems. Pure Python, in the style of `annotate.py`.
- **Transcripts FASTA + regions table**: for transcriptome-only organisms, as in v0. Its
  coordinates are documented as 1-based inclusive and converted once. This avoids the
  v0 BED/seqret off-by-one.
- **Bacterial windows**: −200..+100 nt around each annotated start codon, on the CDS strand,
  allowing for circular contigs.
- **Redundancy**: identical region sequences (isoforms sharing a UTR) are scanned once and
  reported for every transcript. Sites, not transcripts, are the unit for FDR, so isoforms
  do not inflate the counts.

## 4. General logic (`ago` mode, the default)

### 4.1 Site scan (pure Python, exact, fast)

Guide positions are numbered from the 5' end of the small RNA; target position t1 is
opposite guide position 1. With the seed at 2–8 (McGeary et al. 2019; Bartel 2018):

| class | requirement on the target |
|---|---|
| 8mer | match to guide 2–8 + A at t1 |
| 7mer-m8 | match to guide 2–8 |
| 7mer-A1 | match to guide 2–7 + A at t1 |
| 6mer | match to guide 2–7 |
| offset-6mer | match to guide 3–8 |

Each site is assigned only its best class, so an 8mer is not also counted as a 7mer.

- **Other seed ranges** (`legacy: [4, 10]`) are reported as `seed<start>-<end>`, since
  canonical classes are defined only for 2–8.
- **Mismatches and G:U pairs in the seed** are off by default. They can be enabled as
  separate classes (`6mer-GU`, `seed-1mm`), so they never mix with canonical sites.

Each site gets context features that cost nothing to compute and are known to modulate
efficacy (Grimson et al. 2007):

- distance from the region ends (sites < 15 nt after the stop codon are flagged);
- local AU content (±30 nt);
- 3' supplementary pairing, the best match to guide 13–16 within a short offset;
- number and spacing of sites for the same small RNA in the same transcript.

### 4.2 Scoring sites

Only regions that contain at least one site go to the energy predictors.

- **IntaRNA**:
  - Input: the full region, so accessibility is computed in context (`--tAccW`/`--tAccL`).
  - Interactions are constrained to the site with `--tRegion`, and the seed is pinned to
    the guide seed range (`--seedQRange`, `--seedBP`).
  - Kept per site: total energy E, hybridisation energy and the accessibility penalties ED
    of both partners.
  - The query is short, so its accessibility is computed but expected to matter little.
- **RNAhybrid**:
  - Input: a window of site −15/+(guide length + 15) nt.
  - The helix constraint matches the seed range (`-f`), which fixes the v0 bug where seed
    options were not passed.
  - p-values use extreme-value parameters from `RNAcalibrate` on the same target regions
    (`-d xi,theta`) instead of the fly/worm/human presets. That keeps them meaningful for
    non-model taxa.
  - Kept: MFE, p-value, duplex string.
- **miRanda** (optional): score and energy, run on the same windows.

### 4.3 Decoys and FDR

- **Generating decoys.** Each candidate gets `decoys` shuffled sequences that preserve
  dinucleotide composition (Altschul–Erickson; uShuffle, Jiang et al. 2008), in pure
  Python with a fixed seed. Decoys whose seed equals the seed of any real candidate are
  redrawn.
- **Running them.** Decoys go through the identical scan and scoring.
- **Computing q-values.** For each predictor and each site class, real and decoy scores
  are pooled. The FDR at a threshold s is (decoy sites scoring ≥ s / number of decoys
  per candidate) / (real sites scoring ≥ s). The q-value is the minimum FDR over
  thresholds at or below the site's score.
  - Per class, because 8mers and 6mers have very different background frequencies.
  - Pooled over candidates, for stable estimates, with a per-candidate version reported
    alongside.
- **Combining predictors.** `support` is the number of predictors with q ≤ `fdr`. The
  default ranking is by class (8mer first), then by support, then by IntaRNA E.

The FDR measures the specificity of the whole procedure, not the biology. A site that
beats the decoys is a real complementarity signal, not proof of regulation. The report
says so.

### 4.4 Cost

The scan is linear in the size of the target regions. A 3'UTRome of ~10 Mb gives roughly
10³–10⁴ sites per candidate, mostly 6mers. With 30 candidates × 21 sequences (1 + 20
decoys) that makes about 10⁵–10⁶ IntaRNA site evaluations: hours on a few cores, not days.

- **Parallelism.** Snakemake scatters by batches of candidates (real and decoy together).
- **Defaults.** 6mer and offset-6mer go to the energy predictors only with
  `score_weak_sites: true`. Otherwise they are listed with class and context only. This
  is the main knob if the run is too long.
- **Testing.** Before any heavy run, a dry run reports the number of sites per class and
  the expected run time.

## 5. `bacterial` mode (opt-in, per target set)

Used only for target sets declared with `mode: bacterial`. Implemented after the general
logic is working and tested.


- **Regions**: −200..+100 nt around each start codon.
- **Search**: IntaRNA over the whole window. There is no seed range, but IntaRNA's own
  seed requirement applies (default 7 bp, any position). Both partners' accessibility is
  included.
- **Predictors**: RNAhybrid and the site classes are not used (they model AGO guides).
  CopraRNA, which needs several related genomes, is a later option.
- **FDR**: decoys as in 4.3, pooled over all windows.
- **Caveat**: a small RNA library under ~40 nt sees fragments of bacterial sRNAs that are
  50–500 nt long. The query is the representative read by default. With `query: locus`,
  the whole locus span (or a user-supplied FASTA of full-length sRNAs) is used instead.

## 6. Precursor folding

For each candidate locus, two windows on the locus strand cover the candidate in either
arm of a hairpin: `[start − long, end + short]` and `[start − short, end + long]`, with
`flanks: [short, long]` = [15, 50] as in v0. Windows can cross the origin of circular
genomes.

Per window:

- **Folding**: RNAfold MFE structure at the configured temperature.
- **Energy measures**:
  - AMFE = MFE / length × 100;
  - MFEI = AMFE / GC%. Real pre-miRNAs are typically above ~0.85 (Zhang et al. 2006); this
    is used as a descriptor, not a cut-off.
- **Significance**: the fraction of `shuffles` dinucleotide-shuffled windows with MFE ≤ the
  real one, as in randfold (Bonnet et al. 2004).
- **Hairpin check**: whether the structure is a single dominant stem-loop, whether the
  candidate lies in one arm, and the fraction of its bases paired.
- **Output**: dot-bracket string and, optionally, an SVG from RNAplot.

The best of the two windows is reported, and both are kept in the full table.

## 7. Outputs

| file | content |
|---|---|
| `targets/sites.tsv` | one row per candidate × region × site: coordinates (transcript and genomic), class, context features, every predictor's scores, q-values, support |
| `targets/targets.tsv` | one row per candidate × transcript: sites per class, best scores, best q-value, gene name |
| `targets/fdr_curves.tsv` | real vs decoy counts per predictor, class and threshold, for the report |
| `targets/decoys.fasta` | the decoys used, for reproducibility |
| `precursors/precursors.tsv` | both windows per candidate: sequence, structure, MFE, AMFE, MFEI, p-value, hairpin flags |
| `candidates_summary.tsv` | one row per candidate joining module A (expression, ends, isolation, annotation) and module B (sites by class, targets at q ≤ fdr, precursor) |

## 8. Code layout

Same approach as module A: a Python core with no third-party dependencies, external tools
called only in Snakemake rules, and their output parsed in Python.

```
src/smithhunter/targets/
  regions.py     GFF3 / transcripts -> regions, with genomic coordinate maps
  sites.py       seed site scan, classes, context features
  shuffle.py     dinucleotide-preserving shuffle (also used for precursors)
  intarna.py     command builder + CSV parser
  rnahybrid.py   command builder + parser, RNAcalibrate parameters
  miranda.py     optional
  fdr.py         decoy-based q-values
src/smithhunter/precursor.py   windows, RNAfold parsing, MFEI, p-value, hairpin check
workflow/rules/targets.smk, precursors.smk
workflow/envs/intarna.yaml, rnahybrid.yaml, viennarna.yaml, miranda.yaml
```

New subcommands: `regions`, `sites`, `decoys`, `score-sites`, `fdr`, `precursors`.
Module B runs from the same Snakefile with a target `targets` (`snakemake ... targets`),
so module A can be run alone and reviewed first.

## 9. Testing

- **Unit tests, without external tools**:
  - site classes on hand-made targets, including an 8mer that must not also count as a 7mer;
  - minus-strand and spliced region extraction;
  - GFF without UTR features;
  - the shuffle preserves dinucleotide counts exactly;
  - the q-value computation;
  - window coordinates across the origin.
- **Parsers**: tested on small stored outputs of IntaRNA, RNAhybrid and RNAfold.
- **Simulation**, extending `smithhunter simulate`:
  - synthetic nuclear transcripts with 3'UTRs;
  - planted sites of known class for some candidates, and none for others;
  - planted hairpins around some loci.

  `evaluate` then measures the recovery of planted sites per class and whether the q-values
  are calibrated: unplanted candidates should give about `fdr` of their calls.
- **CI**: the simulation run, with a small UTRome and few decoys.

## 10. Implementation order

Progress: step 1 is implemented (`regions`, `decoys`, `sites` subcommands, rule
`targets`). On the simulation all 15 planted sites are found with their class. As
expected, the site-count FDR is close to 1: a seed match is about as frequent for a real
small RNA as for its decoys. Counting sites shows how much of the signal is background
but does not rank targets. Site-level q-values need the scores of step 2.


1. Regions, site scan, decoys, FDR on site counts only (no external tools). This is
   already a usable result: seed sites by class with an empirical FDR.
2. RNAhybrid with calibration, then IntaRNA. Consensus and the tables.
3. Precursor folding.
4. Bacterial mode.
5. Simulation of targets and precursors; CI.
6. miRanda, CopraRNA, report.

## 11. Decisions (accepted as proposed, 2026-10-09)

1. **Default seed**: 2–8 (canonical, Plazzi et al. 2024 for smithRNAs), with 4–10 as
   `legacy`? Proposed: 2–8.
2. **Config schema**: move the interaction mode from the small RNA genome (`targets: ago`)
   to target sets (`target_sets` + `mode`)? Proposed: yes; the old field would be read as
   "all sets of that mode".
3. **Default target regions**: 3'UTR only, or 3'UTR + CDS? v0 used the UTRs in the BED
   file. Proposed: 3'UTR only, with CDS and 5'UTR opt-in.
4. **Folding temperature default**: 37 °C (ViennaRNA default), or required in the config
   with no default, since v0 used 25 °C for bivalves? Proposed: required, so it is never
   chosen silently.
5. **Number of decoys**: 20 per candidate by default (cost ×21)? Fewer makes the FDR
   coarse at strict thresholds.
6. **Weak sites** (6mer, offset-6mer) scored by the energy predictors by default, or
   listed only? Proposed: listed only, scoring opt-in.

## References

- Bartel DP. Metazoan microRNAs. Cell 173:20 (2018).
- McGeary SE et al. The biochemical basis of microRNA targeting efficacy. Science 366:eaav1741 (2019).
- Grimson A et al. MicroRNA targeting specificity in mammals: determinants beyond seed pairing. Mol Cell 27:91 (2007).
- Plazzi F et al. Mitochondrially mediated RNA interference. Heredity 132:156 (2024).
- Marturano G, Carli D et al. SmithHunter. BMC Bioinformatics 25:286 (2024).
- Mann M, Wright PR, Backofen R. IntaRNA 2.0. Nucleic Acids Res 45:W435 (2017); Raden M et al. 2020 (seed and accessibility constraints).
- Krüger J, Rehmsmeier M. RNAhybrid. Nucleic Acids Res 34:W451 (2006).
- Wright PR et al. CopraRNA. Nucleic Acids Res 41:W119 (2013).
- Altschul SF, Erickson BW. Mol Biol Evol 2:526 (1985); Jiang M et al. uShuffle. BMC Bioinformatics 9:192 (2008).
- Bonnet E et al. Evidence that microRNA precursors have lower MFE than shuffled sequences. Bioinformatics 20:2911 (2004).
- Zhang BH et al. Evidence that miRNAs are different from other RNAs. Cell Mol Life Sci 63:246 (2006).
- Lorenz R et al. ViennaRNA Package 2.0. Algorithms Mol Biol 6:26 (2011).
