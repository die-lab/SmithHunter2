# SmithHunter2

**Find the small RNAs produced by every genome in your sample.**

Small RNA-seq libraries mix reads from several genomes: the nuclear and mitochondrial
genomes of the organism, plastids, symbionts, parasites, contaminants. SmithHunter2
lets you declare all of them, choose which ones you want to study, and finds small RNA
loci on those genomes while the others compete for the reads.

- Every read is assigned to the genome it fits **best**, not to the first genome it
  happens to fit. A read from the mitochondrion is not lost because a NUMT in the
  nuclear genome resembles it, and nuclear reads do not inflate mitochondrial loci.
- Reads that fit two genomes equally well are reported as **ambiguous** and handled by
  a policy you choose. A stricter genome choice is available (`assignment.mode`: `best`,
  `severe`, `evil`, see below).
- Circular genomes are handled natively: loci across the origin are found as one locus.
- Loci are called from the 5' ends of reads, where small RNAs pile up and degradation
  fragments do not, so small RNAs inside a degraded rRNA or mRNA keep their own locus.
  They are quantified per sample, scored for the precision of their 5' and 3' ends, and
  annotated with your GFF3.
- A built-in **simulator** writes a realistic library with known truth, so you can see
  what the pipeline recovers before trusting it on real data.

SmithHunter2 is a rewrite of [SmithHunter](https://github.com/ESZlab/SmithHunter)
(Marturano, Carli et al., *BMC Bioinformatics* 25:286, 2024), written to find smithRNAs
(small mitochondrial highly transcribed RNAs) but no longer limited to mitochondria.

> **Status: early development.** Module A (small RNA discovery) works and is tested.
> Module B (target prediction) is in progress: seed sites scored with IntaRNA and RNAhybrid,
> with q-values from decoys, are available; precursor folding is next
> ([design](docs/module_b_design.md)).

## How it works

```
FASTQ (SE or PE)
  │ fastp: adapter and quality trimming (PE: R2 only corrects R1 by overlap)
  ▼
collapse identical reads across samples ── each unique sequence is aligned once
  │
  ▼
bowtie 1 on ONE reference holding all genomes (-v 1, --best --strata, up to 50 hits)
  │
  ▼
origin of each read: assigned │ excluded │ ambiguous │ too_many_hits │ unmapped
  │
  ▼  focus genomes only
blocks of overlapping alignments on the same strand
  │ split at 5' end peaks: small RNAs vs degradation background
  ▼
loci: counts and RPM per sample, expression filter above local background
  │ 5' / 3' end scores (after SmithHunter's sharp_smith)
  │ annotation with GFF3
  ▼
loci.tsv + candidates.fasta
```

Each genome gets a **role**:

| role | what happens to its reads |
|---|---|
| `focus` | loci are defined, quantified, scored and annotated |
| `exclude` | the genome competes for reads, which are counted and then set aside |

A read is `assigned` when all its best alignments fall in one focus genome,
`excluded` when they fall only in exclude genomes, and `ambiguous` when they fall in a
focus genome and in another genome. The `ambiguous` setting decides what happens to
ambiguous reads: `report` (listed, not used for loci), `discard`, or `split` (shared
among all best alignments). Within one genome, reads with several equally good
alignments are shared among them (`multimap: fractional`) or dropped (`unique`).

#### Genome choice: `assignment.mode`

For a read that fits a focus genome and another genome, mismatches are compared genome by
genome (best alignment of the read in each genome):

| mode | rule |
|---|---|
| `best` (default) | the genome with fewer mismatches gets the read; a tie is a double-mapping read (`ambiguous` policy) |
| `severe` | as `best`, but a tie is always discarded |
| `evil` | as `severe`, and the read is also discarded (`evil_rejected`) when an exclude genome fits it at most `assignment.evil_max_diff` mismatches worse than the focus genome, even though the focus genome fits better |

`evil_max_diff` must not exceed `mapping.mismatches`; in `evil` mode bowtie also reports the
worse alignments (no `--strata`/`-m`, up to `max_hits + 1` alignments) so they can be seen.

### How loci are called

Reads that overlap on one strand form a block. Inside a block, 5' end positions are
taken from the most to the least abundant: a position, with the reads starting up to
`peak_window` bases away, becomes a **peak locus** when its read count is significantly
above the local background (Poisson test). The background is the mean number of reads
starting at each position within `background_flank` bases, leaving out stronger peaks.
Reads claimed by no peak (typically degradation fragments) form **background loci**,
which are reported but never pass. A weak peak within `satellite_distance` bases of a
peak with ten times more reads (`satellite_fraction`) is taken as a 5' variant of it and
joins its locus. The **isolation** of a locus is the share of all reads starting within
`satellite_distance` bases of its peak that belong to it: regions with scattered 5'
ends, or minor neighbours of a stronger peak, have low isolation and do not pass. A sample supports a peak locus when the locus is
significantly above that sample's own background and, after subtracting it, reaches
`min_count` and `min_rpm`. With `peaks: false`, every block is one locus.

## Installation

Requirements: Linux and conda (or mamba). No root access is needed.

```bash
git clone https://github.com/die-lab/SmithHunter2.git
cd SmithHunter2
mamba env create -f environment.yml      # fastp, bowtie 1, snakemake, python, IntaRNA, RNAhybrid
mamba activate smithhunter2
```

The Python code has no third-party dependencies. Instead of the single environment,
Snakemake can build one environment per rule: add `--software-deployment-method conda`
to the commands below.

## Quick start

1. Put your FASTQ files and genome FASTA files in a working directory.
2. Copy `config/config.yaml` and `config/samples.tsv` into it and edit them.
3. From the working directory:

```bash
snakemake -s /path/to/SmithHunter2/workflow/Snakefile --cores 8
```

`config/samples.tsv` has one row per sample; leave `fq2` empty for single-end data:

```
sample	fq1	fq2
rep1	data/rep1_1.fastq.gz	data/rep1_2.fastq.gz
rep2	data/rep2_1.fastq.gz	data/rep2_2.fastq.gz
```

## Configuration

The main block declares the genomes:

```yaml
genomes:
  mito:
    fasta: data/organism_mt.fasta
    role: focus
    topology: circular                 # linear (default) or circular
    annotation: data/organism_mt.gff3  # optional, e.g. from MITOS2
  nuclear:
    fasta: data/organism_nuc.fasta
    role: exclude
  symbiont:
    fasta: data/symbiont.fasta
    role: focus
    topology: circular
    length: [18, 40]                   # read lengths for this genome
```

Other settings, with their defaults:

| setting | default | meaning |
|---|---|---|
| `ambiguous` | `report` | `report`, `discard` or `split` (see above) |
| `assignment.mode` / `evil_max_diff` | `best` / 0 | genome choice: `best`, `severe` or `evil`, and the tolerance of `evil` (see above) |
| `trimming.adapter_r1` / `adapter_r2` | TruSeq small RNA | 3' adapters; empty for fastp autodetection |
| `trimming.min_length` / `max_length` | 18 / 35 | read length range kept |
| `mapping.mismatches` | 1 | bowtie `-v` |
| `mapping.max_hits` | 50 | reads with more equally good alignments go to `too_many_hits` |
| `loci.merge_gap` | 0 | join alignments closer than this many bases into one block |
| `loci.peaks` | `true` | split blocks at 5' end peaks |
| `loci.peak_window` | 2 | reads starting this close to a peak belong to it |
| `loci.peak_pvalue` | 0.001 | Poisson test of a peak against the local background |
| `loci.background_flank` | 50 | bases on each side used to estimate the background |
| `loci.satellite_distance` / `satellite_fraction` | 10 / 0.1 | a weaker peak this close, with less than this share of reads, joins the stronger one |
| `loci.multimap` | `fractional` | `fractional` or `unique` within one genome |
| `loci.min_rpm` / `min_count` | 5 / 5 | per-sample thresholds (reads per million and raw reads) |
| `loci.min_samples` | n − 1 | samples that must reach both thresholds |
| `ends.five_filter` | `dominant` | how the 5' end is tested: `dominant`, `jitter` or `score` (below) |
| `ends.min_five_fraction` | 0.5 | `dominant` / `jitter`: a locus passes if this share of its reads starts at (near) the top 5' position |
| `ends.five_jitter` | 1 | `jitter`: bases on each side of the top 5' position that still count |
| `ends.min_five_score` | 0.5 | `score`: a locus passes if its 5' score is above this |
| `ends.min_three_score` | 0.0 | and its 3' score is at least this |
| `ends.min_isolation` | 0.5 | and its isolation is at least this |

The end score of a locus is 1 when one position holds at least half of the reads
(`n_thre`), 0.5 for two adjacent positions, and lower for scattered ends. It only takes
a few values (1, 0.5, 0.476 for two ends 2 nt apart, 0.333, ...), so the 5' test uses the
read shares directly by default:

| `five_filter` | the 5' end passes when |
|---|---|
| `dominant` (default) | reads starting at the top 5' position / reads of the locus ≥ `min_five_fraction` |
| `jitter` | reads starting within `five_jitter` nt of the top position / reads of the locus ≥ `min_five_fraction` |
| `score` | `five_score` > `min_five_score` (the earlier rule; with the defaults it passes the same loci as `dominant`) |

Example (a real locus, 102 reads, 5' ends at +1: 16, +2: 39, +3: 5, +4: 32, 10 elsewhere):
the dominant fraction is 39/102 = 0.38 and fails at 0.5; with `jitter` the share within
1 nt of +2 is (16 + 39 + 5)/102 = 0.59 and passes; within 2 nt it would be 0.91.

## Outputs

All paths are under `results/` (the `outdir` setting).

| file | content |
|---|---|
| `discovery/loci.tsv` | every locus, with filter flags (`pass`) |
| `discovery/candidates.fasta` | representative sequence of each locus that passes |
| `discovery/locus_reads.tsv.gz` | the reads of every locus, for manual review |
| `origin/origin_summary.tsv` | reads per origin class and genome, per sample |
| `origin/length_distribution.tsv` | the same, by read length |
| `origin/read_origin.tsv.gz` | origin class of every unique read |
| `qc/*.fastp.html` | trimming reports |
| `report/read_funnel.tsv` | reads left after each step, per sample |
| `report/provenance.yaml` | config, software and conda pins, git commit and locus counts of the run |

Main columns of `loci.tsv`:

| column | meaning |
|---|---|
| `start`, `end`, `strand`, `span` | 1-based locus coordinates; `wraps_origin` = 1 if it crosses the origin |
| `locus_type` | `peak` (candidate small RNA) or `background` (reads with no 5' peak; never passes) |
| `rep_sequence`, `rep_start`, `rep_end` | the most abundant read of the locus and its position |
| `total_count`, `count_<sample>` | reads (multi-mapping reads count fractionally) |
| `background_count` | reads expected from the local background in the peak window |
| `rpm_<sample>` | reads per million of the trimmed library |
| `rpm_genome_<sample>` | reads per million of the reads assigned to that genome |
| `five_score`, `three_score` | end precision (see above); `*_dominant_fraction` = share of the top position; `five_jitter_fraction` = share within `five_jitter` nt of it |
| `isolation` | share of the reads starting near the peak that belong to the locus |
| `annotation_class`, `annotation_orientation` | most specific overlapping GFF3 feature and its strand relation |
| `expression_pass`, `ends_pass`, `pass` | filter results |

### Plots and QC report

The `report` rule adds plots (R/ggplot2) and a MultiQC report. It needs R and MultiQC,
so it is not part of the default run; with per-rule conda environments:

```bash
snakemake -s /path/to/SmithHunter2/workflow/Snakefile --cores 8 --sdm conda report
```

| file (`results/report/`) | content |
|---|---|
| `reads.pdf` | read funnel (all samples and per sample), genome of origin, read length by origin |
| `genomes.pdf` | coverage of each focus genome contig by strand, passing loci labelled by number |
| `loci.pdf` | profile of each passing locus: depth, 5' and 3' ends per position |
| `multiqc_report.html` | fastp and bowtie reports, read funnel and origin per sample |
| `coverage.tsv.gz` | depth and 5' ends per position of focus genomes (input of `genomes.pdf`) |

The reporting layer (funnel, provenance, ggplot2 plots, MultiQC) follows the one of the
smithhunter24 Snakemake port of SmithHunter v0 module A.

## Target prediction (module B, in progress)

After reviewing the candidates of module A, declare the transcripts to search as
**target sets** and run the `targets` rule:

```yaml
target_sets:
  host:
    genome: data/organism_nuc.fasta
    annotation: data/organism_nuc.gff3
    regions: [three_prime_UTR]       # also five_prime_UTR, CDS, transcript
targets:
  candidates: results/discovery/candidates.fasta
  seed: [2, 8]
  decoys: 20
```

```bash
snakemake -s /path/to/SmithHunter2/workflow/Snakefile --cores 8 targets
```

Every target set uses the general seed-based logic: canonical site classes (8mer,
7mer-m8, 7mer-A1, 6mer, offset-6mer) of the seed at guide positions 2–8. Bacterial mode
(no seed rule, windows around start codons) will be opt-in per target set
(`mode: bacterial`) and is not implemented yet. Target regions come from a genome with
its GFF3 (spliced, both strands; UTRs are derived from the CDS when the GFF has no UTR
features) or from a transcripts FASTA.

| file (`results/targets/<set>/`) | content |
|---|---|
| `sites.tsv` | every site: candidate, transcript, gene, class, position in the region and in the genome, distance from the region ends, AU content of the flanks, 3' supplementary pairing (guide 13–16) |
| `targets.tsv` | one row per candidate and transcript: sites per class, best class |
| `site_fdr.tsv` | sites per class for each candidate and for its dinucleotide-shuffled decoys |

The decoys show how many sites a sequence of the same composition finds by chance.
Expect this to be close to the real count: a seed match on its own does not make a
target. The next step scores each site.

### Scoring sites (IntaRNA and RNAhybrid)

The `targets` rule also runs `score-sites` on every site of the candidates **and of their
decoys**, with the same code and settings:

- **IntaRNA** (3.4): the interaction is anchored on the site. The target is a window of the
  region (the site plus 150 nt each side, so accessibility is computed in context), the
  seed is pinned to the guide seed range and to the seed match on the target
  (`--seedQRange`, `--seedTRange`) and the interaction is confined to the guide's reach
  (`--tRegion`). Kept: total energy `E` (hybridisation plus accessibility cost), the
  accessibility penalties `ED` of target and query, the unpaired probabilities and the duplex.
- **RNAhybrid**: a window of the site plus 15 nt, with the helix constraint (`-f`) set to the
  seed range of the site class. p-values use the extreme-value parameters of
  `RNAcalibrate` (`-d xi,theta`), fitted for each real small RNA on random sequences with
  the dinucleotide frequencies of the target set, with the same helix constraint, as the RNAhybrid
  paper advises (Rehmsmeier et al. 2004). The sequences have RNAcalibrate's default lengths
  (500 +- 300 nt): fitting at the real window length (~55 nt) was tried and fails (negative xi),
  because the method was designed for long targets. RNAhybrid then applies the parameters to the
  short window through its length normalisation, which is an extrapolation. The p-value is the
  one of the best constrained duplex in such a window, conditional on the site, not a
  transcriptome-wide search p-value.
  Without calibration the `3utr_human` preset is used and the p-value is only indicative.
  The q-values use the MFE and the energy, never this p-value.
- **q-values**: for each predictor and site class, the FDR at a score is
  (decoy sites at least as good, per candidate) / (real sites at least as good), and the
  q-value of a site is the lowest FDR over the thresholds it passes. Pooled over candidates
  (`*_q`) and for the candidate alone (`*_q_cand`). `support` is the number of predictors
  with q <= `targets.fdr`. No score is a cut-off: all values stay in the table.
- **Decoys** (`decoy_mode`): `shuffle` (default) uses dinucleotide shuffles of the whole small
  RNA, so a decoy has another seed and the test asks whether this sequence is special;
  `keep_seed` shuffles only the part 3' of the seed, so every decoy has exactly the real seed
  sites and the test isolates the pairing beyond the seed (it controls for seed GC content).
  With `shuffle`, a GC-rich seed tends to get low q-values at all its sites and an AU-rich seed
  high ones, which is worth keeping in mind when reading per-candidate q-values (`*_q_cand`).
  A q-value of 0 means better than every decoy site; its resolution is about one decoy site
  divided by the number of decoys.
- Candidates without decoys (e.g. added to a reviewed candidate list after the decoys step)
  get no q-values and never `support`; pooled q-values weight each decoy site by
  1 / the number of decoys of its candidate. Without any decoy the run stops.
- RNAhybrid duplexes must cover the seed of the site (`rnahybrid_at_site`); otherwise the best
  helix of the window lies on another seed match and the site has no RNAhybrid score (real and
  decoy sites alike; `--rnahybrid-any-duplex` keeps them). Every tool call has a timeout
  (300 s); a failing site is left unscored and the run stops only if more than 5% fail.
- **Two-stage mode** (`intarna_prefilter: 0.15`): RNAhybrid costs about 0.01 s per site and
  IntaRNA about 0.05-0.15 s, and a human 3'UTRome gives ~5,000 strong sites per small RNA, times
  21 sequences with decoys. With a prefilter, RNAhybrid scores every site and IntaRNA only the
  best share of each class by MFE. The cutoff is one pooled threshold per class over real and
  decoy sites, so both meet the same rule and the q-values remain valid (conditional on passing
  it); sites below the cutoff have RNAhybrid values only. Default 1.0 (off).
- 6mer and offset-6mer sites are only listed (`sites.tsv`); scoring them is opt-in
  (`score_weak_sites: true`) because they are ten times more frequent.

```yaml
targets:
  predictors: [intarna, rnahybrid]   # [] = seed sites only
  score_weak_sites: false
  temperature: 37                    # deg C, set it for your organism
  calibrate_samples: 1000            # 0 = RNAhybrid preset
  intarna_prefilter: 1.0             # e.g. 0.15 for large target sets
  decoy_mode: shuffle                # or keep_seed
  fdr: 0.05
```

`smithhunter score-sites --dry-run ...` counts the sites to score before committing to a
long run (about 0.01 s per site for RNAhybrid, 0.05-0.15 s for IntaRNA). Long runs are
**resumable**: the scores are saved in the work directory (`checkpoint_stage1.pkl` after the
RNAhybrid stage, `checkpoint_scores.pkl` at the end) and reused when the small RNAs, decoys,
target sequences and settings are the same (`--no-resume` ignores them). A killed run restarts
from its last stage, and q-values can be recomputed without scoring again. The checkpoints
are Python pickles written by the program: do not point the work directory at files from
elsewhere. Progress and an ETA are printed every 5%.

| file (`results/targets/<set>/`) | content |
|---|---|
| `sites_scored.tsv` | every scored site: class, context, IntaRNA and RNAhybrid scores, q-values, `support` |
| `targets_scored.tsv` | one row per candidate and transcript: best energies and q-values, `support` |
| `fdr_curves.tsv` | real vs decoy sites per predictor, class and threshold |

#### Statistical notes and limits

- **What the FDR is.** An empirical, separate-decoy target-decoy estimate (as in Elias & Gygi
  2007, Nat Methods), with the q-value as the minimum FDR over the thresholds a site passes
  (Storey & Tibshirani 2003, PNAS) and no pi0 correction (which makes it conservative). Such
  estimates are conservative on average but do not provably control the FDR (He et al. 2015;
  Gupta et al. 2011, JASMS; Kertesz-Farkas et al. 2015, J Proteome Res).
- **The null.** A site is "good" if it scores better than the sites of composition-matched random
  small RNAs, not if it is functional. Functional targets also depend on conservation,
  accessibility and the cell, which no shuffle models; energy alone predicts repression weakly
  (Agarwal et al. 2015, eLife). Sites within one target or one small RNA are not independent.
- **The +1.** One pseudo decoy site (weight 1 / decoys of the candidate) is added to every FDR
  numerator, the (D+1)/T estimate (Barber & Candes 2015; Phipson & Smyth 2010), so q is never 0:
  its floor is (1 / number of decoys) / number of real sites passing. With thousands of sites
  per class this changes the pooled q by ~1e-5; for a candidate with few sites it matters (with
  20 decoys a lone site cannot go below q = 0.05). `fdr_pseudocount: false` removes it.
- **Per-candidate q-values** (`*_q_cand`) rest on one small RNA's own sites and its own decoys:
  low resolution (smallest step 1 / decoys), noisier than the pooled ones. Use them as a
  second opinion, not a threshold.
- **Pooled q-values** assume real false sites and decoy sites share a score distribution across
  candidates. They need not when candidates differ in length, GC content or seed; the decoy
  mode matters here (below). Classes are kept apart because 8mers and 7mers have different
  backgrounds; small cells give unstable estimates.
- **Decoy modes** (`decoy_mode`). `shuffle`: whole-sequence dinucleotide shuffles (Altschul &
  Erickson 1985; Clote et al. 2005); their seeds differ from the real one, so seed GC is a
  confounder. `match_seed_gc`: as `shuffle`, with the seed's G+C count matched. `keep_seed`: only
  the part 3' of the seed is shuffled, so decoys have the real seed sites; the null is then
  "pairing beyond the seed", a different question, not an invalid one. On the *Pseudomonas*
  panel `shuffle` and `keep_seed` gave the same supported candidate (see the calibration notes).
- **`support`** is the number of predictors (0-2) with pooled q <= the FDR. IntaRNA and RNAhybrid
  both rank duplex energy and are correlated (Spearman 0.3-0.4 in our runs), so two votes are
  not two independent tests and the count is not a probability. The error rate of
  "support >= 1" can be as high as the sum of the two thresholds.
- **Two-stage mode** applies one pooled, label-blind cutoff per class to real and decoy sites
  alike, so the q-values of IntaRNA are conditional on passing the RNAhybrid prefilter (the
  decoy count and the real count are both taken among the survivors). A real site with a poor
  RNAhybrid MFE and a good IntaRNA energy (accessibility-driven) is lost, and decoys cannot
  measure that loss. On the simulation, prefilters of 30% and 15% recovered all 4 supported
  sites of the one-stage run; this is a small test, not a recall estimate for real data.
- **IntaRNA has not been benchmarked on miRNA data** in the papers we found (Mann et al. 2017,
  NAR; Raden et al. 2020, BMC Bioinformatics: bacterial sRNAs only, and seed constraints
  "need careful tuning" because they cause false negatives). Pinning the seed to the canonical
  site restricts where the seed may sit; `--seedBP` is set per class (7 for 8mer and 7mer-m8,
  6 for 7mer-A1 and 6mer). The accessibility window is the IntaRNA default scale (150 / 100
  nt); a flank of 100 nt around the site was used. The result's sensitivity to the flank has
  not been measured.

`smithhunter validate-targets --targets targets_scored.tsv --truth pairs.tsv` ranks
experimentally validated small RNA -> gene pairs among all genes (rank and percentile by
IntaRNA, RNAhybrid and `support`, AUC over the pairs). With a handful of validated pairs
this is a sanity check, not an estimate of precision.

## Benchmark on simulated data

`smithhunter simulate` writes a complete experiment with known truth: synthetic
mitochondrial (focus, circular, with a GFF3), nuclear and bacterial (exclude) genomes,
FASTQ files, a ready-made configuration, and the list of every simulated small RNA with
the outcome expected from the pipeline. The read model follows the small RNA-seq
literature: precise 5' ends and variable 3' ends, non-templated 3' A/U additions,
negative binomial replicate counts, sequencing errors, rRNA and mRNA degradation
fragments, tRNA fragments (5' fragments and 3' fragments with the CCA tail), piRNA-like
reads, contaminants and adapter dimers. The simulated cases include NUMTs, a small RNA
across the origin, sense/antisense pairs, neighbouring and overlapping small RNAs,
SNPs, heteroplasmy, weakly expressed and sample-specific small RNAs.

```bash
pip install -e .                            # once, to get the smithhunter command
smithhunter simulate --outdir sim
cd sim
snakemake -s ../workflow/Snakefile --cores 2
smithhunter evaluate                        # writes evaluation/summary.txt and tables
```

`smithhunter simulate --help` lists the options (number of samples and reads, genome
sizes, degradation level, error rate, paired-end reads, seed). The first five planted 8mer sites also pair with guide positions 9-21; `evaluate` reports
whether the scored sites are supported (`support` >= 1): in the CI simulation the extended
sites are supported 4/4, the bare seed matches 0/11 (a seed alone is not a target) and 0/170
of the other sites (false support stays below the 5% FDR). `evaluate` reports for
each simulated small RNA whether it was recovered as its own locus, merged with
another, split, missed, or filtered out, and lists passing loci that contain no
simulated small RNA.

## Known limitations

- On simulated data, every small RNA has the expected outcome with up to 30% of the
  reads from degradation. Real degradation is less uniform than the simulated one;
  satellite merging and the isolation filter handle the secondary peaks seen on real
  data, but thresholds may need tuning on deep libraries.
- A peak locus also contains degradation fragments that start at the same position,
  so its `end` can extend beyond the small RNA; `rep_sequence` is the reliable sequence.
- Reads from 3' tRNA fragments carry a non-templated CCA and usually do not align.
- 5' tRNA fragments look like genuine small RNAs: use `annotation_class` to tell them apart.

## Testing

```bash
pytest                                              # unit tests, from the repository root
cd .test && snakemake -s ../workflow/Snakefile --cores 2   # small real dataset
```

The `.test` data are a heavily subset *Ceratitis capitata* example from SmithHunter;
they check that the workflow runs, not that the results are meaningful.

## Roadmap

1. Module B: precursor folding, then the opt-in bacterial mode. (Site scoring with
   IntaRNA and RNAhybrid and decoy q-values is done.)
2. Report: plots of module B results.

## Citation

If you use SmithHunter2, please cite the SmithHunter paper:
Marturano, Carli et al., *BMC Bioinformatics* 25:286 (2024),
and the smithRNA papers: Pozzi et al., *Mol Biol Evol* 34:1960 (2017);
Passamonti et al., *Sci Rep* 10:8219 (2020).

## License

GPL-3.0-or-later, as SmithHunter.
