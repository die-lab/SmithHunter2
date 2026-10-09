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
  a policy you choose.
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
> Module B (target prediction) is being designed and is not yet available.

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
mamba env create -f environment.yml      # fastp, bowtie 1, snakemake, python
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
| `ends.min_five_score` | 0.5 | a locus passes if its 5' score is above this |
| `ends.min_three_score` | 0.0 | and its 3' score is at least this |
| `ends.min_isolation` | 0.5 | and its isolation is at least this |

The end score of a locus is 1 when one position holds at least half of the reads
(`n_thre`), 0.5 for two adjacent positions, and lower for scattered ends.

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
| `five_score`, `three_score` | end precision (see above); `*_dominant_fraction` = share of the top position |
| `isolation` | share of the reads starting near the peak that belong to the locus |
| `annotation_class`, `annotation_orientation` | most specific overlapping GFF3 feature and its strand relation |
| `expression_pass`, `ends_pass`, `pass` | filter results |

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
sizes, degradation level, error rate, paired-end reads, seed). `evaluate` reports for
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

1. Module B: target prediction (seed match classes, IntaRNA and RNAhybrid, empirical
   false discovery rate from shuffled small RNAs; seed rules for Argonaute-guided small
   RNAs, start-codon windows for bacterial genomes) and precursor folding.
2. HTML report.

## Citation

If you use SmithHunter2, please cite the SmithHunter paper:
Marturano, Carli et al., *BMC Bioinformatics* 25:286 (2024),
and the smithRNA papers: Pozzi et al., *Mol Biol Evol* 34:1960 (2017);
Passamonti et al., *Sci Rep* 10:8219 (2020).

## License

GPL-3.0-or-later, as SmithHunter.
