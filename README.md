# SmithHunter2

Discovery of small RNA loci from small RNA-seq data, across all the genomes present in a
sample (nuclear, mitochondrial, bacterial, ...). You declare every genome and choose which
ones to study; reads are assigned to the genome they align to best, reads from the other
genomes are counted and set aside, and reads that fit more than one genome are reported
separately.

SmithHunter2 is a rewrite of [SmithHunter](https://github.com/ESZlab/SmithHunter)
(Marturano, Carli et al., BMC Bioinformatics 25:286, 2024), which identifies candidate
smithRNAs (small mitochondrial highly-transcribed RNAs) and their nuclear targets.

> **Status: early development.** Module A (discovery) is implemented. Module B (target
> prediction) is not yet. See [CLAUDE.md](CLAUDE.md) for the design and the roadmap.

## What module A does

1. **Trimming** with fastp. Paired-end reads are merged, since small RNA inserts are shorter
   than the reads.
2. **Collapsing** of identical reads across samples; each unique sequence is aligned once.
3. **Competitive alignment** with bowtie 1 (`--best --strata`) on one reference holding all
   declared genomes. Circular genomes are extended so reads across the origin align.
4. **Origin assignment**: `assigned` (one focus genome), `excluded`, `ambiguous`,
   `too_many_hits`, `unmapped`, with per-sample and per-length summaries.
5. **Loci**: overlapping alignments on the same strand, quantified per sample as counts, RPM
   on the library and RPM on the genome.
6. **End precision**: 5' and 3' end scores, after `sharp_smith.R` from SmithHunter.
7. **Annotation** with the GFF3 of each genome (e.g. MITOS2 for mitochondria).

## Running

Requirements: Linux, conda/mamba.

```bash
git clone https://github.com/die-lab/SmithHunter2.git
cd SmithHunter2
mamba env create -f environment.yml
mamba activate smithhunter2
```

Edit `config/config.yaml` (genomes, adapters, thresholds) and `config/samples.tsv`, then:

```bash
snakemake -s workflow/Snakefile --cores 8
```

Alternatively, let Snakemake create one environment per rule:
`snakemake -s workflow/Snakefile --cores 8 --software-deployment-method conda`.

### Test run

```bash
cd .test
snakemake -s ../workflow/Snakefile --cores 2
pytest   # from the repository root: unit tests
```

The test data is the SmithHunter example: heavily subset *Ceratitis capitata* reads and
genomes, not meant to give biologically meaningful results.

## Main outputs

| File | Content |
|---|---|
| `discovery/loci.tsv` | every locus with coordinates, representative sequence, counts, RPM, end scores, annotation and filter flags |
| `discovery/candidates.fasta` | representative sequences of loci passing all filters |
| `discovery/locus_reads.tsv.gz` | the reads of every locus, for manual review |
| `origin/origin_summary.tsv` | reads per origin class and genome, per sample |
| `origin/length_distribution.tsv` | the same, split by read length |
| `qc/*.fastp.html` | trimming reports |

## Citation

If you use this work, please cite the SmithHunter paper above, and the smithRNA papers:
Pozzi et al., Mol Biol Evol 34:1960 (2017); Passamonti et al., Sci Rep 10:8219 (2020).

## License

GPL-3.0, as SmithHunter.
