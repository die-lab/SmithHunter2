# Benchmark: Wolbachia wAlbB in *Aedes aegypti* (module A on a symbiont)

Why: on the maintainer's Wolbachia data SmithHunter2 module A returned 300+ loci, most passing,
where v0 found 105 clusters. This benchmark tells apart the possible causes on public data that
includes a **negative control**.

Data (Sci Rep 2022, doi:10.1038/s41598-022-19574-x): small RNA-seq of whole adult females,
line WB2 (wAlbB-infected) and WB2.tet (same line, tetracycline-cleared), 3 pools x 2, 6, 12 days
= 9 libraries per line, HiSeq 2500 SE 50. Genomes: wAlbB NZ_CP031221.1 (focus, circular, NCBI
GFF3), AaegL5.0 GCF_002204515.2 split into nuclear and mito NC_035159.1 (both exclude).

## Run (on the server)

Ask before steps 3-4: they are heavy and the machine is shared.

```bash
source ~/miniforge3/etc/profile.d/conda.sh && conda activate smithhunter2
B=benchmarks/wolbachia_walbb          # from the repo root
export BENCH=~/bench_walbb            # data and results live here, outside the repo

bash $B/01_download.sh                # 1. ~1.3 GB genome + 18 fastq; network only
python3 $B/02_samples.py --threads 8  # 2. splits WB2 / WB2.tet, checks the adapter, writes
                                      #    $BENCH/{infected,tet}/config/; stops if the lines are unclear
nohup env CORES=8 bash $B/03_run.sh > $BENCH/run.log 2>&1 &   # 3. both lines + report
                                      #    bowtie-build AaegL5 once (~1-3 h, ~4-8 GB RAM), reused by tet
cat $BENCH/compare_report.md          # 4. read the report
```

Optional, the original pipeline on the same libraries (about a day per line; set the v0
parameters of the run that gave 105 clusters):

```bash
V0_DIR=~/prog/SmithHunter V0_I=0.95 V0_S=0.5 V0_M=8 bash $B/v0_prepare.sh infected   # prints the command
RUN=1 V0_DIR=... bash $B/v0_prepare.sh infected                                      # runs it
V0_FASTA=$(find $BENCH/v0_infected -name presumptive_smithRNAs.fa) STEP=compare bash $B/03_run.sh
```

`03_run.sh` steps can be run alone: `STEP=infected|tet|compare`.

## Reading `compare_report.md`

| what the report shows | cause | fix |
|---|---|---|
| many passing loci in **tet** | host reads placed on wAlbB (assignment, or contigs shared with the host) | assignment (`severe`/`evil`), check the loci in tet |
| tet clean, infected passing loci mostly on rRNA / tRNA / tmRNA / RNase P | degradation fragments of abundant structural RNAs pass a fixed RPM filter | filter structural RNAs and/or raise `min_rpm` (sweep table) |
| tet clean, passing loci spread over CDS / intergenic, sweep falls fast with `min_rpm` | expression filter too lenient for the library depth | raise `min_rpm` / `min_count`; consider a v0-like relative threshold |
| v0 presumptive count close to "peak loci passing expression" | 105 vs 300 is only peak splitting (v0 clusters by sequence) | none: compare loci, not counts |

Notes: v0 module A has no end filter (that is `sharp_smith.R`), so its presumptive list compares
with "peak loci passing expression", not with "pass = 1". The sweep table recomputes the
expression filter on gross counts (no background subtraction): read it as a trend.
The positive control of the literature, WsnRNA-46 and WsnRNA-49 (Mayoral et al. 2014, PNAS,
wMelPop-CLA), is a different strain; check by sequence if a wAlbB homolog is among the loci.

Tested here only on simulated data (scripts, configs, index reuse, report); download and
the real run are untested.
