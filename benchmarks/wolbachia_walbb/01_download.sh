#!/usr/bin/env bash
# Download the wAlbB benchmark data into $BENCH (default ~/bench_walbb):
#   - small RNA-seq of Aedes aegypti WB2 (wAlbB-infected) and WB2.tet (tetracycline-cleared),
#     whole adult females, 3 pools x 3 ages per line, Illumina HiSeq 2500 SE 50
#     (Sci Rep 2022, doi:10.1038/s41598-022-19574-x)
#   - wAlbB genome NZ_CP031221.1 (+ GFF3), Ae. aegypti AaegL5.0 (nuclear and mito split)
# Light on CPU, heavy on network/disk (~1.3 GB genome + 18 fastq). Safe to re-run: existing
# files are skipped.
set -euo pipefail

BENCH=${BENCH:-$HOME/bench_walbb}
mkdir -p "$BENCH"/{data,fastq}
cd "$BENCH"

ENA="https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS="run_accession,study_accession,sample_accession,sample_title,sample_alias,experiment_title,library_name,library_strategy,library_layout,read_count,base_count,fastq_ftp"
SEED_RUN=SRR12893564   # one of the WB2.tet runs cited in the paper

fetch() {  # url out
    [ -s "$2" ] && { echo "have $2"; return; }
    echo "get  $2"
    curl -fsSL --retry 3 -o "$2.part" "$1" && mv "$2.part" "$2"
}

# 1. run table of the whole study
study=$(curl -fsSL "$ENA?accession=$SEED_RUN&result=read_run&fields=study_accession&format=tsv" \
        | awk -F'\t' 'NR==2{print $2}')
[ -n "$study" ] || { echo "could not resolve the study of $SEED_RUN" >&2; exit 1; }
echo "study: $study"
fetch "$ENA?accession=$study&result=read_run&fields=$FIELDS&format=tsv" runs.tsv
echo "runs in study: $(($(wc -l < runs.tsv) - 1))  (see runs.tsv)"

# 2. fastq (single end: one file per run)
tail -n +2 runs.tsv | while IFS=$'\t' read -r run rest; do
    url=$(awk -F'\t' -v r="$run" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i} $1==r{print $c["fastq_ftp"]}' runs.tsv)
    first=${url%%;*}
    [ -n "$first" ] || { echo "no fastq for $run" >&2; continue; }
    fetch "https://$first" "fastq/${run}.fastq.gz"
done

# 3. genomes
EUTILS="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NZ_CP031221.1"
fetch "$EUTILS&rettype=fasta&retmode=text" data/wAlbB.fasta
fetch "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=NZ_CP031221.1&retmode=text" \
      data/wAlbB.gff3 || echo "WARNING: no GFF3 for wAlbB; remove 'annotation' from the configs" >&2

AAEG="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/204/515/GCF_002204515.2_AaegL5.0/GCF_002204515.2_AaegL5.0_genomic.fna.gz"
fetch "$AAEG" data/AaegL5.fna.gz
if [ ! -s data/AaegL5_mt.fasta ] || [ ! -s data/AaegL5_nuc.fasta ]; then
    echo "split AaegL5 into nuclear and mitochondrial (NC_035159.1)"
    gzip -dc data/AaegL5.fna.gz | awk '
        /^>/ {mt = ($1 == ">NC_035159.1")}
        { if (mt) print > "data/AaegL5_mt.fasta"; else print > "data/AaegL5_nuc.fasta" }'
fi
head -c1 data/wAlbB.fasta | grep -q '>' || { echo "data/wAlbB.fasta is not FASTA: delete it and retry" >&2; exit 1; }
if [ -s data/wAlbB.gff3 ] && ! head -1 data/wAlbB.gff3 | grep -q '^##gff-version'; then
    mv data/wAlbB.gff3 data/wAlbB.gff3.bad
    echo "WARNING: data/wAlbB.gff3 is not GFF3 (kept as .bad); remove 'annotation' from the configs" >&2
fi
grep -c '^>' data/AaegL5_mt.fasta | grep -qx 1 || { echo "mito contig NC_035159.1 not found" >&2; exit 1; }

echo
echo "done. next: python3 $(dirname "$0")/02_samples.py   (BENCH=$BENCH)"
