#!/usr/bin/env bash
# Prepare (and optionally run) the original SmithHunter (ESZlab, "v0") on the same libraries.
#   BENCH=~/bench_walbb V0_DIR=~/prog/SmithHunter bash v0_prepare.sh [infected|tet]
#   RUN=1 ... bash v0_prepare.sh infected      # also launch it (about a day, bowtie2-build of AaegL5)
# Parameters: set them to the values of the run that gave 105 clusters, so the two tools are
# compared, not two parameter sets. Defaults: v0's own defaults for -I and -S, and -M 8 so that
# v0 asks for the same replicates as SmithHunter2 (n - 1 of 9).
# v0 treats the focus genome as "_mit" and filters only against "_nuc": host nuclear and host
# mito are concatenated into _nuc here. v0 module A output (presumptive_smithRNAs.fa) has no
# end filter yet (sharp_smith.R is a separate step): compare its count with "peak loci passing
# expression" in compare_report.md, and with "pass = 1" only after sharp_smith.R.
set -euo pipefail

BENCH=${BENCH:-$HOME/bench_walbb}
GROUP=${1:-infected}
V0_DIR=${V0_DIR:-$HOME/prog/SmithHunter}
V0_I=${V0_I:-0.95}
V0_S=${V0_S:-0.5}
V0_M=${V0_M:-8}
THREADS=${THREADS:-8}
ADAPTER=$(awk '/adapter_r1:/{gsub(/"/,"",$2); print $2}' "$BENCH/$GROUP/config/config.yaml")
W="$BENCH/v0_$GROUP"

mkdir -p "$W/rawdata"
[ -e "$W/WALB_mit.fasta" ] || cp "$BENCH/data/wAlbB.fasta" "$W/WALB_mit.fasta"
[ -s "$W/WALB_nuc.fasta" ] || cat "$BENCH/data/AaegL5_nuc.fasta" "$BENCH/data/AaegL5_mt.fasta" > "$W/WALB_nuc.fasta"
tail -n +2 "$BENCH/$GROUP/config/samples.tsv" | while IFS=$'\t' read -r sample fq1 _; do
    ln -sf "$fq1" "$W/rawdata/${sample}_1.fastq.gz"
done
echo "v0 working dir: $W ($(ls "$W/rawdata" | wc -l) libraries)"

CMD=(bash -i "$V0_DIR/smithHunterA.sh" -O WALB -W "$W" -F "$W/rawdata" -T SE
     -a "${ADAPTER:-TGGAATTCTCGGGTGCCAAGG}" -I "$V0_I" -S "$V0_S" -M "$V0_M" -t "$THREADS")
echo "command: ${CMD[*]}"
echo "output:  find $W -name presumptive_smithRNAs.fa   (then: V0_FASTA=<that file> STEP=compare bash 03_run.sh)"
if [ "${RUN:-0}" = 1 ]; then
    cd "$W" && "${CMD[@]}" > "$W/v0.log" 2>&1
fi
