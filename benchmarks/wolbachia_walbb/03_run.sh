#!/usr/bin/env bash
# Run SmithHunter2 module A on both lines, then the comparison.
#   BENCH=~/bench_walbb CORES=8 bash 03_run.sh            # both lines + compare
#   STEP=infected|tet|compare bash 03_run.sh              # one step only
# Heavy: bowtie-build of AaegL5 (1.3 Gb; ~1-3 h, ~4-8 GB RAM) runs once, in the infected run;
# the tet run reuses that index. The machine is shared: check `htop` first, and prefer
#   nohup bash 03_run.sh > run.log 2>&1 &
set -euo pipefail

BENCH=${BENCH:-$HOME/bench_walbb}
CORES=${CORES:-8}
STEP=${STEP:-all}
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
SNAKE="snakemake -s $REPO/workflow/Snakefile --cores $CORES --rerun-triggers mtime"

if [ -z "${CONDA_DEFAULT_ENV:-}" ] || [ "$CONDA_DEFAULT_ENV" != "smithhunter2" ]; then
    # shellcheck disable=SC1091
    source "$HOME/miniforge3/etc/profile.d/conda.sh" && conda activate smithhunter2
fi

run_line() {  # group
    cd "$BENCH/$1"
    echo "== $1: $(($(wc -l < config/samples.tsv) - 1)) libraries, $(date)"
    $SNAKE -n --quiet rules 2>/dev/null || true
    $SNAKE
    echo "== $1 done $(date)"
}

if [ "$STEP" = all ] || [ "$STEP" = infected ]; then
    run_line infected
fi

if [ "$STEP" = all ] || [ "$STEP" = tet ]; then
    # same genomes -> same reference: reuse the index (timestamps kept)
    if [ -d "$BENCH/infected/results/reference/index" ] && [ ! -d "$BENCH/tet/results/reference" ]; then
        mkdir -p "$BENCH/tet/results"
        cp -a "$BENCH/infected/results/reference" "$BENCH/tet/results/"
    fi
    cd "$BENCH/tet"
    # capture first: with pipefail, grep -q closing the pipe early can make the test false
    plan=$($SNAKE -n 2>/dev/null || true)
    if grep -q "bowtie_index" <<< "$plan"; then
        echo "tet run would rebuild the bowtie index; touching the copied reference" >&2
        $SNAKE --touch "results/reference/combined.fa" "results/reference/index" >/dev/null 2>&1 || true
    fi
    run_line tet
fi

if [ "$STEP" = all ] || [ "$STEP" = compare ]; then
    cd "$BENCH"
    python3 "$HERE/04_compare.py" --infected infected/results --tet tet/results \
        --genome "${FOCUS:-wolbachia}" --genome-fasta "${GENOME_FASTA:-data/wAlbB.fasta}" \
        ${POSITIVE:+--positive "$POSITIVE"} \
        ${V0_FASTA:+--v0-fasta "$V0_FASTA"} | tee compare_report.md
    echo "report: $BENCH/compare_report.md"
fi
