#!/usr/bin/env bash
set -euo pipefail

OUTDIR="output/prelim/plot2_em_mst_accuracy_error"
mkdir -p "$OUTDIR"

python -m benchmark.benchmark_runner \
  --algorithms diploid-em diploid-mst \
  --ploidy 2 \
  --num-variants 500 \
  --num-reads 500 \
  --read-length 60 \
  --missing-rate 0.0 \
  --vary error_rate \
  --vary-values 0 0.002 0.005 0.01 0.02 0.05 \
  --num-runs 5 \
  --outdir "$OUTDIR"

echo "CSV ready: $OUTDIR/benchmark_summary.csv"
