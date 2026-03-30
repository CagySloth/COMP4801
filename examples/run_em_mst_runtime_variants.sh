#!/usr/bin/env bash
set -euo pipefail

mkdir -p output/prelim

COV=20
L=60

for N in 100 250 500 1000 2000 400; do
  R=$(( COV * N / L ))
  echo "Running N=$N, R=$R, L=$L (approx cov=$COV)"

  python -m benchmark.benchmark_runner \
    --algorithms diploid-em diploid-mst \
    --ploidy 2 \
    --num-variants "$N" \
    --num-reads "$R" \
    --read-length "$L" \
    --error-rate 0.01 \
    --missing-rate 0.0 \
    --num-runs 3 \
    --outdir "output/prelim/em_mst_runtime_N${N}_cov${COV}"
done

# merge results for Excel/R
OUT="output/prelim/em_mst_runtime_variants.csv"
head -n 1 "output/prelim/em_mst_runtime_N200_cov${COV}/benchmark_summary.csv" > "$OUT"
for d in output/prelim/em_mst_runtime_N*_cov${COV}; do
  tail -n +2 "$d/benchmark_summary.csv" >> "$OUT"
done

echo "Done. Combined CSV at: $OUT"
