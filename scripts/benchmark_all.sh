#!/bin/bash

# Example benchmark sweep across algorithms and ploidy levels
# Modify these values as needed

PLOIDIES=(2 3 4)
ALGORITHMS=("em" "spectral")
NVAR=1000
NREADS=5000

for P in "${PLOIDIES[@]}"; do
    for ALGO in "${ALGORITHMS[@]}"; do
        echo "Running benchmark for ploidy=${P}, algo=${ALGO}"
        python benchmarking/benchmark_runner.py \
            --ploidy "$P" \
            --algorithm "$ALGO" \
            --num-variants "$NVAR" \
            --num-reads "$NREADS" \
            --output "bench_results/p${P}_${ALGO}.json" \
            --runs 3
    done
done
