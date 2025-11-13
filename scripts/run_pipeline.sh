#!/bin/bash

# Example end-to-end pipeline: simulate → convert → phase → evaluate
# Usage: ./scripts/run_pipeline.sh PREFIX ALGO PLOIDY

PREFIX=$1         # e.g. data/example1
ALGO=$2           # e.g. em, mst, spectral
PLOIDY=$3         # e.g. 2, 3, 4

if [ -z "$PREFIX" ] || [ -z "$ALGO" ] || [ -z "$PLOIDY" ]; then
    echo "Usage: $0 PREFIX ALGO PLOIDY"
    exit 1
fi

echo "Simulating dataset..."
python dataset/simulate.py --ploidy "$PLOIDY" --num-variants 1000 --num-reads 5000 --output "$PREFIX"

echo "Converting to .npz..."
python -m algorithms.cli.convert "$PREFIX.reads.sparse.tsv" --output "$PREFIX.reads.npz"

echo "Running phasing algorithm..."
python -m algorithms.cli.phase "$PREFIX.reads.npz" --algorithm "$ALGO" --ploidy "$PLOIDY" --output "$PREFIX.$ALGO"

echo "Evaluating..."
python benchmarking/benchmark_accuracy.py \
    --truth "$PREFIX.haplotypes.tsv" \
    --pred "$PREFIX.$ALGO.haplotypes.tsv" \
    --output "$PREFIX.$ALGO.accuracy.json"
