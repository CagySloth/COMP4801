#!/usr/bin/env bash
set -euo pipefail

source .venv/bin/activate

OUTDIR="output/prelim/whats_maxcov_sweep_vcf"
mkdir -p "$OUTDIR"

# 1) Simulate ONE dataset (fixed for fairness across max_coverage)
# Choose a high-coverage setting so read selection matters.
N=500
L=50
COV=30
R=$(( COV * N / L ))

ERROR=0.01
MISS=0.0
SEED=42

SIM_PREFIX="$OUTDIR/sim"
python -m dataset.simulate \
  -p 2 -n "$N" -r "$R" -l "$L" \
  -e "$ERROR" -m "$MISS" \
  --seed "$SEED" \
  -o "$SIM_PREFIX"

READS="${SIM_PREFIX}.reads.npz"
TRUTH="${SIM_PREFIX}.haplotypes.tsv"
VCF="${SIM_PREFIX}.vcf"

# Sanity check: make sure VCF exists
if [[ ! -f "$VCF" ]]; then
  echo "ERROR: VCF not found at $VCF"
  exit 1
fi

# 2) Sweep max-coverage (5 steps)
MAXCOVS=(10 15 20 25 30 35 40)

for MC in "${MAXCOVS[@]}"; do
  PREFIX="$OUTDIR/wh_maxcov${MC}"

  echo "== Running WhatsHap (VCF mode) with max_coverage=$MC =="

  python -m algorithms.cli.phase diploid-whats \
    -i "$READS" \
    --vcf "$VCF" \
    --solver hapchat \
    --max-coverage "$MC" \
    --output-prefix "$PREFIX"

  # Compute accuracy against truth haplotypes (your benchmark metric)
  python -m benchmark.benchmark_accuracy \
    --truth "$TRUTH" \
    --pred  "${PREFIX}.haplotypes.tsv" \
    --output "${PREFIX}.accuracy.json"
done

echo "Done. Results in: $OUTDIR"
echo "Next: python scripts/collect_whats_maxcov_sweep.py (update its OUTDIR path if needed)"
