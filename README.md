# 🧬 COMP4801: Haplotype Phasing Benchmark Suite

A research-oriented benchmark suite for **haplotype phasing** on **synthetic diploid and polyploid** data.

The core workflow is:

1. **Simulate** haplotypes + reads (dense `.npz` and sparse `.tsv`)
2. **Phase** reads using one of the implemented algorithms (including a **vendored WhatsHap core**)
3. **Evaluate** phasing accuracy
4. **Sweep** parameters (error rate, missing rate, read depth, etc.) via the benchmark runner

---

## What’s included

- **Synthetic data generator** (`dataset.simulate`)
  - Diploid simulation additionally writes an **unphased VCF** (`.vcf`) derived from truth haplotypes.
- **Phasing algorithms** (`algorithms`)
  - Diploid: `diploid-em`, `diploid-mst`, `diploid-whats`
  - Polyploid: `polyploid-em`, `polyploid-spectral`
- **Vendored WhatsHap core** (`vendor/whatshap_core`)
  - Supports two solvers in VCF-mode:
    - `--solver whatshap` → **PedigreeDPTable** (WhatsHap’s default DP-style solver)
    - `--solver hapchat` → **HapChatCore** (MEC-style solver)
- **Benchmarking** (`benchmark`)
  - `benchmark.benchmark_runner` runs sweeps and collects runtime + accuracy.
  - `benchmark.benchmark_accuracy` computes accuracy metrics from TSV haplotypes.

---

## Installation

### Prerequisites

- Python **3.11+**
- A C/C++ build toolchain (recommended) if you need to rebuild the vendored WhatsHap extension:
  - Linux: `build-essential`, `python3-dev`
  - macOS: Xcode command line tools

### Create a virtual environment

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
```

### Install the project

Recommended (uses `pyproject.toml` dependencies):

```bash
pip install -e .
```

Minimal dependencies (note: does **not** include VCF-mode deps like `pysam/pyfaidx`):

```bash
pip install -r requirements.txt
```

### Install the vendored WhatsHap core

The `diploid-whats` algorithm depends on the vendored WhatsHap Python package living in `vendor/whatshap_core`.

```bash
pip install -e vendor/whatshap_core
```

Quick import check:

```bash
python - <<'PY'
from whatshap import core, readselect
print("OK", core, readselect)
PY
```

---

## Quick start

### 1) Simulate a dataset

```bash
python -m dataset.simulate       -p 2 -n 200 -r 200 -l 60       -e 0.01 -m 0.0       --seed 0       -o output/demo
```

Outputs:

- `output/demo.haplotypes.tsv` (truth)
- `output/demo.reads.sparse.tsv` (sparse fragments)
- `output/demo.reads.npz` (dense matrix, used by algorithms)
- `output/demo.vcf` (diploid only; **unphased** GT from truth)

### 2) Phase with a diploid algorithm

Example (WhatsHap in **VCF-mode**, default solver = **PedigreeDPTable**):

```bash
python -m algorithms.cli.phase diploid-whats       -i output/demo.reads.npz       --vcf output/demo.vcf       --output-prefix output/demo_phased       --solver whatshap
```

This writes:

- `output/demo_phased.haplotypes.tsv`
- `output/demo_phased.assignments.tsv`
- `output/demo_phased.summary.json`
- `output/demo_phased.phased.vcf` (VCF-mode only)

Notes on WhatsHap VCF-mode behavior in this repo:

- Only **heterozygous** variants (from the VCF GT) are phased.
- **PS (phase set)** is assigned by **connectivity** (variants connected by selected reads share PS). PS is reported as **leftmost variant index + 1**.
- The TSV output is “dense” (always 0/1 per site) for compatibility with the benchmark scorer; the VCF output is the authoritative “WhatsHap-like” phasing signal (phased GT uses `|`, unphased uses `/` and PS=`.`).

Example (diploid EM):

```bash
python -m algorithms.cli.phase diploid-em       -i output/demo.reads.npz       --output-prefix output/demo_em
```

### 3) Evaluate accuracy

```bash
python -m benchmark.benchmark_accuracy       --truth output/demo.haplotypes.tsv       --pred  output/demo_phased.haplotypes.tsv       --output output/demo_phased.accuracy.json
```

---

## Run benchmark sweeps

The benchmark runner simulates datasets and runs selected algorithms across multiple settings.

Example: sweep `error-rate` over three values, with 2 runs each:

```bash
python -m benchmark.benchmark_runner       --algorithms diploid-em diploid-mst diploid-whats       --ploidy 2       --num-variants 500       --num-reads 500       --read-length 60       --missing-rate 0.0       --vary error_rate       --vary-values 0.005 0.01 0.02       --num-runs 2       --outdir output/bench_error_sweep
```

Output:

- per-run, per-algorithm result prefixes in `output/bench_error_sweep/`
- `output/bench_error_sweep/benchmark_summary.json`

The runner automatically enables VCF-mode for `diploid-whats` when a simulated `*.vcf` exists.

---

## Repository layout

```text
COMP4801/
├── algorithms/          # Diploid & polyploid algorithms + CLI
├── dataset/             # Synthetic data generators
├── benchmark/           # Accuracy and performance benchmarks
├── vendor/              # Vendored WhatsHap core (C++/Cython extension)
├── scripts/             # Convenience shell wrappers (may require editing)
├── docs/                # Additional documentation
├── tests/               # Unit + integration tests (pytest)
├── pyproject.toml / requirements.txt
└── README.md
```

---

## Testing

```bash
pytest tests/
```

Notable integration tests:

- `tests/test_vcf_mode_pipeline.py` — end-to-end simulate → phase (VCF-mode) → evaluate
- `tests/test_pedigreedptable_vcf_mode.py` — unit regression for **PedigreeDPTable** path and PS semantics

---

## Troubleshooting

- **`Phasing did not write phased VCF`**
  - Ensure you are running `diploid-whats` with `--vcf ...` (VCF-mode).
- **`VCF variant count != reads matrix N`**
  - In this repo, VCF-mode assumes the simulator’s VCF records align 1:1 with columns in `reads.npz`.
- **Import errors in `whatshap`**
  - Make sure you installed the vendored package:
    `pip install -e vendor/whatshap_core`
