# 🧬 COMP4801: Haplotype Phasing Benchmark Suite

A research-oriented benchmark suite for **haplotype phasing** on **synthetic diploid and polyploid** data.

This repo supports two complementary workflows:

1. **Matrix track (fast + controlled)**  
   Simulate variant-level read matrices (`.reads.npz`) → phase → evaluate TSV accuracy → run sweeps.

2. **Long-read end-to-end track (WhatsHap-like)**  
   Generate reference + diploid truth → simulate long reads (FASTQ) → align (BAM) → call variants (VCF) → phase (phased VCF) → evaluate.

---

## What’s included

- **Matrix simulator** (`dataset.simulate`)
  - Diploid simulation additionally writes a **minimal unphased VCF** (`.vcf`) derived from truth haplotypes.
- **Long-read pipeline generators** (`dataset.longread.*`)
  - `dataset.longread.reference` → reference FASTA (+ optional region metadata)
  - `dataset.longread.truth` → truth phased VCF + haplotype FASTAs
  - `dataset.longread.readsim` → long reads FASTQ (+ per-read truth labels)
- **Phasing algorithms** (`algorithms`)
  - Diploid: `diploid-em`, `diploid-mst`, `diploid-whats`, `diploid-whats-bam`
  - Polyploid: `polyploid-em`, `polyploid-spectral`
- **Vendored WhatsHap core** (`vendor/whatshap_core`)
  - Used by `diploid-whats` and `diploid-whats-bam`.
  - VCF-mode supports two solver backends:
    - `--solver whatshap` → **PedigreeDPTable** (WhatsHap DP-style solver)
    - `--solver hapchat` → **HapChatCore** (MEC-style solver)
- **Benchmarking & evaluation** (`benchmark`)
  - `benchmark.benchmark_runner` runs matrix-track sweeps (simulate → phase → score).
  - `benchmark.benchmark_accuracy` computes TSV accuracy.
  - `benchmark.vcf_phase_eval` evaluates **phased VCF vs phased truth VCF** (block-flip aware).
  - `benchmark.longread_pipeline_runner` runs the full long-read pipeline (Steps 1–7).

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

### Install the project (Python deps)

Recommended:

```bash
pip install -e .
```

Minimal (may not include VCF-mode deps):

```bash
pip install -r requirements.txt
```

### Vendored WhatsHap core

The `diploid-whats` and `diploid-whats-bam` implementations use the vendored WhatsHap package in `vendor/whatshap_core`.

Install the vendored package (NOT the PyPI whatshap):

```bash
pip install -e vendor/whatshap_core
```

Quick import check (paths should point under `vendor/whatshap_core`):

```bash
python - <<'PY'
import whatshap as wh
from whatshap import core, readselect
print("whatshap:", wh.__file__)
print("core:", core.__file__)
print("readselect:", readselect.__file__)
PY
```

### System tools (required for long-read pipeline)

The long-read end-to-end track uses external command-line tools:

- `minimap2`
- `samtools`
- `bcftools`

On Ubuntu/Debian:

```bash
sudo apt-get update
sudo apt-get install -y minimap2 samtools bcftools
```

Check versions:

```bash
minimap2 --version
samtools --version
bcftools --version
```

---

## Quick start (matrix track)

### 1) Simulate a dataset

```bash
python -m dataset.simulate \
  -p 2 -n 200 -r 200 -l 60 \
  -e 0.01 -m 0.0 \
  --seed 0 \
  -o output/demo
```

Outputs:

- `output/demo.haplotypes.tsv` (truth)
- `output/demo.reads.sparse.tsv` (sparse fragments)
- `output/demo.reads.npz` (dense matrix)
- `output/demo.vcf` (diploid only; **unphased** GT from truth)

### 2) Phase (diploid-whats in VCF-mode)

```bash
python -m algorithms.cli.phase diploid-whats \
  -i output/demo.reads.npz \
  --vcf output/demo.vcf \
  --output-prefix output/demo_phased \
  --solver whatshap
```

Writes:

- `output/demo_phased.haplotypes.tsv`
- `output/demo_phased.assignments.tsv`
- `output/demo_phased.summary.json`
- `output/demo_phased.phased.vcf` (VCF-mode only)

Notes on WhatsHap VCF-mode behavior in this repo:

- Only **heterozygous** variants (from VCF GT) are phased.
- **PS (phase set)** is assigned by **connectivity** (variants connected by selected reads share PS).
- TSV output is “dense” (0/1 at every site) for the TSV scorer; the phased VCF is the authoritative WhatsHap-like output.

### 3) Evaluate TSV accuracy

```bash
python -m benchmark.benchmark_accuracy \
  --truth output/demo.haplotypes.tsv \
  --pred  output/demo_phased.haplotypes.tsv \
  --output output/demo_phased.accuracy.json
```

---

## Quick start (long-read end-to-end track)

### Run Steps 1–7 automatically

```bash
python -m benchmark.longread_pipeline_runner \
  --prefix output/lr_demo \
  --seed 0 \
  --ref-length 20000 \
  --num-snps 200 \
  --het-rate 0.8 \
  --num-reads 200 \
  --min-len 2000 \
  --max-len 6000
```

Key outputs (all under the chosen prefix):

- `.ref.fasta` (reference)
- `.truth.vcf` (+ optional `.gz/.tbi`)
- `.hap1.fasta`, `.hap2.fasta`
- `.reads.fastq`, `.reads.truth.tsv`
- `.bam` (+ `.bai`)
- `.called.vcf.gz` (+ `.tbi`)
- `.ws.phased.vcf` (phased output)
- `.ws.eval.json` (evaluation vs truth)
- `.pipeline.json` (end-to-end report)

---

## Benchmark sweeps (matrix track)

Example: sweep `error-rate` over three values, 2 runs each:

```bash
python -m benchmark.benchmark_runner \
  --algorithms diploid-em diploid-mst diploid-whats \
  --ploidy 2 \
  --num-variants 500 \
  --num-reads 500 \
  --read-length 60 \
  --missing-rate 0.0 \
  --vary error_rate \
  --vary-values 0.005 0.01 0.02 \
  --num-runs 2 \
  --outdir output/bench_error_sweep
```

---

## Repository layout

```text
COMP4801/
├── algorithms/          # Diploid & polyploid algorithms + CLI
├── dataset/             # Dataset generators (matrix + longread)
├── benchmark/           # Benchmark runners + evaluators
├── vendor/              # Vendored WhatsHap core (C++/Cython extension)
├── scripts/             # Convenience wrappers (may require editing)
├── docs/                # Documentation
├── tests/               # Unit + integration tests (pytest)
├── pyproject.toml / requirements.txt
└── README.md
```

---

## Testing

```bash
pytest tests/
```

Notable tests:

- `tests/test_vcf_mode_pipeline.py` — simulate → phase (VCF-mode) → evaluate
- `tests/test_pedigreedptable_vcf_mode.py` — regression for PedigreeDPTable path + PS semantics

---

## Troubleshooting

- **Missing minimap2/samtools/bcftools**
  - Install via apt (see Installation → System tools)
- **VCF variant count != reads matrix N** (matrix VCF-mode)
  - VCF-mode assumes simulator’s VCF records align 1:1 with columns in `reads.npz`.
- **Imported the wrong whatshap**
  - Ensure you installed the vendored package: `pip install -e vendor/whatshap_core`
  - Confirm `whatshap.__file__` points under `vendor/whatshap_core`.
