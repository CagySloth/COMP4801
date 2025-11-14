# 🧬 COMP4801: Haplotype Phasing Benchmark Suite

This project benchmarks multiple algorithms for phasing haplotypes using synthetic read data. It supports diploid and polyploid models, sparse/dense read formats, and synthetic simulations.

---

## 📦 Installation

### ✅ Prerequisites

* Python 3.11+
* Recommended: Create a virtual environment

```bash
python3.11 -m venv .venv
source .venv/bin/activate
```

### 📥 Install Dependencies

```bash
pip install -r requirements.txt
```

Alternatively, if using `pyproject.toml`:

```bash
pip install .
```

---

## 🚀 Quick Start

### 🔁 Run a Full Benchmark (All Algorithms)

```bash
python benchmarking/benchmark_runner.py \
  --algorithms diploid_em diploid_mst polyploid_em polyploid_spectral \
  --ploidy 4 \
  --num-variants 1000 \
  --num-reads 3000 \
  --vary error_rate \
  --vary-values 0.005 0.01 0.02 \
  --num-runs 3 \
  --outdir results/benchmark_run
```

This will:

1. Simulate datasets
2. Run each algorithm
3. Score accuracy
4. Save `benchmark_summary.json`

---

## 🧪 Example Pipelines

### Simulate and Phase with Diploid EM

```bash
python dataset/simulate.py -p 2 -n 1000 -r 5000 -l 50 -e 0.01 -m 0.05 -o sim/test

python -m algorithms.cli.convert -i sim/test.reads.sparse.tsv --to-npz sim/test.reads.npz

python -m algorithms.cli.phase diploid_em \
  -i sim/test.reads.npz -o sim/test.out
```

### Evaluate Accuracy

```bash
python benchmarking/benchmark_accuracy.py \
  --truth sim/test.haplotypes.tsv \
  --pred sim/test.out.haplotypes.tsv \
  --output sim/test.out.accuracy.json
```

---

## 🧠 Scripts Overview

Each directory has its own documentation under [`docs/scripts.md`](docs/scripts.md), but here are the highlights:

| Script                               | Description                                      |
| ------------------------------------ | ------------------------------------------------ |
| `dataset/simulate.py`                | Generates synthetic haplotypes and reads         |
| `algorithms/cli/phase.py`            | Runs phasing with selected algorithm             |
| `algorithms/cli/convert.py`          | Converts TSV ↔ NPZ formats                       |
| `benchmarking/benchmark_runner.py`   | Benchmarks all algorithms over a parameter sweep |
| `benchmarking/benchmark_accuracy.py` | Scores phasing accuracy                          |

---

## 🧪 Testing

Run all tests using `pytest`:

```bash
pytest tests/
```

Tests cover I/O, algorithms, and pipeline functionality.

---

## 📁 Repo Structure

```
COMP4801/
├── algorithms/          # Diploid & polyploid algorithms + CLI
├── dataset/             # Synthetic data generators
├── benchmarking/        # Accuracy and performance benchmarks
├── scripts/             # Shell wrappers for quick usage
├── docs/                # Script usage documentation
├── tests/               # Unit and integration tests
├── README.md
├── requirements.txt / pyproject.toml
```