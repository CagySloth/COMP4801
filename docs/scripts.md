## 📂 `scripts` – Utility Wrappers for Command-Line Use

This directory contains shell scripts that simplify frequent tasks like running simulations, benchmarking all algorithms, or executing the full pipeline in one go.

---

### 📄 `simulate.sh`

Wrapper for generating synthetic read and haplotype data using `dataset/simulate.py`.

#### Usage

```bash
./scripts/simulate.sh
```

#### Behavior

* Defines ploidy, variant count, error rate, etc.
* Saves output files with prefix `sim_output`

You can modify the script parameters directly inside the file.

---

### 📄 `benchmark_all.sh`

Runs `benchmark_runner.py` with a full sweep over common parameters and all algorithms.

#### Usage

```bash
./scripts/benchmark_all.sh
```

#### Behavior

* Benchmarks all diploid and polyploid algorithms
* Sweeps over number of reads
* Stores output in `results/` directory

Edit the script to change the `--vary` parameter or algorithms list.

---

### 📄 `run_pipeline.sh`

Runs the complete pipeline: simulate → convert → phase → evaluate

#### Usage

```bash
./scripts/run_pipeline.sh
```

#### Behavior

1. Simulates one dataset
2. Converts sparse `.tsv` to `.npz`
3. Phases with one algorithm (edit script to change)
4. Evaluates accuracy

This is a good sanity check or demo for your repo.