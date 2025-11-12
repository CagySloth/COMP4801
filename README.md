# COMP4801
source .venv/bin/activate

## 📜 Scripts Documentation

This section describes the command-line tools included in this project. These scripts allow you to simulate haplotype sequencing data, convert TSV files to NumPy formats, and run various phasing algorithms.

---

### 🔧 `simulate_phasing_data.py` – Synthetic Data Generator

**Location:** `dataset/1/simulate_phasing_data.py`

Generates synthetic haplotype data and sparse sequencing reads for benchmarking phasing algorithms.

#### **Usage**

```bash
python dataset/1/simulate_phasing_data.py [options]
```

#### **Example**

```bash
python dataset/1/simulate_phasing_data.py -p 2 -n 1000 -r 5000 -l 50 -e 0.02 -m 0.1 -o output/sim1
```

#### **Parameters**

| Flag                        | Type  | Default | Description                                   |
| --------------------------- | ----- | ------- | --------------------------------------------- |
| `-p`, `--ploidy`            | int   | 2       | Number of haplotypes (e.g. 2 for diploid)     |
| `-n`, `--num-variants`      | int   | 1000    | Number of variant positions                   |
| `-r`, `--num-reads`         | int   | 20000   | Number of reads to simulate                   |
| `-l`, `--read-length`       | int   | 25      | Length of each read in variant sites          |
| `-e`, `--error-rate`        | float | 0.01    | Per-site sequencing error rate                |
| `-m`, `--missing-rate`      | float | 0.05    | Per-site missing data rate                    |
| `--maf-alpha`, `--maf-beta` | float | 0.4     | Beta distribution params for allele frequency |
| `--allow-monomorphic`       | flag  | off     | Include sites with no variation               |
| `-s`, `--seed`              | int   | None    | Random seed for reproducibility               |
| `-o`, `--output-prefix`     | str   | "sim"   | Prefix for output files                       |

#### **Outputs**

* `<prefix>.haplotypes.tsv`
* `<prefix>.haplotypes.npz`
* `<prefix>.reads.sparse.tsv`
* `<prefix>.summary.json`

---

### 📎 `read_tsv.py` – TSV Reader and Converter

**Location:** `algorithms/read_tsv.py`

Reads sparse or dense TSV read data and optionally converts it to `.npz` format for use in phasing algorithms.

#### **Usage**

```bash
python algorithms/read_tsv.py -i <input.tsv> [--mode auto|sparse|dense] [--to-npz file] [--to-dense-npz file] [--num-variants N] [--print-first N]
```

#### **Example**

```bash
python algorithms/read_tsv.py -i dataset/0\ -\ diploid\ only/demo.reads.tsv --to-npz demo_reads.npz --print-first 3
```

#### **Parameters**

| Flag             | Type | Default    | Description                                                     |
| ---------------- | ---- | ---------- | --------------------------------------------------------------- |
| `-i`, `--input`  | str  | *Required* | Input TSV file (dense or sparse format)                         |
| `--mode`         | str  | `auto`     | Force parsing mode (`dense`, `sparse`, or `auto`)               |
| `--to-npz`       | str  | None       | Path to output `.npz` file                                      |
| `--to-dense-npz` | str  | None       | Output `.npz` file after sparse-to-dense conversion             |
| `--num-variants` | int  | None       | Number of variants (required for sparse reads if not inferable) |
| `--print-first`  | int  | 0          | Print N example reads for verification                          |

#### **Output**

* Summary printed to console
* Optional `.npz` files (sparse and/or dense)

---

### 🧬 `diploid_mst_phase.py` – MST-Based Diploid Phasing

**Location:** `algorithms/diploid_mst_phase.py`

Phases diploid reads using a graph-based maximum spanning tree approach.

#### **Usage**

```bash
python algorithms/diploid_mst_phase.py -i reads.npz -o output_prefix [--min-overlap N] [--min-het-minor N]
```

#### **Example**

```bash
python algorithms/diploid_mst_phase.py -i demo_reads.npz -o output/mst_phase
```

#### **Parameters**

| Flag                    | Type | Default    | Description                                               |
| ----------------------- | ---- | ---------- | --------------------------------------------------------- |
| `-i`, `--input`         | str  | *Required* | Input `.npz` reads file                                   |
| `-o`, `--output-prefix` | str  | *Required* | Prefix for output files                                   |
| `--min-overlap`         | int  | 3          | Minimum shared reads between variant pairs                |
| `--min-het-minor`       | int  | 1          | Minimum minor allele count to call a variant heterozygous |

#### **Outputs**

* `<prefix>.haplotypes.tsv`
* `<prefix>.assignments.tsv`
* `<prefix>.summary.json`

---

### 🧪 `diploid_em_cluster.py` – EM-Based Diploid Phasing

**Location:** `algorithms/diploid_em_cluster.py`

Performs hard Expectation-Maximization clustering to phase diploid reads.

#### **Usage**

```bash
python algorithms/diploid_em_cluster.py -i reads.npz -o output_prefix [--max-iters K] [--tol-iters T] [-s SEED]
```

#### **Example**

```bash
python algorithms/diploid_em_cluster.py -i demo_reads.npz -o output/em_phase
```

#### **Parameters**

| Flag                    | Type | Default    | Description                                    |
| ----------------------- | ---- | ---------- | ---------------------------------------------- |
| `-i`, `--input`         | str  | *Required* | Input `.npz` reads file                        |
| `-o`, `--output-prefix` | str  | *Required* | Prefix for output files                        |
| `--max-iters`           | int  | 30         | Maximum number of EM iterations                |
| `--tol-iters`           | int  | 2          | Number of stable iterations before convergence |
| `-s`, `--seed`          | int  | None       | Random seed for initialization                 |

#### **Outputs**

* `<prefix>.haplotypes.tsv`
* `<prefix>.assignments.tsv`
* `<prefix>.summary.json`

---

### 🌈 `polyploid_spectral_cluster.py` – Spectral Clustering for Polyploids

**Location:** `algorithms/polyploid_spectral_cluster.py`

Performs spectral clustering of reads into K haplotype groups using read-read similarity.

#### **Usage**

```bash
python algorithms/polyploid_spectral_cluster.py -i reads.npz -k PLOIDY -o output_prefix [--min-overlap N] [-s SEED]
```

#### **Example**

```bash
python algorithms/polyploid_spectral_cluster.py -i demo_reads.npz -k 4 -o output/spectral_phase
```

#### **Parameters**

| Flag                    | Type | Default    | Description                             |
| ----------------------- | ---- | ---------- | --------------------------------------- |
| `-i`, `--input`         | str  | *Required* | Input `.npz` reads file                 |
| `-k`, `--ploidy`        | int  | *Required* | Number of haplotypes (K)                |
| `-o`, `--output-prefix` | str  | *Required* | Prefix for output files                 |
| `--min-overlap`         | int  | 3          | Minimum shared positions for similarity |
| `-s`, `--seed`          | int  | None       | Random seed for k-means step            |

#### **Outputs**

* `<prefix>.haplotypes.tsv`
* `<prefix>.assignments.tsv`
* `<prefix>.summary.json`

---

### 🧪 `polyploid_em_cluster.py` – Hard EM for Polyploid Phasing

**Location:** `algorithms/polyploid_em_cluster.py`

Phases polyploid reads using general hard EM clustering for K haplotypes.

#### **Usage**

```bash
python algorithms/polyploid_em_cluster.py -i reads.npz -k PLOIDY -o output_prefix [--max-iters K] [--tol-iters T] [-s SEED]
```

#### **Example**

```bash
python algorithms/polyploid_em_cluster.py -i demo_reads.npz -k 3 -o output/em3
```

#### **Parameters**

| Flag                    | Type | Default    | Description                    |
| ----------------------- | ---- | ---------- | ------------------------------ |
| `-i`, `--input`         | str  | *Required* | Input `.npz` reads file        |
| `-k`, `--ploidy`        | int  | *Required* | Number of haplotypes (K)       |
| `-o`, `--output-prefix` | str  | *Required* | Prefix for output files        |
| `--max-iters`           | int  | 40         | Maximum number of iterations   |
| `--tol-iters`           | int  | 3          | Tolerance for convergence      |
| `-s`, `--seed`          | int  | None       | Random seed for initialization |

#### **Outputs**

* `<prefix>.haplotypes.tsv`
* `<prefix>.assignments.tsv`
* `<prefix>.summary.json`

### 📊 `accuracy.py` – Evaluate Phasing Accuracy

**Location:** `benchmarking/accuracy.py`

Compares predicted haplotypes to ground truth, computes best label alignment, and reports overall accuracy.

#### **Usage**

```bash
python benchmarking/accuracy.py --truth truth.haplotypes.tsv --pred predicted.haplotypes.tsv [--output accuracy.json]
```

#### **Example**

```bash
python benchmarking/accuracy.py \
  --truth dataset/1/sim.haplotypes.tsv \
  --pred results/em_phase_1.haplotypes.tsv \
  --output results/em_phase_1.accuracy.json
```

#### **Parameters**

| Flag       | Type | Required   | Description                                    |
| ---------- | ---- | ---------- | ---------------------------------------------- |
| `--truth`  | str  | ✅ Yes      | Path to ground truth haplotypes `.tsv` file    |
| `--pred`   | str  | ✅ Yes      | Path to predicted haplotypes `.tsv` file       |
| `--output` | str  | ❌ Optional | Path to save evaluation results in JSON format |

#### **Output**

* **Console:** Prints best possible accuracy (after matching labels) and the optimal permutation of haplotypes.
* **Optional JSON:** If `--output` is specified, writes a report like:

```json
{
  "accuracy": 0.975,
  "label_permutation": [1, 0],
  "truth_shape": [2, 1000]
}
```

## 🧪 Benchmarking Pipeline

The benchmarking tool automates simulation, execution, and evaluation of phasing algorithms under different conditions. It enables reproducible experiments to assess how accuracy and runtime vary with input parameters like number of reads, error rate, or ploidy.

### 📌 Script

```bash
benchmarking/benchmark_runner.py
```

### ✅ Features

* Automatically selects the appropriate simulation script for diploid vs polyploid
* Supports multiple runs with different seeds
* Collects runtime, accuracy, and algorithm metadata
* Supports parameter sweeping (e.g. vary number of reads or error rate)
* Outputs results in a structured JSON file (`benchmark_summary.json`)

---

### 🔧 Usage

```bash
python benchmarking/benchmark_runner.py \
  --algorithms diploid-em diploid-mst \
  --ploidy 2 \
  --num-reads 5000 \
  --error-rate 0.01 \
  --num-runs 3 \
  --outdir results/benchmark_test
```

This will:

* Simulate 3 diploid datasets with 5000 reads and 1% error
* Run both `diploid-em` and `diploid-mst` on each
* Save outputs and accuracy reports
* Aggregate results into `results/benchmark_test/benchmark_summary.json`

---

### 🔁 Parameter Sweeping

To benchmark an algorithm under varying conditions:

```bash
python benchmarking/benchmark_runner.py \
  --algorithms diploid-em \
  --ploidy 2 \
  --error-rate 0.01 \
  --num-runs 3 \
  --vary num_reads \
  --vary-values 1000 5000 10000 20000 \
  --outdir results/sweep_reads
```

This will:

* Sweep over number of reads (1000, 5000, etc.)
* For each value, simulate and benchmark 3 runs
* Run `diploid-em` and evaluate accuracy

---

### 📁 Output

Each run produces:

* `.reads.sparse.tsv` + `.reads.npz` (simulated reads)
* `.haplotypes.tsv` (ground truth + predicted)
* `.assignments.tsv` (read-to-haplotype map)
* `.accuracy.json` (accuracy + label permutation)
* `benchmark_summary.json` (full aggregated log for plotting)


## 🚀 Unified Phasing CLI

The repository provides a unified command-line interface for running all supported haplotype phasing algorithms. Each algorithm is accessible via a subcommand.

### 📄 Script

```bash
python -m algorithms.cli.phase <algorithm> [options]
```

---

### 🧩 Supported Algorithms

| Subcommand           | Description                                 |
| -------------------- | ------------------------------------------- |
| `diploid-em`         | Diploid phasing using hard EM               |
| `diploid-mst`        | Diploid phasing using MST heuristic         |
| `polyploid-em`       | Polyploid phasing using hard EM             |
| `polyploid-spectral` | Polyploid phasing using spectral clustering |

---

### 🧪 Examples

**Diploid EM**

```bash
python -m algorithms.cli.phase diploid-em \
  -i data/sample.reads.npz \
  -o results/sample_em \
  --max-iters 30 \
  --tol-iters 2 \
  --seed 42
```

**Polyploid EM**

```bash
python -m algorithms.cli.phase polyploid-em \
  -i data/sample.reads.npz \
  -k 4 \
  -o results/poly_em \
  --max-iters 40 \
  --tol-iters 3
```

**Diploid MST**

```bash
python -m algorithms.cli.phase diploid-mst \
  -i data/sample.reads.npz \
  -o results/sample_mst \
  --min-overlap 3
```

**Polyploid Spectral**

```bash
python -m algorithms.cli.phase polyploid-spectral \
  -i data/sample.reads.npz \
  -k 4 \
  -o results/poly_spectral \
  --min-overlap 5
```

---

### 🧵 Common Parameters

| Flag            | Description                                  |
| --------------- | -------------------------------------------- |
| `-i, --input`   | Input `.npz` file (required)                 |
| `-o, --output`  | Output prefix (required)                     |
| `-k, --ploidy`  | Ploidy (required for polyploid algorithms)   |
| `--max-iters`   | Max EM iterations (default varies by method) |
| `--tol-iters`   | Early stopping threshold                     |
| `--min-overlap` | Minimum variant overlap to join fragments    |
| `--seed`        | Random seed for reproducibility              |
