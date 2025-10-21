# COMP4801
source .venv/bin/activate

Certainly! Below is the full **Section 3: Scripts Documentation** written in clean **Markdown format**, suitable for direct pasting into your `README.md`. It includes all scripts from the `algorithms/` and `dataset/` directories, with detailed descriptions, parameters, usage examples, and suggested values.

---

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

---

Let me know if you'd like this content as a `.md` file download or want to auto-generate usage info from your argparse definitions.
