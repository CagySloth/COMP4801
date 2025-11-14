## 📂 `dataset` – Simulated Data Generator

This module contains the main script for generating synthetic read and haplotype data used in benchmarking.

---

### 📄 `simulate.py`

Generates random haplotypes and sequencing reads for diploid or polyploid organisms with configurable parameters.

#### Usage

```bash
python dataset/simulate.py [options]
```

#### Options

| Option                 | Type  | Default | Description                                                        |
| ---------------------- | ----- | ------- | ------------------------------------------------------------------ |
| `-p`, `--ploidy`       | int   | 2       | Number of haplotypes (e.g., 2 for diploid, ≥3 for polyploid)       |
| `-n`, `--num-variants` | int   | 1000    | Number of variant positions                                        |
| `-r`, `--num-reads`    | int   | 5000    | Number of sequencing reads                                         |
| `-l`, `--read-length`  | int   | 50      | Length of each read                                                |
| `-e`, `--error-rate`   | float | 0.01    | Per-base sequencing error rate                                     |
| `-m`, `--missing-rate` | float | 0.05    | Fraction of missing bases per read                                 |
| `--maf-alpha`          | float | 0.4     | Alpha parameter of Beta distribution for allele frequency sampling |
| `--maf-beta`           | float | 0.4     | Beta parameter of Beta distribution for allele frequency sampling  |
| `--allow-monomorphic`  | flag  | False   | If set, allows non-polymorphic sites (otherwise filters them out)  |
| `-s`, `--seed`         | int   | 42      | Random seed for reproducibility                                    |
| `-o`, `--output`       | str   | —       | Output path prefix (no extension)                                  |

#### Output Files

Given an output prefix like `sim_data`, it generates:

* `sim_data.haplotypes.tsv`: ground-truth haplotypes (ploidy × variants)
* `sim_data.reads.sparse.tsv`: sparse-format read matrix
* `sim_data.truth.assignments.tsv`: true read-to-haplotype mapping

#### Example

```bash
python dataset/simulate.py \
    -p 4 \
    -n 800 \
    -r 6000 \
    -l 40 \
    -e 0.02 \
    -m 0.05 \
    --maf-alpha 0.3 \
    --maf-beta 0.5 \
    --allow-monomorphic \
    -s 123 \
    -o dataset/sim_poly_example
```

This would generate a polyploid dataset with 4 haplotypes, 800 SNPs, and 6000 noisy reads.