## 📂 `algorithms` – Phasing Algorithms

This module implements multiple haplotype phasing strategies for diploid and polyploid data.

---

### 📁 `algorithms/diploid/`

Implements phasing methods designed specifically for diploid organisms (ploidy = 2).

#### 📄 `em.py` — EM Algorithm for Diploid Phasing

An Expectation-Maximization method assuming two haplotypes and probabilistic assignment of reads.

```bash
python -m algorithms.cli.phase diploid_em -i input.npz -o result_prefix
```

#### 📄 `mst.py` — Minimum Spanning Tree (MST) Phasing

Graph-based approach using co-occurrence of alleles to build a minimum spanning tree.

```bash
python -m algorithms.cli.phase diploid_mst -i input.npz -o result_prefix
```

---

### 📁 `algorithms/polyploid/`

Phasing strategies for data with ≥3 haplotypes (polyploid genomes).

#### 📄 `em.py` — EM Algorithm for Polyploid Phasing

Extended EM method supporting arbitrary ploidy.

```bash
python -m algorithms.cli.phase polyploid_em -i input.npz -o result_prefix -k 4
```

#### 📄 `spectral.py` — Spectral Clustering for Polyploid Phasing

Clustering-based phasing using the similarity graph of reads and spectral decomposition.

```bash
python -m algorithms.cli.phase polyploid_spectral -i input.npz -o result_prefix -k 4
```

---

### 📁 `algorithms/io/`

Handles I/O utilities like loading `.tsv` and `.npz` files, writing results, and standardizing data structures.

| File            | Responsibility                                      |
| --------------- | --------------------------------------------------- |
| `parser.py`     | Parses reads from dense/sparse `.tsv` or `.npz`     |
| `writer.py`     | Outputs `.haplotypes.tsv`, `.assignments.tsv` files |
| `reads_data.py` | Defines the `ReadsData` structure for algorithms    |

These modules are used internally by CLI scripts and algorithm runners.