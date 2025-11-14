## 📂 `algorithms/io` – Input/Output Utilities

This module standardizes file I/O, including reads/haplotypes parsing, sparse format support, and output writing. It also includes a `ReadsData` container for consistent downstream processing.

---

### 📄 `parser.py`

Handles reading `.tsv` (dense/sparse) and `.npz` files, converting them into a unified `ReadsData` object.

#### Functions

| Function                 | Description                                                                   |
| ------------------------ | ----------------------------------------------------------------------------- |
| `load_reads(path)`       | Load reads from `.npz` or `.tsv` (dense or sparse) into a `ReadsData` object. |
| `parse_sparse_tsv(path)` | Read a sparse TSV file into list of fragments (id, indices, values).          |
| `parse_dense_tsv(path)`  | Read a dense TSV into a NumPy matrix.                                         |
| `is_sparse_tsv(path)`    | Heuristic: detect if a `.tsv` file is sparse by checking for colons.          |

#### Example

```python
from algorithms.io.parser import load_reads

reads = load_reads("dataset/example.reads.npz")
print(reads.reads.shape, reads.positions.shape, reads.num_variants)
```

---

### 📄 `writer.py`

Functions to write outputs such as predicted haplotypes, read assignments, and result summaries.

#### Functions

| Function                                   | Description                                      |
| ------------------------------------------ | ------------------------------------------------ |
| `write_haplotypes_tsv(path, haplotypes)`   | Save haplotypes as TSV rows (variants × ploidy). |
| `write_assignments_tsv(path, assignments)` | Save inferred read → haplotype assignments.      |
| `write_summary(path, metrics: dict)`       | Save benchmark/evaluation metrics to JSON.       |

#### Example

```python
from algorithms.io.writer import write_haplotypes_tsv, write_assignments_tsv

write_haplotypes_tsv("results/hap.tsv", haplotypes)
write_assignments_tsv("results/assignments.tsv", assignments)
```

---

### 📄 `reads_data.py`

A dataclass container for normalized read datasets, compatible across I/O formats.

#### Class

```python
@dataclass
class ReadsData:
    reads: np.ndarray          # (num_reads, read_length)
    positions: np.ndarray      # (num_reads, read_length)
    num_variants: int

    @staticmethod
    def from_fragments(fragments: list[dict]) -> "ReadsData":
        ...
```

#### Example

```python
from algorithms.io.reads_data import ReadsData

reads_obj = ReadsData(reads=..., positions=..., num_variants=...)
```

Or build from sparse fragments:

```python
from algorithms.io.parser import parse_sparse_tsv
from algorithms.io.reads_data import ReadsData

fragments = parse_sparse_tsv("reads.sparse.tsv")
reads_obj = ReadsData.from_fragments(fragments)
```
