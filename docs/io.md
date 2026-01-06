## 📂 `algorithms/io` — Input/Output utilities

This package standardizes how simulated data and algorithm outputs are stored and loaded.

---

## Read formats

### Sparse TSV (`*.reads.sparse.tsv`)

One read per line:

```text
<read_id>\t<variant_index>:<allele>\t<variant_index>:<allele>...
```

Example:

```text
0    5:0    7:1    9:1
1    6:1    7:1
```

- Alleles are typically `0/1`
- Missing observations are simply absent from the row

### Dense NPZ (`*.reads.npz`)

A compressed numpy archive containing:

- `reads`: `R x N` matrix with values `{0, 1, -1}` (`-1` = missing)
- `positions`: `R x N` matrix of positions
  - For simulated data this is usually a tiled `0..N-1`.
  - Some drivers may override it (e.g., WhatsHap VCF-mode preserves original indices after subsetting).

The helper class `ReadsData` wraps these matrices.

---

## Key modules

- `parser.py`
  - `parse_sparse_tsv(path)` → list of sparse fragments
  - `parse_dense_tsv(path)` → dense matrix
  - `load_reads(path)` → `ReadsData` from `.npz` or `.tsv`
- `reads_data.py`
  - `ReadsData.from_fragments(fragments, num_variants=...)`
  - `ReadsData.from_npz(path)`
  - `ReadsData.to_npz(path)`
  - `ReadsData.to_sparse_tsv(path)`
- `writer.py`
  - `write_haplotypes_tsv(path, haplotypes)`
  - `write_assignments_tsv(path, assignments)`
  - `write_summary_json(path, summary)`

---

## Haplotype TSV (`*.haplotypes.tsv`)

Two common shapes:

- Diploid truth/prediction: 2 lines, each a `0/1` string of length `N`
- Polyploid: `k` lines

The benchmark scorer expects dense `0/1` strings for each haplotype line.
