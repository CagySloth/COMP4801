## 📂 `algorithms/io` — Input/Output utilities

This package standardizes how simulated data and algorithm outputs are stored and loaded for the **matrix track**.

---

## Matrix-track read formats

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
  - For simulated matrix data this is usually a tiled `0..N-1`.

The helper class `ReadsData` wraps these matrices.

---

## Matrix-track outputs

### Haplotype TSV (`*.haplotypes.tsv`)

Two common shapes:

- Diploid truth/prediction: 2 lines, each a `0/1` string of length `N`
- Polyploid: `k` lines

The benchmark scorer expects dense `0/1` strings for each haplotype line.

### Assignments TSV (`*.assignments.tsv`)

Algorithm-dependent read → haplotype assignments (mainly for debugging/analysis).

### Summary JSON (`*.summary.json`)

Algorithm metadata (runtime breakdowns, selected reads, solver name, etc.).

---

## Long-read pipeline I/O (not handled by `algorithms/io`)

The long-read track uses standard genomics formats:

- Reference: `*.ref.fasta`
- Haplotypes: `*.hap1.fasta`, `*.hap2.fasta`
- Reads: `*.reads.fastq`
- Alignments: `*.bam` (+ `.bai`)
- Variants: `*.truth.vcf(.gz)`, `*.called.vcf.gz`
- Phasing output: `*.phased.vcf`
