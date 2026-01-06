## 📂 `algorithms` — Phasing Algorithms

This package contains diploid and polyploid phasing implementations and their shared adapters/utilities.

The easiest way to run algorithms is via the CLI:

```bash
python -m algorithms.cli.phase <subcommand> ...
```

---

## Diploid algorithms (`algorithms/diploid`)

### `diploid-em` — Expectation Maximization

- Probabilistic assignment of reads to two haplotypes.
- Useful baseline for accuracy vs. runtime tradeoffs.

```bash
python -m algorithms.cli.phase diploid-em       -i input.reads.npz       --output-prefix output/dip_em
```

### `diploid-mst` — Minimum Spanning Tree (graph-based)

- Builds a variant graph from co-occurrence on reads and derives a phasing via spanning-tree logic.

```bash
python -m algorithms.cli.phase diploid-mst       -i input.reads.npz       --output-prefix output/dip_mst
```

### `diploid-whats` — Vendored WhatsHap core

This integrates a vendored WhatsHap core (`vendor/whatshap_core`) and supports two modes:

**Matrix-mode (no VCF)**  
- Input: only `reads.npz`
- Solver: uses `HapChatCore` on the selected reads
- Output: TSV haplotypes + assignments + JSON summary

**VCF-mode (`--vcf ...`)**  
- Input: `reads.npz` + a VCF containing unphased GTs
- Only **heterozygous** variants (from GT) are phased; homozygous sites are fixed by GT.
- Solvers:
  - `--solver whatshap` (default) → **PedigreeDPTable** (WhatsHap default DP solver)
  - `--solver hapchat` → **HapChatCore** (MEC-style solver)
- Writes a phased VCF: `{output-prefix}.phased.vcf`

Example:

```bash
python -m algorithms.cli.phase diploid-whats       -i output/demo.reads.npz       --vcf output/demo.vcf       --output-prefix output/demo_phased       --solver whatshap
```

WhatsHap-style details implemented in this repo:

- **Read selection** via `whatshap.readselect.readselection(max_coverage=...)`
- **Phase sets (PS)** computed from connected components of the selected read graph  
  (PS is reported as **leftmost variant index + 1**).

---

## Polyploid algorithms (`algorithms/polyploid`)

### `polyploid-em`

- EM-style baseline extended to `k` haplotypes.

```bash
python -m algorithms.cli.phase polyploid-em       -i input.reads.npz       --ploidy 4       --output-prefix output/poly_em
```

### `polyploid-spectral`

- Spectral / matrix-factorization inspired approach for polyploid phasing.

```bash
python -m algorithms.cli.phase polyploid-spectral       -i input.reads.npz       --ploidy 4       --output-prefix output/poly_spec
```
