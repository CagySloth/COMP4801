## 📂 `algorithms` — Phasing Algorithms

This package contains diploid and polyploid phasing implementations and their shared adapters/utilities.

Run algorithms via the CLI:

```bash
python -m algorithms.cli.phase <subcommand> ...
```

---

## Diploid algorithms (`algorithms/diploid`)

### `diploid-em` — Expectation Maximization

```bash
python -m algorithms.cli.phase diploid-em \
  -i input.reads.npz \
  --output-prefix output/dip_em
```

### `diploid-mst` — Minimum Spanning Tree (graph-based)

```bash
python -m algorithms.cli.phase diploid-mst \
  -i input.reads.npz \
  --output-prefix output/dip_mst
```

### `diploid-whats` — Vendored WhatsHap core (matrix / VCF-mode)

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
python -m algorithms.cli.phase diploid-whats \
  -i output/demo.reads.npz \
  --vcf output/demo.vcf \
  --output-prefix output/demo_phased \
  --solver whatshap
```

### `diploid-whats-bam` — WhatsHap-like BAM + VCF phasing

This mode is designed to mimic a more practical WhatsHap workflow:

- Input: **BAM** (aligned reads) + **VCF** (called variants with GT)
- Internally extracts allele observations from BAM at VCF sites, then runs:
  - read selection
  - solver (vendored WhatsHap core)
- Output: phased VCF + summary JSON

Example:

```bash
python -m algorithms.cli.phase diploid-whats-bam \
  --bam output/lr_demo.bam \
  --vcf output/lr_demo.called.vcf.gz \
  --output-prefix output/lr_demo.ws \
  --output-vcf output/lr_demo.ws.phased.vcf \
  --max-coverage 15 \
  --min-mapq 20 \
  --min-baseq 20
```

---

## Polyploid algorithms (`algorithms/polyploid`)

### `polyploid-em`

```bash
python -m algorithms.cli.phase polyploid-em \
  -i input.reads.npz \
  --ploidy 4 \
  --output-prefix output/poly_em
```

### `polyploid-spectral`

```bash
python -m algorithms.cli.phase polyploid-spectral \
  -i input.reads.npz \
  --ploidy 4 \
  --output-prefix output/poly_spec
```
