## 📂 `algorithms` — Phasing Algorithms

Run algorithms via:

```bash
python -m algorithms.cli.phase <subcommand> ...
```

---

## Diploid algorithms

### `diploid-em`
Input: `.reads.npz`  
Output: `<prefix>.haplotypes.tsv` + summary

```bash
python -m algorithms.cli.phase diploid-em \
  -i input.reads.npz \
  --output-prefix output/dip_em
```

### `diploid-mst`
Graph-based heuristic.

```bash
python -m algorithms.cli.phase diploid-mst \
  -i input.reads.npz \
  --output-prefix output/dip_mst
```

### `diploid-whats` (matrix / VCF-mode adapter)
Uses **vendored WhatsHap core** on NPZ matrices.

Modes:
- Matrix-only: phase all sites (mostly for compatibility testing)
- VCF-mode: if `--vcf` is provided, phase **heterozygous sites only** (WhatsHap-like)

```bash
python -m algorithms.cli.phase diploid-whats \
  -i input.reads.npz \
  --vcf input.vcf \
  --output-prefix output/dip_whats \
  --solver whatshap
```

Solvers:
- `--solver whatshap`: DP-style (PedigreeDPTable)
- `--solver hapchat`: MEC-style (HapChatCore)

Outputs:
- `<prefix>.haplotypes.tsv` (dense; for TSV scorer)
- `<prefix>.phased.vcf` (authoritative phasing output for VCF-mode)
- `<prefix>.summary.json`

### `diploid-whats-bam` (BAM + VCF → phased VCF)
This is the **WhatsHap-like long-read phasing** path.

Input:
- BAM (aligned reads)
- VCF (variants + genotypes)

Output:
- phased VCF + summary JSON

```bash
python -m algorithms.cli.phase diploid-whats-bam \
  --bam input.bam \
  --vcf input.vcf.gz \
  --output-prefix output/ws \
  --output-vcf output/ws.phased.vcf \
  --max-coverage 15 --min-mapq 20 --min-baseq 20
```

Notes:
- Only het genotypes are phased
- PS is assigned by connectivity across selected reads
- In indel experiments, prefer SNP-only VCF inputs (`--phase-snps-only` in the pipeline runner)

---

## Polyploid algorithms

### `polyploid-em`, `polyploid-spectral`
These operate on `.reads.npz` and output polyploid haplotypes TSV.

See `python -m algorithms.cli.phase --help` for details.
