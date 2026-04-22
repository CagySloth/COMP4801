# Algorithms and phasing backends

All phasing algorithms are exposed through:

```bash
python -m algorithms.cli.phase <subcommand> ...
```

---

## Diploid backends

### `diploid-em`
Diploid phasing on matrix-track inputs using a hard-EM style approach.

```bash
python -m algorithms.cli.phase diploid-em \
  -i input.reads.npz \
  --output-prefix output/dip_em
```

### `diploid-mst`
Graph-based diploid phasing on matrix-track inputs.

```bash
python -m algorithms.cli.phase diploid-mst \
  -i input.reads.npz \
  --output-prefix output/dip_mst
```

### `diploid-whats`
Vendored WhatsHap adapter for matrix / VCF-mode phasing.

```bash
python -m algorithms.cli.phase diploid-whats \
  -i input.reads.npz \
  --vcf input.vcf \
  --output-prefix output/dip_whats \
  --solver whatshap
```

#### Notes
- If `--vcf` is provided, only heterozygous sites are phased and homozygous sites are fixed by genotype.
- Supports `--solver {whatshap,hapchat}`.
- Can write `<prefix>.phased.vcf`.

### `diploid-whats-bam`
The main long-read phasing path in this project.

```bash
python -m algorithms.cli.phase diploid-whats-bam \
  --bam input.bam \
  --vcf input.vcf.gz \
  --output-prefix output/ws \
  --output-vcf output/ws.phased.vcf \
  --max-coverage 15 \
  --min-mapq 20 \
  --min-baseq 20
```

#### Inputs
- sorted, indexed BAM
- called or oracle VCF/VCF.GZ

#### Outputs
- phased VCF
- summary JSON

#### Main parameters
- `--max-coverage`
- `--min-mapq`
- `--min-baseq`
- `--solver`
- `--recomb-rate`
- `--sample`

This backend uses the vendored WhatsHap core under `vendor/whatshap_core`.

---

## Polyploid backends

### `polyploid-em`
Polyploid phasing on matrix-track inputs.

### `polyploid-spectral`
Polyploid spectral-clustering baseline on matrix-track inputs.

These are legacy / supporting parts of the repository and are not the main focus of the long-read benchmarking work.

---

## Supporting components

### `algorithms.diploid.whatshap_adapter`
Adapter that converts the project’s internal read representation into a WhatsHap-compatible `ReadSet` for the matrix / VCF-mode path.

### `algorithms.diploid.whatshap_bam_driver`
Long-read BAM + VCF driver used by the end-to-end pipeline.

### `algorithms.vendor.whatshap_vendor`
Helpers for importing the vendored WhatsHap implementation consistently.

---

## Practical recommendation

For the final-year project’s main experiments, the most relevant phasing backend is:

- `diploid-whats-bam`

The matrix-track backends remain useful for debugging, algorithm tests, and compatibility checks.
