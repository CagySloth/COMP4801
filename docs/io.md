# I/O and artifact contracts

The project relies heavily on stable, prefix-based output contracts.

This is especially important in the long-read pipeline because every experiment is traced through explicit files rather than only logs.

---

## Matrix-track artifacts

### Truth haplotypes TSV
`<prefix>.haplotypes.tsv`

Plain-text truth haplotypes used by TSV-based scoring.

### Dense read matrix NPZ
`<prefix>.reads.npz`

Typical contents:
- allele matrix
- optional position metadata

### Sparse read TSV
`<prefix>.reads.sparse.tsv`

Sparse fragment-style representation.

### Compatibility VCF
`<prefix>.vcf`

Unphased diploid VCF used by the matrix / VCF-mode WhatsHap adapter.

---

## Long-read pipeline artifacts

### Reference
- `<prefix>.ref.fasta`
- `<prefix>.ref.meta.json`

The metadata records information such as duplication events and region annotations.

### Truth and oracle
- `<prefix>.truth.vcf(.gz)`
- `<prefix>.truth.meta.json`
- `<prefix>.oracle.vcf(.gz)`
- `<prefix>.hap1.fasta`
- `<prefix>.hap2.fasta`

### Reads
- `<prefix>.reads.fastq`
- `<prefix>.reads.truth.tsv`
- `<prefix>.reads.meta.json`

### Alignment
- `<prefix>.bam`
- `<prefix>.bam.bai`

### Variant calling
- `<prefix>.called.vcf.gz`
- `<prefix>.called.vcf.gz.tbi`

### Phasing
- `<prefix>.ws*.phased.vcf`
- `<prefix>.ws*.summary.json`
- `<prefix>.ws*.eval.json`

### Run summary
- `<prefix>.pipeline.json`

This is the most important run-level artifact for downstream aggregation.

---

## Aggregated experiment artifacts

Within experiment directories, the standard high-level outputs are:
- `aggregate.csv`
- `plots/`

This is the interface used for report figures and result summaries.

---

## Why prefix-based naming matters

The project intentionally groups all outputs from one run under one prefix.

Benefits:
- easy traceability
- easy cleanup
- easier aggregation
- simpler debugging
- clear linkage between final plots and raw run outputs

---

## SNP-only filtering

When indels are enabled, the runner can create SNP-only phasing/evaluation inputs.

Equivalent manual filter:

```bash
bcftools view -v snps -m2 -M2 -Oz -o out.snps.vcf.gz in.vcf.gz
bcftools index -t out.snps.vcf.gz
```

This is the recommended way to preserve meaningful SNP phasing evaluation when indels are acting only as realism stressors.
