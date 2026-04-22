# Long-read end-to-end pipeline

This is the main research workflow of the project.

The pipeline follows the practical structure:

**Reference FASTA → Truth / Oracle VCF → Reads FASTQ → Alignment BAM → Called VCF → Phased VCF → Evaluation reports**

It is orchestrated through:

```bash
python -m benchmark.longread_pipeline_runner
```

---

## Why this pipeline exists

The project is not only evaluating a phasing algorithm in isolation. It is evaluating how WhatsHap behaves inside a realistic long-read workflow, where losses may happen in:

- alignment
- variant calling
- phasing connectivity
- evidence quality

This is why the pipeline supports both:
- **oracle** phasing input
- **called** phasing input

Comparing them enables attribution of phaser-limited versus caller-limited behavior.

---

## Pipeline stages

### 1. Reference generation
Generates a synthetic reference FASTA and metadata.

Main output:
- `<prefix>.ref.fasta`
- `<prefix>.ref.meta.json`

### 2. Truth generation
Generates truth SNPs / indels, phased truth haplotypes, and truth metadata.

Main output:
- `<prefix>.truth.vcf(.gz)`
- `<prefix>.hap1.fasta`
- `<prefix>.hap2.fasta`
- `<prefix>.truth.meta.json`

### 3. Oracle VCF derivation
The pipeline runner derives an oracle VCF from the truth sites for phaser-limited benchmarking.

Main output:
- `<prefix>.oracle.vcf(.gz)`

### 4. Read simulation
Simulates ONT-like reads with configurable realism knobs.

Main output:
- `<prefix>.reads.fastq`
- `<prefix>.reads.truth.tsv`
- `<prefix>.reads.meta.json`

### 5. Alignment
Aligns reads to the reference.

External tools used:
- `minimap2`
- `samtools`

Main output:
- `<prefix>.bam`
- `<prefix>.bam.bai`

### 6. Variant calling
Calls variants from the BAM.

External tool used:
- `bcftools`

Main output:
- `<prefix>.called.vcf.gz`
- `<prefix>.called.vcf.gz.tbi`

### 7. Phasing
Runs vendored WhatsHap phasing on the oracle VCF, called VCF, or both.

Main output:
- `<prefix>.ws*.phased.vcf`
- `<prefix>.ws*.summary.json`

### 8. Evaluation
Evaluates phased output against truth.

Main output:
- `<prefix>.ws*.eval.json`
- `<prefix>.pipeline.json`

---

## Common runner parameters

### Global
- `--prefix`
- `--seed`

### Reference and truth
- `--ref-preset`
- `--ref-length`
- `--dup-segments`, `--dup-len`, `--dup-min-gap`
- `--num-snps`, `--het-rate`
- `--avoid-regions`
- `--num-indels`, `--indel-min-len`, `--indel-max-len`, `--indel-het-rate`

### Read simulation
- `--num-reads`
- `--min-len`, `--max-len`
- `--platform`
- `--ont-profile`
- `--len-model`, `--ln-mean`, `--ln-sigma`
- `--start-model`, `--dropout-fraction`, `--dropout-block-len`
- `--burst-prob`, `--burst-count`, `--burst-len`, `--burst-mult`

### Alignment and calling
- `--map-preset`
- `--call-min-mapq`
- `--call-min-baseq`

### Phasing and evaluation
- `--vcf-source {called,oracle,both}`
- `--max-coverage`
- `--min-mapq`
- `--min-baseq`
- `--phase-snps-only`
- `--eval-snps-only`

---

## Example run

```bash
python -m benchmark.longread_pipeline_runner \
  --prefix output/lr_demo \
  --seed 0 \
  --ref-length 80000 \
  --num-snps 800 \
  --het-rate 0.8 \
  --num-reads 200 \
  --min-len 2000 \
  --max-len 6000 \
  --platform ont \
  --ont-profile q20 \
  --vcf-source both
```

---

## Indel-enabled runs

If indels are introduced into truth, the project’s recommended mode is:

```bash
--phase-snps-only --eval-snps-only
```

This preserves meaningful SNP phasing evaluation while still allowing indels to stress the upstream pipeline.

The experiment driver automatically turns these on when `num_indels > 0`.

---

## Why prefix-based artifacts matter

All artifacts of a run share one prefix. This makes the workflow:
- reproducible
- traceable
- easy to aggregate
- easy to debug

For the final report, this is also what allows each result figure to be traced back to exact run-level JSON summaries.
