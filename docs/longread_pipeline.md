## Long-read end-to-end pipeline (Steps 1–7)

Goal: mimic a practical WhatsHap-like workflow **from BAM+VCF input to phased VCF output**, while keeping full control for benchmarking.

Pipeline stages:

1. **Reference**: `dataset.longread.reference` → `<prefix>.ref.fasta` (+ `<prefix>.ref.meta.json`)
2. **Truth**: `dataset.longread.truth` → `<prefix>.truth.vcf` (+ hap FASTAs)
3. **Reads**: `dataset.longread.readsim` → `<prefix>.reads.fastq` (+ read truth/meta)
4. **Alignment**: minimap2 + samtools → `<prefix>.bam`
5. **Variant calling**: bcftools mpileup/call → `<prefix>.called.vcf.gz`
6. **Phasing**: `diploid-whats-bam` → phased VCF + summary
7. **Evaluation**: `benchmark.vcf_phase_eval` → JSON metrics

External tools required:
- `minimap2`, `samtools`, `bcftools`

---

# Option A (recommended): run everything automatically

```bash
python -m benchmark.longread_pipeline_runner \
  --prefix output/lr_demo \
  --seed 0 \
  --ref-length 20000 \
  --num-snps 200 \
  --het-rate 0.8 \
  --num-reads 200 \
  --min-len 2000 \
  --max-len 6000 \
  --platform ont --ont-profile q20 \
  --vcf-source both
```

This writes `<prefix>.pipeline.json` which includes:
- paths to all generated artifacts
- callset precision/recall (truth vs called SNPs)
- phasing evaluation for `called` and/or `oracle`
- derived end-to-end metrics (effective phased recall, phasing rate, fragmentation)

---

# Option B: run steps manually

## Step 1: reference

```bash
python -m dataset.longread.reference \
  -o output/demo \
  --length 20000 \
  --seed 0 \
  --preset plain
```

Outputs:
- `output/demo.ref.fasta`
- `output/demo.ref.meta.json` (regions + duplications)

## Step 2: truth (+ optional indels)

```bash
python -m dataset.longread.truth \
  --ref output/demo.ref.fasta \
  -o output/demo \
  --seed 0 \
  --num-snps 200 \
  --het-rate 0.8 \
  --phased-truth --random-phase
```

If you want indels:
```bash
python -m dataset.longread.truth \
  --ref output/demo.ref.fasta \
  -o output/demo \
  --seed 0 \
  --num-snps 800 \
  --num-indels 80 --indel-min-len 1 --indel-max-len 5 --indel-het-rate 0.5 \
  --phased-truth --random-phase
```

Outputs:
- `output/demo.truth.vcf` (+ `.meta.json`)
- `output/demo.hap1.fasta`, `output/demo.hap2.fasta`

## Step 3: read simulation (ONT)

```bash
python -m dataset.longread.readsim \
  --hap1 output/demo.hap1.fasta \
  --hap2 output/demo.hap2.fasta \
  -o output/demo \
  --seed 0 \
  --num-reads 200 \
  --min-len 2000 --max-len 6000 \
  --platform ont --ont-profile q20
```

Optional realism knobs:
- `--len-model lognormal --ln-mean 8.3 --ln-sigma 0.6`
- `--start-model dropout --dropout-fraction 0.1 --dropout-block-len 1000`
- `--burst-prob 0.6 --burst-count 2 --burst-len 300 --burst-mult 8`

Outputs:
- `output/demo.reads.fastq`
- `output/demo.reads.truth.tsv`
- `output/demo.reads.meta.json`

## Step 4: alignment (BAM)

```bash
minimap2 -a -x map-ont output/demo.ref.fasta output/demo.reads.fastq \
  | samtools sort -o output/demo.bam
samtools index output/demo.bam
```

## Step 5: calling (called VCF)

```bash
bcftools mpileup -Ou -f output/demo.ref.fasta -q 20 -Q 15 output/demo.bam \
  | bcftools call -mv -Oz -o output/demo.called.vcf.gz
bcftools index -t output/demo.called.vcf.gz
```

## Step 6: phasing (WhatsHap core, BAM+VCF)

```bash
python -m algorithms.cli.phase diploid-whats-bam \
  --bam output/demo.bam \
  --vcf output/demo.called.vcf.gz \
  --output-prefix output/demo.ws \
  --output-vcf output/demo.ws.phased.vcf \
  --max-coverage 15 --min-mapq 20 --min-baseq 20
```

## Step 7: evaluation

```bash
python -m benchmark.vcf_phase_eval \
  --truth output/demo.truth.vcf.gz \
  --pred  output/demo.ws.phased.vcf \
  --out   output/demo.ws.eval.json
```

---

# Oracle vs Called runs

`benchmark.longread_pipeline_runner` supports `--vcf-source`:

- `called`: phase the **called VCF** → measures end-to-end effects (calling + phasing)
- `oracle`: phase an **oracle VCF** derived from truth (GT unphased, PS removed) → isolates phasing quality from calling
- `both`: run both (recommended)

---

# Indels: use SNP-only phasing/eval (current recommended mode)

When `--num-indels > 0`, use:

- `--phase-snps-only` on the pipeline runner, or manually filter:
  - `bcftools view -v snps -m2 -M2 ...`

Rationale: indels stress alignment/calling realism, but current evaluation focuses on **biallelic SNP phasing**.

---

# Metrics to report (from pipeline.json)

Common report fields:
- callset precision/recall (truth SNPs vs called SNPs)
- `phase_accuracy_blockflip` (block-flip aware)
- `switch_error_rate`
- `num_phase_sets` (fragmentation)
- derived:
  - `shared_het_recall`
  - `phasing_rate_on_shared_het`
  - `effective_phased_recall` (recommended end-to-end headline)

See `docs/benchmark.md` and `docs/validation.md` for sanity expectations.
