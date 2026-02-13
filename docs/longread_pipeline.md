## Long-read end-to-end pipeline (Steps 1–7)

Goal: mimic a practical WhatsHap-like workflow **from VCF input to phased output**, but with full control for benchmarking.

This pipeline uses:
- reference FASTA
- simulated haplotypes → truth phased VCF
- reads FASTQ
- minimap2/samtools alignment → BAM
- bcftools calling → called VCF
- `diploid-whats-bam` → phased VCF
- `vcf_phase_eval` → evaluation JSON

System requirements:
- `minimap2`
- `samtools`
- `bcftools`

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
  --max-len 6000
```

---

# Option B: run steps manually

## Step 1: reference
```bash
python -m dataset.longread.reference -o output/lr_demo --length 20000 --seed 0
```

## Step 2: truth (phased VCF + hap FASTAs)
```bash
python -m dataset.longread.truth \
  --ref output/lr_demo.ref.fasta \
  -o output/lr_demo \
  --seed 0 \
  --num-snps 200 \
  --het-rate 0.8 \
  --phased-truth \
  --random-phase

bcftools view -Oz -o output/lr_demo.truth.vcf.gz output/lr_demo.truth.vcf
bcftools index -t output/lr_demo.truth.vcf.gz
```

## Step 3: reads FASTQ
```bash
python -m dataset.longread.readsim \
  --hap1 output/lr_demo.hap1.fasta \
  --hap2 output/lr_demo.hap2.fasta \
  -o output/lr_demo \
  --seed 0 \
  --num-reads 200 \
  --min-len 2000 \
  --max-len 6000
```

## Step 4: align reads → BAM
```bash
minimap2 -a -x map-ont output/lr_demo.ref.fasta output/lr_demo.reads.fastq \
  | samtools sort -o output/lr_demo.bam
samtools index output/lr_demo.bam
```

## Step 5: call variants → called VCF
```bash
bcftools mpileup -Ou -f output/lr_demo.ref.fasta output/lr_demo.bam \
  | bcftools call -mv -Oz -o output/lr_demo.called.vcf.gz
bcftools index -t output/lr_demo.called.vcf.gz
```

## Step 6: phase called VCF using BAM
```bash
python -m algorithms.cli.phase diploid-whats-bam \
  --bam output/lr_demo.bam \
  --vcf output/lr_demo.called.vcf.gz \
  --output-prefix output/lr_demo.ws \
  --output-vcf output/lr_demo.ws.phased.vcf \
  --max-coverage 15 --min-mapq 20 --min-baseq 20
```

## Step 7: evaluate phased VCF vs truth VCF
```bash
python -m benchmark.vcf_phase_eval \
  --truth output/lr_demo.truth.vcf.gz \
  --pred  output/lr_demo.ws.phased.vcf \
  --out   output/lr_demo.ws.eval.json
```

---

# Notes / current limitations

- The current read simulator generates “perfect” sequences (no explicit substitution/indel error model yet).
- `vcf_phase_eval` currently targets biallelic SNPs.
