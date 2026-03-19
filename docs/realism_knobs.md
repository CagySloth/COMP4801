## Realism knobs (long-read track)

This repo intentionally separates **pipeline correctness** from **data realism**.  
Realism knobs let you stress different failure modes while keeping the pipeline reproducible.

All knobs are exposed through `benchmark.longread_pipeline_runner`.

---

# Reference realism (Step 1)

## `--ref-preset`
Controls reference complexity:

- `plain`: random DNA (good baseline)
- `toy`: small special patterns (debug)
- `realistic`: adds regions to mimic repeats/low-complexity

The reference meta file (`<prefix>.ref.meta.json`) stores:
- `regions` (e.g., low-complexity / homopolymer regions)
- `duplications` (see below)

## Segment duplications (repeat-like)
- `--dup-segments K`
- `--dup-len LEN`
- `--dup-min-gap GAP`

Adds duplicated segments into the reference to create ambiguous mapping / harder calling/phasing.

Expected effects:
- called recall ↓ (some variants become harder to call)
- called switch error ↑ (wrong constraints / mapping ambiguity)
- phase sets ↑ (fragmentation)

---

# Truth realism (Step 2)

## SNPs
- `--num-snps`, `--het-rate`, `--min-distance`
- optional: `--avoid-regions --ref-meta <prefix>.ref.meta.json`

## Indels (small)
- `--num-indels`, `--indel-min-len`, `--indel-max-len`, `--indel-het-rate`

Indels stress alignment and calling, and may change SNP callability.  
**Recommended mode**: evaluate SNP phasing with `--phase-snps-only` (see below).

---

# Read simulation realism (Step 3)

## Length distribution
- `--len-model uniform` (default)  
- `--len-model lognormal --ln-mean 8.3 --ln-sigma 0.6` (more realistic tail)

Expected effects:
- longer reads → more connectivity → fewer phase sets

## Start distribution / dropout
- `--start-model uniform` (default)
- `--start-model dropout --dropout-fraction F --dropout-block-len B`

Creates “coverage deserts” (mimics mappability issues, sample dropout, etc.).
Expected effects:
- fewer shared het sites
- phase sets ↑ (fragmentation)

## Correlated error bursts (ONT)
- `--burst-prob P`
- `--burst-count C`
- `--burst-len L`
- `--burst-mult M`

Local, correlated high-error segments within reads.
Expected effects:
- calling recall ↓
- shared het recall ↓
- end-to-end effective phased recall can shift depending on which sites survive calling

## Platform & base error profile
- `--platform ont` (default) with `--ont-profile q20` or `classic`
- (Optional perfect mode: `--platform perfect`)

---

# Phasing/evaluation mode under indels

When indels are enabled, use:

```bash
python -m benchmark.longread_pipeline_runner ... --phase-snps-only
```

This:
- filters phasing input VCF to biallelic SNPs (`bcftools view -v snps -m2 -M2`)
- filters truth evaluation VCF to biallelic SNPs

This keeps the evaluation meaningful while indels still stress mapping/calling.

---

# Suggested “hard realistic ONT” setting (starter)

```bash
python -m benchmark.longread_pipeline_runner \
  --prefix output/acc_hard \
  --seed 0 \
  --ref-length 120000 --num-snps 1200 \
  --ref-preset realistic \
  --dup-segments 5 --dup-len 3000 --dup-min-gap 500 \
  --start-model dropout --dropout-fraction 0.1 --dropout-block-len 1000 \
  --burst-prob 0.6 --burst-count 2 --burst-len 300 --burst-mult 8 \
  --num-indels 120 --indel-min-len 1 --indel-max-len 5 --indel-het-rate 0.5 \
  --num-reads 300 --min-len 2000 --max-len 6000 \
  --platform ont --ont-profile q20 \
  --vcf-source both \
  --phase-snps-only
```
