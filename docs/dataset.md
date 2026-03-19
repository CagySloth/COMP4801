## 📂 `dataset` — Synthetic Data Generators

This repo provides two dataset generators:

1) **Matrix generator** (`dataset.simulate`)  
   Outputs the variant-level representation used by most algorithms in this repo: truth haplotypes TSV + reads in NPZ/TSV.

2) **Long-read generator** (`dataset.longread.*`)  
   Outputs a WhatsHap-like workflow: reference FASTA → truth VCF + haplotypes FASTA → reads FASTQ.

---

# 1) Matrix generator: `python -m dataset.simulate`

Generate a synthetic dataset (diploid or polyploid) and write it to disk.

Common parameters:
- `-p, --ploidy`: number of haplotypes (2=diploid)
- `-n, --num-variants`: number of variant sites
- `-r, --num-reads`: number of reads/fragments
- `-l, --read-length`: mean length (in variants, not bp)
- `-e, --error-rate`: allele flips
- `-m, --missing-rate`: missing alleles
- `--seed`: reproducibility
- `-o, --output-prefix`: output prefix

Outputs:
- `<prefix>.haplotypes.tsv` (truth)
- `<prefix>.reads.npz` (dense reads×variants allele matrix)
- `<prefix>.reads.sparse.tsv` (sparse fragments)

Diploid-only extra:
- `<prefix>.vcf` (unphased GT; used for VCF-mode phasing tests)

---

# 2) Long-read generator (Steps 1–3)

## Step 1: `dataset.longread.reference`

```bash
python -m dataset.longread.reference \
  -o output/demo \
  --length 80000 \
  --seed 0 \
  --preset plain \
  --dup-segments 0
```

Key options:
- `--preset {plain,toy,realistic}`
- `--dup-segments/--dup-len/--dup-min-gap` (repeat-like duplications)

Outputs:
- `<prefix>.ref.fasta`
- `<prefix>.ref.meta.json` (regions + duplications)

## Step 2: `dataset.longread.truth`

```bash
python -m dataset.longread.truth \
  --ref output/demo.ref.fasta \
  -o output/demo \
  --seed 0 \
  --num-snps 800 --het-rate 0.8 \
  --phased-truth --random-phase
```

Indels:
```bash
python -m dataset.longread.truth \
  --ref output/demo.ref.fasta \
  -o output/demo \
  --seed 0 \
  --num-snps 800 \
  --num-indels 80 --indel-min-len 1 --indel-max-len 5 --indel-het-rate 0.5 \
  --phased-truth --random-phase
```

Region-avoidance:
- `--avoid-regions --ref-meta <prefix>.ref.meta.json`

Outputs:
- `<prefix>.truth.vcf` (+ `.truth.meta.json`)
- `<prefix>.hap1.fasta`, `<prefix>.hap2.fasta`

## Step 3: `dataset.longread.readsim`

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

Realism knobs:
- length: `--len-model lognormal --ln-mean 8.3 --ln-sigma 0.6`
- dropout: `--start-model dropout --dropout-fraction 0.1 --dropout-block-len 1000`
- bursts: `--burst-prob 0.6 --burst-count 2 --burst-len 300 --burst-mult 8`

Outputs:
- `<prefix>.reads.fastq`
- `<prefix>.reads.truth.tsv`
- `<prefix>.reads.meta.json`

See `docs/realism_knobs.md` for expected effects.
