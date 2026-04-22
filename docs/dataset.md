# Dataset generation

This repository provides two dataset-generation tracks:

1. **Legacy matrix-track generation**
2. **Long-read generation for the end-to-end pipeline**

---

## 1. Matrix-track generator

Run:

```bash
python -m dataset.simulate \
  -p 2 -n 200 -r 200 -l 60 \
  -e 0.01 -m 0.0 \
  --seed 0 \
  -o output/demo
```

### Main parameters
- `-p, --ploidy` — number of haplotypes
- `-n, --num-variants` — number of variant sites
- `-r, --num-reads` — number of reads / fragments
- `-l, --read-length` — mean fragment length (in variant units)
- `-e, --error-rate` — allele-flip noise
- `-m, --missing-rate` — missing observations
- `--seed` — reproducibility
- `-o, --output-prefix` — output prefix

### Outputs
- `<prefix>.haplotypes.tsv`
- `<prefix>.reads.npz`
- `<prefix>.reads.sparse.tsv`
- `<prefix>.vcf` (diploid compatibility aid)

---

## 2. Long-read generation

The long-read track is split into three explicit generation stages.

### Step 1: synthetic reference

```bash
python -m dataset.longread.reference \
  -o output/demo \
  --length 80000 \
  --seed 0 \
  --preset plain
```

#### Main parameters
- `--length`
- `--contig`
- `--seed`
- `--preset {plain,toy,realistic}`
- `--dup-segments`, `--dup-len`, `--dup-min-gap`

#### Additional fine-grained realism controls
The reference generator also supports optional:
- homopolymer insertion
- STR insertion
- GC-window bias construction

These are available as direct CLI parameters and recorded in the output metadata.

#### Outputs
- `<prefix>.ref.fasta`
- `<prefix>.ref.meta.json`

---

### Step 2: truth generation

```bash
python -m dataset.longread.truth \
  --ref output/demo.ref.fasta \
  -o output/demo \
  --seed 0 \
  --num-snps 800 \
  --het-rate 0.8 \
  --phased-truth \
  --random-phase
```

#### Main parameters
- `--num-snps`
- `--het-rate`
- `--min-distance`
- `--ref-meta`
- `--avoid-regions`
- `--num-indels`
- `--indel-min-len`
- `--indel-max-len`
- `--indel-het-rate`
- `--sample`

#### Outputs
- `<prefix>.truth.vcf` (and compressed/indexed versions in the pipeline)
- `<prefix>.oracle.vcf` (created by the pipeline runner from truth)
- `<prefix>.hap1.fasta`
- `<prefix>.hap2.fasta`
- `<prefix>.truth.meta.json`

---

### Step 3: ONT-like read simulation

```bash
python -m dataset.longread.readsim \
  --hap1 output/demo.hap1.fasta \
  --hap2 output/demo.hap2.fasta \
  -o output/demo \
  --seed 0 \
  --num-reads 200 \
  --min-len 2000 \
  --max-len 6000 \
  --platform ont \
  --ont-profile q20
```

#### Main parameters
- `--num-reads`
- `--min-len`, `--max-len`
- `--platform {perfect,ont}`
- `--ont-profile {classic,q20}`
- `--len-model {uniform,lognormal}`
- `--ln-mean`, `--ln-sigma`
- `--start-model {uniform,dropout}`
- `--dropout-fraction`, `--dropout-block-len`
- `--burst-prob`, `--burst-count`, `--burst-len`, `--burst-mult`

#### Additional direct controls
If you run the simulator directly, you can also set:
- explicit substitution / insertion / deletion rates
- quality-model parameters
- homopolymer error amplification
- `--hap1-frac`

#### Outputs
- `<prefix>.reads.fastq`
- `<prefix>.reads.truth.tsv`
- `<prefix>.reads.meta.json`

---

## Why the generation is split into stages

The long-read generation path is intentionally stage-local:
- reference-level realism belongs in `reference.py`
- truth-level realism belongs in `truth.py`
- read-level realism belongs in `readsim.py`

This makes the pipeline easier to debug, easier to document, and easier to sweep in controlled experiments.
