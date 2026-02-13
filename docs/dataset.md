## 📂 `dataset` — Synthetic Data Generators

This repo provides **two different dataset generators**:

1) **Matrix generator** (`dataset.simulate`)  
   Outputs the variant-level representation used by most algorithms in this repo: truth haplotypes TSV + reads in NPZ/TSV.

2) **Long-read generator** (`dataset.longread.*`)  
   Outputs a WhatsHap-like workflow: reference FASTA + diploid truth (phased VCF + haplotypes FASTA) + reads FASTQ.

---

# 1) Matrix generator: `python -m dataset.simulate`

Generate a synthetic dataset (diploid or polyploid) and write it to disk.

## Parameters

- `-p, --ploidy` (int, required): number of haplotypes (2 = diploid)
- `-n, --num-variants` (int, required): number of variant sites
- `-r, --num-reads` (int, required): number of reads/fragments
- `-l, --read-length` (int, required): length of each read (in variant sites, not bp)
- `-e, --error-rate` (float): flip probability for an observed allele (default `0.01`)
- `-m, --missing-rate` (float): probability that an allele is missing (`-1`) (default `0.0`)
- `--maf-alpha`, `--maf-beta` (float): Beta distribution parameters for minor allele frequency (default `1.0, 1.0`)
- `--allow-monomorphic`: allow sites that are always 0 or always 1
- `--seed` (int): random seed for reproducibility
- `-o, --output-prefix` (str, required): output prefix (directory will be created if needed)

## Output files

For all ploidies:

- `{prefix}.haplotypes.tsv` — ground truth haplotypes
- `{prefix}.reads.sparse.tsv` — sparse read fragments (`read_id \t idx:allele \t idx:allele ...`)
- `{prefix}.reads.npz` — dense reads matrix for algorithm drivers

For **diploid only** (`--ploidy 2`), additionally:

- `{prefix}.vcf` — minimal, single-sample VCF with **unphased GT** derived from truth haplotypes
  - Homozygous sites: `0/0` or `1/1`
  - Heterozygous sites: `0/1`

## Example

```bash
python -m dataset.simulate \
  -p 2 -n 100 -r 50 -l 30 \
  -e 0.01 -m 0.0 \
  --seed 0 \
  -o output/demo
```

---

# 2) Long-read generator modules

These are building blocks for an end-to-end “FASTA → FASTQ → BAM → VCF → phase” workflow.

## 2.1 `python -m dataset.longread.reference`

Generates a reference FASTA:

- `{prefix}.ref.fasta`
- `{prefix}.ref.meta.json` (region metadata for “avoid regions” logic)

Example:

```bash
python -m dataset.longread.reference \
  -o output/lr_demo \
  --length 20000 \
  --seed 0
```

## 2.2 `python -m dataset.longread.truth`

Generates a *diploid truth* on top of the reference:

- `{prefix}.truth.vcf` (optionally phased)
- `{prefix}.hap1.fasta`
- `{prefix}.hap2.fasta`

Example:

```bash
python -m dataset.longread.truth \
  --ref output/lr_demo.ref.fasta \
  -o output/lr_demo \
  --seed 0 \
  --num-snps 200 \
  --het-rate 0.8 \
  --phased-truth \
  --random-phase
```

Optional “avoid certain reference regions”:

```bash
python -m dataset.longread.truth \
  --ref output/lr_demo.ref.fasta \
  --ref-meta output/lr_demo.ref.meta.json \
  --avoid-regions \
  -o output/lr_demo \
  --seed 0 \
  --num-snps 200 \
  --het-rate 0.8 \
  --phased-truth \
  --random-phase
```

## 2.3 `python -m dataset.longread.readsim`

Simulates long reads from the two haplotype FASTAs:

- `{prefix}.reads.fastq`
- `{prefix}.reads.truth.tsv` (which haplotype each read came from)

Example:

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

---

# Notes on formats

The matrix `.reads.npz` stores:

- `reads`: an `R x N` dense matrix with values `{0, 1, -1}` (`-1` = missing)
- `positions`: an `R x N` matrix of positions
  - For simulated matrix data this is usually a tiled `0..N-1` index grid.
  - Some drivers (e.g., WhatsHap VCF-mode) override it to preserve original coordinates after subsetting.
