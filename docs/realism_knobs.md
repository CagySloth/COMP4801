# Realism knobs

The project uses **realism knobs** as controlled stressors.

They are not intended to reproduce every aspect of biological sequencing perfectly. Instead, they introduce practical difficulty modes into a synthetic benchmark so that:

- failure modes become measurable
- attribution becomes possible
- optimization can be studied in a controlled way

---

## Main realism stressors in the project

### 1. Duplicated regions
Reference-level realism.

Main parameters:
- `--dup-segments`
- `--dup-len`
- `--dup-min-gap`

What it simulates:
- repeated / ambiguous genomic regions
- weaker mapping confidence
- noisier downstream evidence for calling and phasing

Expected effects:
- called recall may fall
- switch error may rise
- phasing completeness may degrade under stronger settings

---

### 2. Coverage dropout
Read-start / coverage realism.

Main parameters:
- `--start-model dropout`
- `--dropout-fraction`
- `--dropout-block-len`

What it simulates:
- local coverage deserts
- missing evidence between nearby heterozygous sites

Expected effects:
- more phase sets
- lower phasing completeness
- reduced called overlap under stronger dropout

---

### 3. Correlated error bursts
Read-level structured noise.

Main parameters:
- `--burst-prob`
- `--burst-count`
- `--burst-len`
- `--burst-mult`

What it simulates:
- clusters of local sequencing errors rather than uniform independent noise

Expected effects:
- noisier allele observations
- possible calling degradation
- possible phasing correctness instability

---

### 4. Read length variation
Read-length realism.

Main parameters:
- `--len-model`
- `--ln-mean`
- `--ln-sigma`
- `--min-len`
- `--max-len`

What it simulates:
- more realistic read length distributions
- changes in how far a read can link variants

Expected effects:
- longer reads can improve continuity
- end-to-end effects may still be mixed if correctness changes as well

---

### 5. Truth indels
Truth-level realism.

Main parameters:
- `--num-indels`
- `--indel-min-len`
- `--indel-max-len`
- `--indel-het-rate`

What it simulates:
- more realistic variant structure
- representation mismatch and alignment/calling stress around indel regions

Recommended handling:
- use `--phase-snps-only`
- use `--eval-snps-only`

This keeps SNP phasing evaluation meaningful while still allowing indels to act as realism stressors.

---

## Supporting realism controls

### Reference complexity preset
Main parameter:
- `--ref-preset {plain,toy,realistic}`

This changes the overall synthetic reference complexity and can be combined with more explicit stressors.

### Oracle vs called regime
Main parameter:
- `--vcf-source {called,oracle,both}`

This is not a biological stressor, but it is central to the benchmarking design because it separates:
- phaser-limited behavior
- realistic end-to-end behavior

---

## Why the realism knobs matter

Without realism knobs, the benchmark would mostly show ideal-case algorithm behavior.

With realism knobs, the benchmark can answer questions such as:
- Is a performance drop caused by calling or by phasing?
- Which stressor hurts most?
- Do stressors interact non-linearly?
- Which parameter settings remain strong under difficult conditions?

This is what makes the project useful for optimization research rather than only accuracy reporting.
