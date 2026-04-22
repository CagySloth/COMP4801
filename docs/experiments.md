# Experiments

This document summarizes the main experiment families supported by the repository and the final report.

The long-read study is organized around a few clear questions:

1. What is the baseline attribution picture?
2. Which realism stressors matter most?
3. What happens when stresses combine?
4. Which practical tuning choices improve end-to-end phasing?

---

## Baseline setup used in the report

Common defaults in the final report include:
- `ref_length = 80000`
- `num_snps = 800`
- `het_rate = 0.8`
- `num_reads = 200` unless varied
- ONT-like `q20` profile
- read length range `2000–6000`

These values define the default synthetic scenario for most experiments unless a section explicitly varies one parameter.

---

## Main experiment families

### 1. Baseline attribution: depth scaling
Purpose:
- compare oracle vs called across increasing depth
- identify whether the bottleneck is mostly in calling or phasing

Typical sweep:
- `num_reads ∈ {50, 100, 200, 400}`

Headline metrics:
- `call_recall`
- `oracle_effective_phased_recall`
- `called_effective_phased_recall`
- phase-set counts

---

### 2. Isolated realism stressors
Purpose:
- vary one stressor at a time while holding others fixed

Main stressor families:
- duplicated regions
- coverage dropout
- correlated error bursts
- read length model
- truth indels with SNP-only policy

This is the best place to measure which stressors are mild, strong, noisy, or primarily methodological.

---

### 3. Interaction study
Purpose:
- test whether multiple stressors combine more severely than isolated sweeps suggest

Main combined study in the report:
- duplication × dropout

---

### 4. Optimization under hard conditions
Purpose:
- choose one composite hard scenario
- screen practical tuning knobs
- compare representative configurations

Main practical knobs studied in the report include:
- `call_min_baseq`
- `call_min_mapq`
- `min_baseq`
- `min_mapq`
- `max_coverage`

---

### 5. Robustness and scaling
Purpose:
- check whether tuned settings generalize across scenarios
- examine solve-stage scaling with increasing SNP density

Examples:
- robustness across baseline / dropout / interaction / hard conditions
- solve-time scaling under increasing `num_snps`

---

## Using the experiment driver

Run everything:

```bash
python -m benchmark.experiment_driver \
  --outdir output/experiments_final \
  --seeds 0 1 2
```

Run selected sections only:

```bash
python -m benchmark.experiment_driver \
  --outdir output/experiments_final \
  --seeds 0 1 2 \
  --only depth,dropout,dup,interaction,optimize
```

Other available sections include:
- `bursts`
- `indels`
- `lenmodel`
- `reality`

---

## Reporting strategy

For most plots, report both:
- **oracle** metrics
- **called** metrics

This makes attribution much easier.

Recommended headline metric:
- `called_effective_phased_recall`

Recommended supporting metrics:
- `call_recall`
- `called_num_phase_sets`
- `called_switch_error`
- `called_shared_het_recall`

---

## Practical result summary from the final report

The report’s high-level findings are:
- oracle vs called comparison shows that upstream variant recovery is often the main bottleneck
- coverage dropout is the strongest isolated realism stressor
- duplication + dropout is the clearest compounded weakness
- moderate caller-side and phasing-side tuning gives the best practical improvement

This document should therefore be read as both an experiment-planning guide and a map of the final report’s experiment structure.
