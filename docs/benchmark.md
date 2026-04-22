# Benchmarking, evaluation, and experiment orchestration

The `benchmark/` package contains the runnable benchmarking logic for both tracks, with the long-read end-to-end pipeline as the main research workflow.

---

## Main components

### `benchmark.longread_pipeline_runner`
Canonical single-run long-read workflow.

Responsibilities:
- generate reference and truth
- simulate reads
- align reads with `minimap2` + `samtools`
- call variants with `bcftools`
- phase oracle and/or called VCFs
- evaluate phased output
- write `<prefix>.pipeline.json`

This is the main entrypoint for reproducible long-read runs.

### `benchmark.vcf_phase_eval`
Evaluates phased VCF output against phased truth.

Key properties:
- block-flip aware
- matches records by `(CHROM, POS, REF, ALT)`
- focused on diploid SNP phasing
- used by both oracle and called workflows

### `benchmark.aggregate_pipeline_reports`
Collects many `*.pipeline.json` files under a root directory and writes a single `aggregate.csv`.

This CSV is the standard interface for:
- experiment summaries
- plots
- report tables

### `benchmark.experiment_driver`
Runs complete experiment suites by repeatedly calling the single-run pipeline and then aggregating the results.

Supported experiment families include:
- smoke
- depth
- dropout
- duplication
- bursts
- indels
- read-length model
- duplication × dropout interaction
- optimization
- reality check

### Legacy / supporting utilities
- `benchmark.benchmark_runner`
- `benchmark.benchmark_accuracy`
- `benchmark.longread_baseline_grid`
- plotting helpers such as `plot_baseline_results.py`

---

## Core long-read report structure

Each long-read run writes `<prefix>.pipeline.json`.

Typical contents include:
- input parameters
- output artifact paths
- callset comparison metrics
- phasing results for `oracle` and/or `called`
- derived metrics used in experiment plots

This file is the **source of truth** for downstream aggregation.

---

## Headline metrics in long-read experiments

The most important metrics for the long-read research are:

- `call_recall`
- `call_precision`
- `oracle_effective_phased_recall`
- `called_effective_phased_recall`
- `oracle_num_phase_sets`
- `called_num_phase_sets`
- `oracle_switch_error`
- `called_switch_error`

The most important attribution idea is:

- **oracle** = phaser-limited behavior
- **called** = end-to-end behavior

The gap between oracle and called effective phased recall is therefore a key signal for where performance is being lost.

---

## Typical workflow

### 1. Run many prefixes

```bash
python -m benchmark.experiment_driver \
  --outdir output/experiments \
  --seeds 0 1 2
```

### 2. Aggregate run-level JSON files

```bash
python -m benchmark.aggregate_pipeline_reports \
  --root output/experiments/01_depth \
  --out  output/experiments/01_depth/aggregate.csv
```

### 3. Plot or inspect aggregate results

Use the generated `aggregate.csv` files directly in report or poster plotting scripts.

---

## Why the design matters

The benchmark package is built around a simple rule:

- **one canonical runner for a single run**
- **one canonical aggregator for many runs**

This makes the project:
- easier to reproduce
- easier to validate
- easier to trace from final figure back to raw outputs
