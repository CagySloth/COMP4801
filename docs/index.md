# Documentation index

This repository combines a legacy matrix-track phasing suite with a newer long-read end-to-end benchmarking pipeline.

## Recommended reading order

### Core project docs
1. [Long-read pipeline](longread_pipeline.md) — canonical end-to-end workflow
2. [Realism knobs](realism_knobs.md) — stressors and what they simulate
3. [Benchmarking](benchmark.md) — runners, evaluator, aggregation, experiment driver
4. [Experiments](experiments.md) — experiment families and reporting strategy
5. [Validation](validation.md) — acceptance checks and testing guidance
6. [I/O formats](io.md) — artifact contracts and output files

### Supporting docs
- [Dataset generation](dataset.md)
- [Algorithms](algorithms.md)
- [CLIs](cli.md)
- [Examples](examples.md)

## Project structure at a glance

- `dataset/` — synthetic data generation
- `algorithms/` — phasing backends and adapters
- `benchmark/` — orchestration, evaluation, aggregation, experiments
- `tests/` — unit and integration tests
- `report/` — final report source and build assets

## Which track should I use?

- Use the **matrix track** if you want a fast, controlled algorithm test with `.reads.npz` inputs.
- Use the **long-read track** if you want the main project workflow: FASTA → FASTQ → BAM → VCF → phased VCF → evaluation.
