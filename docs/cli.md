# CLI guide

The repository is designed to be run through Python module entrypoints.

## Matrix-track entrypoints

- `python -m dataset.simulate`
- `python -m algorithms.cli.phase`
- `python -m benchmark.benchmark_runner`
- `python -m benchmark.benchmark_accuracy`

## Long-read entrypoints

- `python -m dataset.longread.reference`
- `python -m dataset.longread.truth`
- `python -m dataset.longread.readsim`
- `python -m benchmark.longread_pipeline_runner`
- `python -m benchmark.vcf_phase_eval`
- `python -m benchmark.aggregate_pipeline_reports`
- `python -m benchmark.experiment_driver`

## Inspect available phasing backends

```bash
python -m algorithms.cli.phase --help
```

Current subcommands include:

### Diploid
- `diploid-em`
- `diploid-mst`
- `diploid-whats`
- `diploid-whats-bam`

### Polyploid
- `polyploid-em`
- `polyploid-spectral`

## Recommended style

Prefer commands such as:

```bash
python -m benchmark.longread_pipeline_runner --prefix output/demo ...
```

rather than ad-hoc scripts. This keeps runs reproducible and easier to document.

## Help and self-documentation

All major entrypoints expose `--help`, for example:

```bash
python -m benchmark.longread_pipeline_runner --help
python -m dataset.longread.readsim --help
python -m algorithms.cli.phase diploid-whats-bam --help
```
