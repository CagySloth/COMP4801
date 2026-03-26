## CLIs — Running the repo reproducibly

Recommended style: run everything as Python modules.

Matrix track:
- `python -m dataset.simulate`
- `python -m algorithms.cli.phase`
- `python -m benchmark.benchmark_runner`
- `python -m benchmark.benchmark_accuracy`

Long-read track:
- `python -m dataset.longread.reference`
- `python -m dataset.longread.truth`
- `python -m dataset.longread.readsim`
- `python -m benchmark.longread_pipeline_runner`
- `python -m benchmark.vcf_phase_eval`

---

## `python -m algorithms.cli.phase`

Run one phasing algorithm.

Diploid subcommands:
- `diploid-em`
- `diploid-mst`
- `diploid-whats`
- `diploid-whats-bam`

Example (long-read phasing on BAM+VCF):
```bash
python -m algorithms.cli.phase diploid-whats-bam \
  --bam output/demo.bam \
  --vcf output/demo.called.vcf.gz \
  --output-prefix output/demo.ws \
  --output-vcf output/demo.ws.phased.vcf
```

---

## `python -m benchmark.longread_pipeline_runner`

End-to-end automation (Steps 1–7).

Key options:
- `--vcf-source {called,oracle,both}`
- `--ref-preset {plain,toy,realistic}`
- `--dup-segments/--dup-len/--dup-min-gap`
- `--len-model {uniform,lognormal}`
- `--start-model {uniform,dropout}`
- `--burst-*`
- `--num-indels ...`
- `--phase-snps-only` (recommended when indels are enabled)

Run `--help` to see the full list.
