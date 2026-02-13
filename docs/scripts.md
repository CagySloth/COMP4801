## 📂 `scripts` — Convenience shell wrappers

This directory contains shell scripts intended as *convenience wrappers*.

Some scripts may be older examples and might require editing to match the latest CLI/module names.
When in doubt, prefer the module commands shown in the README and docs.

---

## Matrix track: recommended snippet

```bash
PREFIX=output/demo

python -m dataset.simulate -p 2 -n 200 -r 200 -l 60 -o "$PREFIX"

python -m algorithms.cli.phase diploid-whats \
  -i "$PREFIX.reads.npz" \
  --vcf "$PREFIX.vcf" \
  --output-prefix "${PREFIX}_phased" \
  --solver whatshap

python -m benchmark.benchmark_accuracy \
  --truth "$PREFIX.haplotypes.tsv" \
  --pred  "${PREFIX}_phased.haplotypes.tsv" \
  --output "${PREFIX}_phased.accuracy.json"
```

---

## Long-read track: recommended automation

The most up-to-date way to run the long-read end-to-end pipeline is:

```bash
python -m benchmark.longread_pipeline_runner \
  --prefix output/lr_demo \
  --seed 0
```

This runs Steps 1–7 and writes a pipeline report JSON alongside all intermediate artifacts.

If you prefer the manual steps, see `docs/longread_pipeline.md`.
