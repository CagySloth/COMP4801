## 📂 `scripts` — Convenience shell wrappers

This directory contains shell scripts intended as *convenience wrappers*.

Some scripts may be older examples and might require editing to match the latest CLI/module names.
When in doubt, prefer the module commands shown in the README:

- `python -m dataset.simulate`
- `python -m algorithms.cli.phase`
- `python -m benchmark.benchmark_runner`

---

## `simulate.sh`

A thin wrapper over the simulator:

```bash
./scripts/simulate.sh -p 2 -n 200 -r 200 -l 60 -o output/demo
```

---

## Recommended “pipeline” snippet

If you want an end-to-end run in bash, this is the most up-to-date shape:

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
