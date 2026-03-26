## 📂 `scripts` — Convenience wrappers

This directory contains shell snippets intended as convenience wrappers.  
Prefer the module commands shown in the README and docs when possible.

---

## Matrix track snippet

```bash
PREFIX=output/demo

python -m dataset.simulate -p 2 -n 200 -r 200 -l 60 -e 0.01 -m 0.0 --seed 0 -o "$PREFIX"

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

## Long-read pipeline snippet

```bash
python -m benchmark.longread_pipeline_runner \
  --prefix output/lr_demo \
  --seed 0 \
  --ref-length 80000 --num-snps 800 \
  --num-reads 200 --min-len 2000 --max-len 6000 \
  --platform ont --ont-profile q20 \
  --vcf-source both
```

Indels mode:
```bash
python -m benchmark.longread_pipeline_runner ... --num-indels 80 ... --phase-snps-only
```
