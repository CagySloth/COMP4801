## Pipeline validation (acceptance suite)

You do not need to “prove realism” for every knob before experiments.  
Instead, validate that:

1) each knob affects metrics in the expected direction  
2) knobs do not break the end-to-end pipeline  
3) oracle-vs-called behavior is consistent  
4) runs are reproducible under fixed seeds

---

# What to check in `.pipeline.json`

## Callset sanity
- `callset.precision` should usually be near 1.0 in synthetic settings
- `callset.recall` should decrease when realism gets harder (duplications, bursts, dropout, indels)

## Phasing sanity (oracle vs called)
- oracle phasing should be significantly better than called (isolates phasing from calling errors)
- `num_phase_sets` increases with dropout / reduced connectivity
- `switch_error_rate` increases with ambiguous mapping / wrong constraints

## Indels mode
If `num_indels > 0`, run with `--phase-snps-only`.
Check that the runner writes:
- `<prefix>.truth.snps.vcf.gz`
- `<prefix>.<called|oracle>.snps.vcf.gz`
and oracle phasing stays high on SNPs.

---

# Recommended acceptance suite (fast)

Baseline:

```bash
python -m benchmark.longread_pipeline_runner \
  --prefix output/acc_baseline \
  --seed 0 \
  --ref-length 80000 --num-snps 800 \
  --num-reads 200 --min-len 2000 --max-len 6000 \
  --platform ont --ont-profile q20 \
  --vcf-source both
```

Knobs one-by-one:
- duplications: `--dup-segments 5 ...`
- dropout: `--start-model dropout --dropout-fraction 0.1 ...`
- bursts: `--burst-prob 0.6 --burst-mult 8 ...`
- indels: `--num-indels 80 ... --phase-snps-only`

Integration “hard case”:
- combine multiple knobs; ensure the pipeline completes and oracle remains reasonable.

---

# Reproducibility check

Run the same command twice with the same `--seed` and compare key metrics in `.pipeline.json`.
Small differences are acceptable only if you intentionally include timestamps; otherwise, metrics should match.
