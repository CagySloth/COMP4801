## 📂 `benchmark` — Benchmarking and evaluation

This folder contains:

1) **Matrix-track benchmarking**
   - `benchmark.benchmark_runner`: parameter sweeps (simulate → phase → score)
   - `benchmark.benchmark_accuracy`: TSV haplotype accuracy

2) **VCF-track evaluation**
   - `benchmark.vcf_phase_eval`: phased VCF vs phased truth VCF (block-flip aware)
     - Matches variants by **(CHROM,POS,REF,ALT)**
     - Focuses on **biallelic SNPs** by default

3) **Long-read end-to-end automation**
   - `benchmark.longread_pipeline_runner`: runs Steps 1–7 and writes `<prefix>.pipeline.json`

4) **Aggregation + plotting**
   - `benchmark.aggregate_pipeline_reports`: gather many pipeline JSONs into a single CSV
   - `benchmark.longread_baseline_grid`: run baseline read-depth curves
   - `benchmark.plot_baseline_results`, `benchmark.plot_refpreset_compare`: plotting helpers

---

# Long-read pipeline report (`.pipeline.json`)

The pipeline report includes:

- `callset`: truth vs called SNP set comparison (precision/recall)
- `phasing_runs.called` and/or `phasing_runs.oracle`
  - `eval`: phase accuracy, switch error, num phase sets
  - `derived`: effective phased recall, phasing rate, shared het recall, etc.

Recommended headline metric for end-to-end performance:
- **called_effective_phased_recall**

---

# Aggregating experiments

After running multiple prefixes under a folder:

```bash
python -m benchmark.aggregate_pipeline_reports \
  --root output/exp_folder \
  --out  output/exp_folder/aggregate.csv
```

Then plot mean ± std vs your sweep variable.

---

# Baseline grid (example)

```bash
python -m benchmark.longread_baseline_grid \
  --outdir output/exp_baseline_q20 \
  --ont-profile q20
python -m benchmark.plot_baseline_results \
  --csv output/exp_baseline_q20/aggregate.csv \
  --outdir output/exp_baseline_q20/plots
```

---

# Realism comparison script

`benchmark/realism_comparison.sh` runs a preset comparison and writes an aggregate CSV.  
Use it as an example of “batch + aggregate + plot” workflow.
