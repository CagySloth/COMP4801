## Experiment planning (recommended structure)

Once the acceptance suite passes, **freeze the pipeline** and move to experiments.  
A good final-report plan is:

---

# A) Baseline curves (ONT q20)

Purpose: show how performance scales with read count.

Sweep:
- `num_reads`: 50, 100, 200, 400
Fixed:
- `ref_length`, `num_snps`, read length range, `max_coverage`, `min_mapq/baseq`

Report:
- call recall
- called effective phased recall
- called num phase sets
- oracle effective phased recall (control)

Tools:
- `benchmark.longread_baseline_grid.py`
- `benchmark.plot_baseline_results.py`

---

# B) One knob at a time (impact analysis)

Recommended knobs:
1) duplications: `dup_segments` {0, 5}
2) dropout: `dropout_fraction` {0, 0.05, 0.1, 0.2}
3) bursts: `burst_prob` {0, 0.3, 0.6}
4) indels: `num_indels` {0, 80, 200} (use `--phase-snps-only`)

For each sweep:
- 3–10 seeds depending on time
- aggregate with `benchmark.aggregate_pipeline_reports`
- plot means ± std

---

# C) “Hard realistic” setting + WhatsHap parameter tuning

Choose one combined scenario (e.g., duplications + dropout + bursts + indels).  
Then sweep WhatsHap parameters:

- `max_coverage`: 10, 15, 25
- `min_mapq`: 0, 10, 20
- `min_baseq`: 0, 10, 20

Goal: show tradeoffs:
- higher coverage improves recall but can increase runtime
- stricter filtering may reduce errors but lose connectivity

---

# D) How to report results

In most plots, report BOTH:
- called (end-to-end, realistic)
- oracle (upper bound for phasing given perfect variant set)

Suggested “headline metric”:
- `called_effective_phased_recall` (captures calling+phasing)
And supporting metrics:
- call recall
- switch error rate
- num phase sets (fragmentation)

---

# Notes on scope

Current evaluation focuses on **biallelic SNP phasing**.  
Indels are included to stress alignment/calling realism; SNP-only phasing/eval mode keeps results interpretable.
