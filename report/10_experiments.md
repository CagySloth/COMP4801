## 10. Experimental Evaluation

### 10.1 Experimental methodology and common setup

#### 10.1.1 Environment and toolchain

- OS / CPU: macOS / Apple M2 Max
- Python version: 3.12.12
- minimap2 version: 2.30-r1287
- samtools version: 1.23
- bcftools version: 1.23
- Vendored WhatsHap core path: `vendor/whatshap_core/whatshap/core.<extension>`

#### 10.1.2 Pipeline regimes and evaluation policy

- **Oracle regime:** Phase using `*.oracle.vcf.gz` to measure phaser-limited performance under correct variant sites.
- **Called regime:** Phase using `*.called.vcf.gz` to measure end-to-end performance including calling limitations.
- **Indel policy:** When truth indels are enabled, use SNP-only phasing/evaluation (`phase_snps_only`, `eval_snps_only`) to avoid representation-induced distortion in SNP phasing metrics.

A custom adaptor layer was used to construct the `ReadSet` consumed by the vendored WhatsHap core. This adaptor was necessary to support oracle-vs-called benchmarking, unified JSON provenance, and controlled SNP-only phasing/evaluation under indel-containing truth sets. It is therefore part of the project’s benchmarking infrastructure rather than a direct reproduction of the standard production-facing WhatsHap front-end.

#### 10.1.3 Baseline defaults

Unless stated otherwise, experiments use the baseline defaults defined by the experiment driver.

- Reference / truth:
  - `ref_length = 80 kb`
  - `num_snps = 800`
  - `het_rate = 0.8`
  - No duplications, no dropout, no bursts, no truth indels
- Read simulation:
  - ONT-like simulation profile (`ont_profile = q20`)
  - `num_reads = 200` unless varied
  - Read length range `2–6 kb`
  - `start_model = uniform` unless varied
- Calling / phasing:
  - `vcf_source = both`
  - Baseline calling and phasing thresholds as defined by the experiment driver
- Seeds:
  - All reported results are aggregated across multiple seeds into `aggregate.csv`

#### 10.1.4 Metrics and interpretation

Metrics are computed per run from `*.pipeline.json` and `*.eval.json`, then aggregated across seeds into `aggregate.csv`.

- **Calling quality**
  - `call_precision`
  - `call_recall`
  - `shared_snps`
- **Phasing completeness / usable output**
  - `called_shared_het_recall`
  - `called_phasing_rate_shared_het`
  - `oracle_effective_phased_recall`
  - `called_effective_phased_recall`
- **Phasing correctness**
  - `oracle_switch_error`
  - `called_switch_error`
  - `oracle_phase_accuracy`
  - `called_phase_accuracy`
- **Fragmentation**
  - `oracle_num_phase_sets`
  - `called_num_phase_sets`
- **Runtime**
  - `time_total_sec`
  - `called_time_*`

Oracle metrics indicate phaser-limited performance when correct sites are supplied. Called metrics indicate end-to-end performance after alignment, variant calling, and phasing. The gap between oracle and called effective phased recall is therefore a key attribution signal: it measures how much performance is lost because the correct heterozygous sites are not recovered and preserved into the phasing stage.

---

### 10.2 Baseline attribution study: depth scaling

#### Aim
Establish baseline performance curves as sequencing depth increases, and attribute performance limitations by comparing the oracle (phaser-limited) regime against the called (end-to-end) regime.

#### Setup

- Varied: `num_reads ∈ {50, 100, 200, 400}` (3 seeds each)
- Fixed: ONT-like simulation profile (`q20`), reference length `80 kb`, `800` truth SNPs (`het_rate = 0.8`), read length range `2–6 kb`, and the baseline calling/phasing thresholds

**Figure 10.2.3 — Called effective phased recall vs sequencing depth (end-to-end performance).**
Called effective phased recall measures the fraction of truth heterozygous SNPs that end up both present in the called set and correctly phased (after best block flip). The curve is substantially lower than the oracle upper bound at low depth, showing that end-to-end performance is dominated early by limited variant recovery and callset overlap.

#### Results summary (means over seeds)

Calling (called framework):

- `call_recall`: 0.224 → 0.399 → 0.662 → 0.877
- `call_precision`: 1.000 at all depths in this synthetic setup

Oracle phasing (phaser-limited):

- `oracle_effective_phased_recall`: 0.513 → 0.750 → 0.921 → 0.969
- `oracle_switch_error`: near 0 across all depths
- `oracle_num_phase_sets`: 10.7 → 5.7 → 1.7 → 1.0

Called phasing (end-to-end):

- `called_effective_phased_recall`: 0.040 → 0.140 → 0.426 → 0.845
- `called_shared_het_recall`: 0.129 → 0.281 → 0.590 → 0.854
- `called_phasing_rate_shared_het`: 0.355 → 0.517 → 0.808 → 0.991
- `called_switch_error`: high variance at 100–200 reads, near 0 at 400 reads
- `called_num_phase_sets`: 5.0 → 7.7 → 3.0 → 1.0

These trends are reflected by Figures 10.2.1–10.2.6, especially Figure 10.2.1 for the depth-limited rise in calling recall, Figure 10.2.3 for the large oracle-vs-called gap at low depth, and Figures 10.2.4–10.2.5 for the reduction in fragmentation with increasing coverage

#### Key observations

- **O1** Calling is depth-limited: Recall of calling increased significantly from 0.224 to 0.877 with `num_reads`. On the other hand, precision is maintained at 1.0 in this setup, indicating that errors mainly appear as missed sites instead of false positives.

- **O2** WhatsHap is accurate given correct variant sites: In the oracle framework, a near 0 switch error rate and a near 1 phase accuracy are maintained, even at low depth. The major limitation at `num_reads` 50-100 is the connectivity gaps, leaving some heterozygous sites unphased.

- **O3** Residual high-depth loss is attributed by variant caller: At `num_reads` 400, in the called framework, phasing is nearly perfect on shared heterozygous sites, achieving `called_phasing_rate_shared_het` ≈ 0.99 and accuracy ≈ 1.0. The remaining effective phased recall gap between the oracle and the called framework is mostly due to missing heterozygous sites in the called VCFs, with `called_shared_het_recall` ≈ 0.85.

- **O4** Mid-depth instability in the called regime: At `num_reads` 100-200, occasional high switch error can be observed in the called framework in some seeds. This could likely be due to seed sensitivity and limited number of informative heterozygous sites.

#### Takeaway
The baseline depth sweep helps confirm that the pipeline is behaving sensibly and as expected, where recall in called regime increases with depth, and phasing becomes both accurate and contiguous in the oracle regime. While end-to-end (called regime) phasing performance improves sharply with depth, but its accuracy remains primarily limited by the overlap of correctly called heterozygous variants. On available sites, WhatsHap's phasing is almost perfect at a high level of depth, therefore, instead of algorithmic improvement on WhatsHap, improvement from variant calling overlap would improve phasing performance more significantly.

---

### 10.3 Isolated realism stress studies

#### Aim
Evaluate the impact of individual realism knobs on calling recall, phasing completeness, correctness, and fragmentation, while the other baseline settings are kept constant.

#### 10.3.1 Duplicated regions

##### Setup

- Varied: `dup_segments ∈ {0, 1, 3, 5}`
- Fixed:
  - `dup_len = 3000`
  - `dup_min_gap = 500`
  - Baseline depth `num_reads = 200`
  - No dropout, no bursts, no truth indels
  - Phasing input source `vcf_source = both`

##### Results summary

- `call_recall`: 0.663 → 0.671 → 0.661 → 0.666
- `called_effective_phased_recall`: 0.436 → 0.473 → 0.479 → 0.372
- `called_shared_het_recall`: 0.590 → 0.599 → 0.586 → 0.593
- `called_phasing_rate_shared_het`: 0.819 → 0.836 → 0.848 → 0.778
- `called_phase_accuracy`: 0.890 → 0.920 → 0.960 → 0.802
- `called_switch_error`: 0.122 → 0.073 → 0.040 → 0.209
- `called_num_phase_sets`: 3.4 → 3.4 → 3.4 → 2.6
- `oracle_effective_phased_recall`: 0.925 → 0.932 → 0.924 → 0.918

These trends are reflected by Figures 10.3.1.1–10.3.1.5, especially Figure 10.3.1.2 for the drop in called effective phased recall only at the strongest duplication setting and Figure 10.3.1.3 for the corresponding switch-error increase.

##### Key observations

- **O1** Recall rate in variant-caller is mostly unaffected by duplication: Calling recall remains nearly unchanged across `dup_segments ∈ {0,1,3,5}`, and call precision remains at 1.0.

- **O2** Oracle phasing remains stable under duplication: Effective phased recall remains high and nearly unchanged in the oracle regime, and the switch error rate stays near 0 across the sweep.

- **O3** Degradation in the called regime mostly appears at the highest duplication level: Effective phased recall in the called regime drops from 0.479 to 0.372, and switch error rate increases from 0.040 to 0.209 at `dup_segments` = 5. The main cause of the degradation is not due to loss of shared heterozygous sites, as `called_shared_het_recall` remains nearly unchanged.

- **O4** Main loss is in phasing completeness and correctness on surviving sites: The strongest duplication setting reduces `called_phasing_rate_shared_het` and `called_phase_accuracy`, indicating that duplicated regions mainly weaken the consistency of phasing evidence in the called regime.

##### Takeaway\\
Under the current parameterization, duplicated regions are a relatively mild stressor for variant recovery and for oracle phasing. Their main impact appears only at the highest duplication level, where end-to-end phasing becomes less complete and less accurate despite largely unchanged callset overlap.

#### 10.3.2 Coverage dropout

##### Setup

- Varied: `dropout_fraction ∈ {0.00, 0.05, 0.10, 0.20}`
- Fixed:
  - `dropout_block_len = 1000`
  - Baseline depth `num_reads = 200`
  - No duplications, no bursts, no truth indels
  - `start_model = dropout`
  - Phasing input source `vcf_source = both`

\begin{center}
\includegraphics[width=0.82\linewidth]{figures/10_experiments/fig_10_3_2_4_called_effective_phased_recall_dropout.png}
\end{center}

**Figure 10.3.2.4 — Called effective phased recall vs dropout fraction.**
Called effective phased recall declines much more steeply than oracle effective phased recall, showing that dropout affects both phasing connectivity and, at stronger dropout levels, variant recovery in the called regime.

##### Results summary

- `call_recall`: 0.663 → 0.635 → 0.605 → 0.432
- `oracle_effective_phased_recall`: 0.925 → 0.865 → 0.857 → 0.637
- `oracle_num_phase_sets`: 1.8 → 2.8 → 3.2 → 10.0
- `called_effective_phased_recall`: 0.436 → 0.359 → 0.269 → 0.081
- `called_shared_het_recall`: 0.590 → 0.565 → 0.529 → 0.362
- `called_phasing_rate_shared_het`: 0.819 → 0.704 → 0.514 → 0.239
- `called_num_phase_sets`: 3.4 → 5.0 → 6.2 → 4.0

These trends are reflected by Figures 10.3.2.1–10.3.2.5, especially Figure 10.3.2.4 for the steep decline in called effective phased recall and Figure 10.3.2.1 for the strong rise in oracle phase-set count.

##### Key observations

- **O1** Coverage dropout strongly increases fragmentation: In the oracle regime, the number of phase sets rises from 1.8 at `dropout_fraction` = 0.00, to 10.0 at 0.20, while oracle switch error remains near zero.

- **O2** Dropout directly reduces phasing completeness even with correct variant sites: As dropout increases, effective phased recall drops from 0.925 to 0.637 in the oracle regime, indicating that weaker read connectivity alone can lower phasing completeness.

- **O3** End-to-end performance degrades more than oracle performance: Effective phased recall falls from 0.436 to 0.081 in the called regime, indicating that dropout hurts end-to-end phasing performance through both lossing connectivity, and reducing overlap of correctly called heterozygous variants.

- **O4** Moderate dropout is mainly a problem of connectivity and completeness: The major degradation in the called regime from `dropout_fraction` = 0.00 to 0.10, is a drop in `called_phasing_rate_shared_het` and increase in fragmentation.

- **O5** Severe dropout becomes a mixed bottleneck in both variant calling and phasing: At `dropout_fraction` = 0.20, `call_recall` drops sharply to 0.432, `called_shared_het_recall` falls to 0.362, and `called_phasing_rate_shared_het` falls to 0.239.

##### Takeaway\\
In this study, coverage dropout is that most impactful isolated stressor. The primary effect of coverage dropout is the creation of low-coverage gaps, preventing reads from connecting heterozygous sites into extended phase blocks, thus reducing phasing completeness, even in the oracle regime. In the called regime, this effect is intensified, both phase connectivity and variant recovery declined significantly at higher levels of dropout.

#### 10.3.3 Correlated error bursts

##### Setup

- Varied: `burst_prob ∈ {0.0, 0.3, 0.6}`
- Fixed:
  - `burst_count = 2`, `burst_len = 300`, `burst_mult = 8.0`
  - Baseline depth `num_reads = 200`
  - No duplications, no coverage dropout, no truth indels
  - `start_model = uniform`
  - Phasing input source `vcf_source = both`

##### Results summary

- `call_recall`: 0.663 → 0.681 → 0.676
- `oracle_effective_phased_recall`: 0.925 → 0.926 → 0.921
- `called_effective_phased_recall`: 0.436 → 0.369 → 0.455
- `called_phase_accuracy`: 0.890 → 0.766 → 0.887
- `called_switch_error`: 0.122 → 0.142 → 0.108

These trends are reflected by Figures 10.3.3.1–10.3.3.5, especially Figures 10.3.3.1 and 10.3.3.2, where the intermediate-burst dip is visible but non-monotonic.

##### Key observations

- **O1** Recall in the called regime remains mostly unaffected, while phasing in the oracle regime remains stable.
- **O2** The main observed effect is a drop in phasing correctness at `burst_prob` = 0.3 in the called regime.
- **O3** The effect is not monotonic as almost all metric returned to near the baseline at `burst_prob` = 0.6.

##### Takeaway
Under the current parameterization, correlated error bursts are a weak and somewhat noisy called-regime correctness stressor rather than a strong practical failure mode.


#### 10.3.4 Read length model

##### Setup

- Varied: `len_model ∈ {uniform, lognormal}`
- Fixed:
  - Baseline depth `num_reads = 200`
  - Lognormal parameters `ln_mean = 8.3`, `ln_sigma = 0.6`
  - No duplications, no coverage dropout, no correlated bursts, no truth indels
  - `start_model = uniform`
  - Phasing input source `vcf_source = both`

##### Results summary

- `call_recall`: 0.663 → 0.701
- `oracle_effective_phased_recall`: 0.925 → 0.931
- `oracle_num_phase_sets`: 1.8 → 1.4
- `called_effective_phased_recall`: 0.436 → 0.443
- `called_shared_het_recall`: 0.590 → 0.633
- `called_phase_accuracy`: 0.890 → 0.855
- `called_switch_error`: 0.122 → 0.165
- `called_num_phase_sets`: 3.4 → 1.6

These trends are reflected by Figures 10.3.4.1–10.3.4.5, especially Figure 10.3.4.2 for the continuity gain and Figure 10.3.4.4 for the small net change in called effective phased recall.

##### Key observations

- **O1** The lognormal model improves callset overlap and phase-block continuity.
- **O2** End-to-end effective phased recall changes only slightly.
- **O3** The gains in continuity and overlap are largely offset by slightly worse called-regime phasing correctness.

##### Takeaway\\
The read length model is informative for mechanism analysis, but under the present parameterization, it is a mixed continuity and correctness trade-off rather than a clear optimization lever.


#### 10.3.5 Truth indels and SNP-only policy

##### Setup

- Varied: `num_indels ∈ {0, 80, 200}`
- Fixed:
  - `phase_snps_only = true`
  - `eval_snps_only = true`
  - Baseline depth `num_reads = 200`
  - No duplications, no dropout, no bursts
  - Phasing input source `vcf_source = both`

##### Results summary

- `call_recall`: 0.663 → 0.691 → 0.654
- `call_precision`: 1.000 → 0.999 → 0.997
- `oracle_effective_phased_recall`: 0.925 → 0.925 → 0.921
- `called_effective_phased_recall`: 0.436 → 0.573 → 0.493
- `called_shared_het_recall`: 0.590 → 0.623 → 0.577
- `called_num_phase_sets`: 3.4 → 2.4 → 4.0

These trends are reflected by Figures 10.3.5.1–10.3.5.4, which show that the main result is evaluation robustness rather than a monotonic degradation trend.

##### Key observations

- **O1** Oracle SNP-phasing remains stable when truth indels are introduced.
- **O2** Called SNP-phasing does not collapse in the presence of indels.
- **O3** At `num_indels` = 80, phasing performance is slightly better than baseline on several called metrics, but this should be interpreted as ordinary stochastic variation and indirect alignment and calling effects rather than evidence that indels improve SNP phasing.
- **O4** Even stronger indel conditions remain manageable under the SNP-only policy.

##### Takeaway\\
The SNP-only policy successfully preserves interpretable SNP phasing evaluation in the presence of truth indels. This is primarily a methodology-validity result rather than a performance stress study.

---

### 10.4 Interaction study: duplicated regions × coverage dropout

#### Aim
Test whether duplicated regions and coverage dropout compound non-linearly, and determine whether the resulting weakness is mainly caller-limited, connectivity-limited, or mixed.

#### Setup

- `dup_segments ∈ {0, 5}`
- `dropout_fraction ∈ {0.0, 0.1}`
- Fixed:
  - `dropout_block_len = 1000`
  - `dup_len = 3000`
  - `dup_min_gap = 500`
  - Baseline depth `num_reads = 200`
  - No bursts, no truth indels
  - Phasing input source `vcf_source = both`

#### Results summary

- `call_recall`:
  - `dup=0, dropout=0.0`: 0.663
  - `dup=0, dropout=0.1`: 0.605
  - `dup=5, dropout=0.0`: 0.666
  - `dup=5, dropout=0.1`: 0.582
- `oracle_effective_phased_recall`:
  - 0.925 → 0.857 → 0.918 → 0.847
- `oracle_num_phase_sets`:
  - 1.8 → 3.2 → 1.8 → 4.0
- `called_effective_phased_recall`:
  - 0.436 → 0.269 → 0.373 → 0.197
- `called_shared_het_recall`:
  - 0.590 → 0.529 → 0.593 → 0.501
- `called_phasing_rate_shared_het`:
  - 0.819 → 0.514 → 0.778 → 0.410
- `called_num_phase_sets`:
  - 3.4 → 6.2 → 2.6 → 6.2

These trends are reflected by Figures 10.4.1–10.4.4, especially Figure 10.4.2 for the compounded drop in called effective phased recall and Figure 10.4.3 for the dropout-driven fragmentation pattern.

\begin{center}
\includegraphics[width=0.82\linewidth]{figures/10_experiments/fig_10_4_2_called_effective_phased_recall_interaction.png}
\end{center}

**Figure 10.4.2 — Called effective phased recall across duplication × dropout conditions.**
The combined duplication + dropout condition produces the lowest end-to-end phased recall, showing that the two stressors interact to create a substantially more difficult phasing problem than either one alone.

#### Key observations

- **O1** The combined condition is the worst end-to-end setting: Called effective phased recall decreases from 0.436 in the baseline to 0.269 with dropout alone, 0.373 with duplication alone, and 0.197 when both stressors are combined.
- **O2** Dropout remains the main driver of phasing fragmentation: In the oracle regime, duplication alone has little effect, whereas dropout increases fragmentation and the combined condition increases it further.
- **O3** The interaction is strongest in the called regime: Relative to dropout alone, the combined condition further reduces `call_recall`, `called_shared_het_recall`, and `called_phasing_rate_shared_het`.
- **O4** The compounded weakness is mixed rather than purely phaser-limited: Oracle effective phased recall decreases only modestly when duplication is added on top of dropout, but called effective phased recall decreases more substantially.

#### Optimization implication
The duplication × dropout interaction suggests that the most important optimization opportunities under compounded stress are not limited to WhatsHap read selection. Preserving phasing evidence remains important, but upstream variant recovery in ambiguous and low-coverage regions is also a major limiting factor.

#### Takeaway
Duplicated regions and coverage dropout interact to produce a more severe end-to-end phasing weakness than either stressor alone. Dropout drives fragmentation, while duplication further reduces callable overlap and phasing completeness in the called regime.

---

### 10.5 Optimization under hard conditions

#### 10.5.1 Hard-scenario setup

Optimization sweeps were conducted on a fixed composite stress condition:

- Duplicated regions: `dup_segments = 5`
- Coverage dropout: `dropout_fraction = 0.1`, `dropout_block_len = 1000`
- Correlated error bursts: `burst_prob = 0.6`, `burst_count = 2`, `burst_len = 300`, `burst_mult = 8.0`
- Truth indels: `num_indels = 120` with SNP-only phasing/evaluation enabled
- `ref_length = 120 kb`
- `num_snps = 1200`
- `num_reads = 300`

The purpose of this setting was to test whether practical parameter tuning can recover performance under combined realistic stresses rather than isolated perturbations.

#### 10.5.2 Screening of individual knobs

##### Max coverage

- Increasing `max_coverage` from 10 to 40 does not improve called effective phased recall, switch error, or phase-set count, but increases runtime substantially.
- A lower-coverage sweep shows that performance plateaus by `max_coverage` = 8–10; `max_coverage` = 4 is slightly too aggressive, while 15 is unnecessarily expensive.

These trends are reflected by Figures 10.5.2.1–10.5.2.4 for the coarse sweep and Figures 10.5.2.13–10.5.2.16 for the lower-coverage refinement, which together show a clear performance plateau with increasing runtime at higher retained coverage.

##### Phasing thresholds

- The phasing quality-threshold grid is separated into two regions:
  - At `min_baseq` = 0-10: In called regime, effective phased recall is better, switch error is lower, less phase sets are constructed.
  - At `min_baseq` = 20: In both called and oracle frameworks, phasing is less accurate, fragmentation is increased.
- Across the tested range, `min_mapq` has negligible impact.
- Through a more detailed sweep, it is observed that phasing performance is effectively flat for `min_baseq` = 0–15 and degrades only at 20.

These trends are reflected by Figures 10.5.2.5–10.5.2.8 for the coarse phasing-threshold grid, and Figures 10.5.2.17–10.5.2.20 for the more detailed `min_baseq` sweep. Together, these figures indicate that, across moderate settings, phasing performance is stable. Degradation occurs only when filtering becomes too strict at high `min_baseq`.

##### Variant-caller thresholds

- The strongest caller-side optimization lever is `call_min_baseq`.
- `call_recall`, `called_shared_het_recall`, and `called_effective_phased_recall` improves significantly when `call_min_baseq` is lowered from 20 to 10 or 5, while call precision stayed around 0.9995.
- `call_min_mapq` has only limited impacts, and does not offer much potential for optimization.

These trends are reflected by Figures 10.5.2.9–10.5.2.12 for the caller-threshold grid, and Figures 10.5.2.25–10.5.2.28 for the `call_min_mapq` rule-out sweep. These figures indicate that `call_min_baseq` is the stronger optimization lever with more optimization potential.

##### Local performance plateau confirmation

- Through a more focused local search around `call_min_baseq ∈ {5, 10, 15}` and `min_baseq ∈ {5, 10, ` `15}`, it is observed that the recommended settings lie on a stable local optimum.
- Caller base quality remains to be the more decision parameter. Without overly strict filtering, phasing base-quality is effectively irrelevant.

These trends are reflected by Figures 10.5.2.21–10.5.2.24, indicating that the recommended threshold pair is lying on a stable local performance plateau.

##### Takeaway\\
Optimization headroom exists, but it is uneven across knobs. The strongest useful signal is moderate caller-side base-quality filtering (`call_min_baseq` = 10), followed by avoiding overly strict phasing base-quality filtering (`min_baseq` = 10). By contrast, `min_mapq`, `call_min_mapq`, and high `max_coverage` values above the plateau are weak or negligible knobs.

#### 10.5.3 Representative configuration comparison

Representative hard-scenario configurations:

- **`default`**
  - `call_min_mapq = 20`
  - `call_min_baseq = 15`
  - `max_coverage = 15`
  - `min_mapq = 20`
  - `min_baseq = 20`
- **`caller_only`**
  - `call_min_mapq = 20`
  - `call_min_baseq = 10`
  - `max_coverage = 15`
  - `min_mapq = 20`
  - `min_baseq = 20`
- **`phasing_only`**
  - `call_min_mapq = 20`
  - `call_min_baseq = 15`
  - `max_coverage = 8`
  - `min_mapq = 20`
  - `min_baseq = 10`
- **`optimized`**
  - `call_min_mapq = 20`
  - `call_min_baseq = 10`
  - `max_coverage = 8`
  - `min_mapq = 20`
  - `min_baseq = 10`
- **`runtime`**
  - `call_min_mapq = 20`
  - `call_min_baseq = 10`
  - `max_coverage = 6`
  - `min_mapq = 20`
  - `min_baseq = 10`

#### Results summary

- `default`: `called_effective_phased_recall = 0.2919`, `called_shared_het_recall = 0.5098`, `called_num_phase_sets = 9.6`, `time_total_sec = 3.9338`
- `caller_only`: `0.2994`, `0.5181`, `9.0`, `4.0540`
- `phasing_only`: `0.2964`, `0.5098`, `6.6`, `3.7970`
- `optimized`: `0.3041`, `0.5181`, `6.6`, `3.8682`
- `runtime`: `0.3039`, `0.5181`, `6.6`, `3.8504`

These trends are reflected by Figures 10.5.3.1–10.5.3.4, and are reinforced by the final default-vs-optimized confirmation in Figures 10.5.3.5–10.5.3.8.

\begin{center}
\includegraphics[width=0.82\linewidth]{figures/10_experiments/fig_10_5_3_1_called_effective_phased_recall_configs.png}
\end{center}

**Figure 10.5.3.1 — Called effective phased recall across representative configurations.**
The optimized configuration achieves the highest end-to-end phased recall, while the runtime-biased configuration performs almost identically. This shows that combining caller-side and phasing-side tuning yields the strongest practical performance.

#### Key observations

- **O1** The optimized configuration provides the best overall trade-off: `optimized` achieves the highest called effective phased recall while also maintaining low fragmentation, very low switch error, and lower runtime than the default configuration.
- **O2** Default is dominated by the tuned configurations: The default configuration is worse than `optimized` and `runtime` on called effective phased recall, called shared heterozygous recall, switch error, and fragmentation.
- **O3** Caller-only and phasing-only tuning recover complementary parts of the gain: `caller_only` improves overlap-related metrics but leaves fragmentation high, whereas `phasing_only` strongly improves continuity and oracle phasing metrics but does not recover additional called-site overlap.
- **O4** The runtime-biased configuration is a viable alternative: `runtime` performs almost identically to `optimized` on the main called metrics while using slightly lower retained coverage.

#### Takeaway
By combining moderate caller-side and phaser-side parameter tuning, the best practical performance can be achieved. The optimized configuration is the preferred setup in general, however, the runtime-oriented configuration can serve as a viable alternative when prioritizing computational efficiency.

#### 10.5.4 Robustness across scenarios

Configurations compared:

- `default`
- `optimized`
- `runtime`

Scenarios compared:

- baseline
- dropout
- interaction
- hard

#### Results summary

- **Baseline**
  - `default`: `called_effective_phased_recall = 0.4364`, `called_shared_het_recall = 0.5902`, `called_num_phase_sets = 3.4`
  - `optimized`: `0.5448`, `0.6097`, `1.4`
  - `runtime`: `0.5448`, `0.6097`, `1.4`
- **Dropout**
  - `default`: `0.2694`, `0.5292`, `6.2`
  - `optimized`: `0.2908`, `0.5420`, `4.2`
  - `runtime`: `0.2908`, `0.5420`, `4.2`
- **Interaction**
  - `default`: `0.1972`, `0.5014`, `6.2`
  - `optimized`: `0.2118`, `0.5094`, `5.4`
  - `runtime`: `0.2118`, `0.5094`, `5.4`
- **Hard**
  - `default`: `0.2919`, `0.5098`, `9.6`
  - `optimized`: `0.3041`, `0.5181`, `6.6`
  - `runtime`: `0.3039`, `0.5181`, `6.6`

These trends are reflected by Figures 10.5.4.1–10.5.4.4, where the tuned configurations remain above the default across all four scenarios.

\begin{center}
\includegraphics[width=0.88\linewidth]{figures/10_experiments/fig_10_5_4_1_called_effective_phased_recall_robustness.png}
\end{center}

**Figure 10.5.4.1 — Called effective phased recall across scenarios and representative configurations.**
Across all tested scenarios, the optimized and runtime configurations achieve higher end-to-end phased recall than the default configuration, showing that the selected tuning generalizes beyond the hard optimization scenario.

#### Key observations

- **O1** The tuned configurations generalize across scenarios: In all four scenarios, both `optimized` and `runtime` outperform `default` on called effective phased recall and called shared heterozygous recall.
- **O2** The optimized configuration is the safest general recommendation: `optimized` is never worse than `default` and is either best or tied for best on the main called metrics across all scenarios.
- **O3** The runtime-biased configuration is highly competitive: `runtime` is effectively identical to `optimized` on the main called metrics and remains competitive on runtime.
- **O4** The gains transfer to both easy and hard cases: The tuned configurations improve performance not only in the hard scenario, but also in the baseline, dropout, and interaction scenarios.
- **O5** The performance difference between the `default` configuration, and the optimized (`optimized` and `runtime`) configurations, is smaller as the phasing difficulty increases.

#### Takeaway
The recommended tuned settings are robust across the representative scenarios tested in this study. By loosening the filtering thresholds, inaccurate variant site information is passed into WhatsHap in harsher scenarios, thus reducing the phasing performance edge of the optimized configurations. The optimzied configuration is the cleanest general recommendation because it consistently improves end-to-end phased recall and reduces fragmentation relative to the default.

#### 10.5.5 Runtime breakdown and DP-scaling interpretation

A small runtime breakdown of the called phasing stage was examined for representative hard-scenario configurations:

| Configuration | Build ReadSet (s) | Read Selection (s) | Solve (s) | Called Phasing Total (s) |
|---|---:|---:|---:|---:|
| Default | 0.3734 | 0.0023 | 0.0610 | 0.4433 |
| optimized | 0.3679 | 0.0061 | 0.0056 | 0.3854 |
| Runtime (efficiency-oriented) | 0.3685 | 0.0050 | 0.0030 | 0.3827 |

This breakdown shows that, in the current research pipeline, the largest share of called phasing runtime is spent in readset construction, while read selection and the downstream solve stage contribute smaller but configuration-dependent fractions. This result should be interpreted cautiously, because the present study uses a custom adaptor rather than the standard production-facing WhatsHap front-end.

A separate SNP-density scaling study compared `default` and `optimized` under the same hard scenario while varying `num_snps ∈ {400, 800, 1200, 1600}`.

#### Results summary

- **optimized**
  - `called_time_solve_sec`: `0.0026 → 0.0046 → 0.0051 → 0.0067`
  - `oracle_time_solve_sec`: `0.0044 → 0.0066 → 0.0091 → 0.0114`
  - `called_effective_phased_recall`: `0.3071 → 0.3705 → 0.3295 → 0.3280`
  - `called_num_phase_sets`: `6.0 → 5.7 → 7.3 → 6.3`
- **Default**
  - `called_time_solve_sec`: `0.0053 → 0.0352 → 0.0545 → 0.1002`
  - `oracle_time_solve_sec`: `0.0114 → 0.0533 → 0.0863 → 0.1407`
  - `called_effective_phased_recall`: `0.2935 → 0.3539 → 0.3164 → 0.3141`
  - `called_num_phase_sets`: `8.3 → 10.7 → 10.0 → 7.7`

These trends are reflected by Figures 10.5.5.1–10.5.5.4, especially Figure 10.5.5.1 for the steeper default solve-time growth and Figure 10.5.5.4 for the corresponding fragmentation difference.

\begin{center}
\includegraphics[width=0.82\linewidth]{figures/10_experiments/fig_10_5_5_1_called_time_solve_sec_dp_scaling.png}
\end{center}

**Figure 10.5.5.1 — Called solve time vs SNP count under default and optimized configurations.**
Called solve time increases with SNP density under both configurations, but grows much more steeply under the default setting. This shows that the tuned configuration makes the residual phasing problem easier to solve as variant density increases.

#### Key observations

- **O1** Solve time increases with SNP density: In both oracle and called regimes, solve time increases as `num_snps` increases from 400 to 1600.
- **O2** The optimized configuration scales much better than default: The growth in solve time is much steeper under the default configuration than under the optimized configuration.
- **O3** Optimized remains better on phasing continuity as the problem grows: Across the scaling range, the optimized configuration generally maintains fewer phase sets and slightly better phased recall than the default configuration.
- **O4** Solve time still does not dominate total runtime end-to-end: Although solve time grows with SNP density, total phasing runtime remains dominated by preprocessing / readset construction in the current research pipeline.
- **O5** Computational workload is transferred to other processes: Despite the significant improvement in runtime scalability, there is a slight increase in runtime for read selection.

#### Takeaway
Optimization improves not only final phased output but also the difficulty of the residual phasing problem faced by the solver. In the current research pipeline, end-to-end runtime remains dominated by preprocessing due to a readset construction adapter, but the tuned configuration substantially improves solve-stage scaling as SNP density increases. However, it is observed that, by loosening the filtering thresholds, some of the computational worklaod is transferred from the solver, to other processes, such as read selection and very possibly variant calling. In order to pin down the exact end-to-end runtime improvement from the optimized configurations, further analysis on runtime profiling in other processes, including alignment and variant calling, need to be done.

---

### 10.6 Cross-experiment synthesis

#### 10.6.1 Attribution summary (oracle vs called)

Taken together, the oracle-vs-called comparison provides a consistent attribution picture across the stress studies.

- **Duplications:** primarily a called-regime phasing / evidence-quality stressor rather than a strong caller-limited bottleneck
- **Coverage dropout:** primarily a fragmentation-dominated stressor, becoming mixed at high severity
- **Correlated error bursts:** a weak and somewhat noisy called-regime correctness stressor under the current parameterization
- **Read length model:** a mixed continuity / overlap trade-off rather than a clear optimization knob
- **Truth indels with SNP-only policy:** primarily a methodology-validity result, not a performance stressor
- **Duplication × dropout interaction:** a mixed compounded weakness, strongest in the called regime

Overall, whenever oracle performance remains strong but called performance degrades, the dominant bottleneck lies in the overlap and quality of called heterozygous sites. When oracle performance also degrades, the stressor is directly harming phase-block connectivity or phasing completeness.

#### 10.6.2 Most impactful stressors

Across the isolated and interaction studies, the realism knobs differ substantially in practical impact.

- **Most impactful overall:** coverage dropout
- **Most impactful compounded weakness:** duplication × dropout interaction
- **Moderate impact:** strong duplication and the read-length trade-off
- **Weak / noisy impact:** correlated error bursts
- **Methodological rather than stress-inducing:** truth indels under SNP-only policy

If stressors are ranked by effect on the main end-to-end metric (`called_effective_phased_recall`), the most important findings are:

1. coverage dropout
2. duplication × dropout interaction
3. strong duplication / read-length trade-off effects
4. bursts (weak / noisy)

#### 10.6.3 Optimization summary

Optimization sweeps under the composite hard scenario showed that optimization headroom exists, but is uneven across knobs.

- **Useful caller-side knob:** `call_min_baseq`. Lowering the caller base-quality threshold from 20 to a moderate value such as 10 substantially improves calling recall, shared heterozygous overlap, and end-to-end phased recall without materially harming precision.
- **Useful phasing-side knob:** `min_baseq`, mainly as a rule-out of overly strict filtering. Phasing performance is stable for `min_baseq` = 0–15, but degrades at 20.
- **Useful runtime knob:** `max_coverage`. Increasing `max_coverage` beyond 8–10 does not improve phasing performance, while lower values such as 8 preserve full performance with lower runtime.
- **Weak / negligible knobs:** Caller `call_min_mapq`, phasing `min_mapq`, and high `max_coverage` values above the optimum.
- **Robustness of tuning:** The optimized configurations are applicable across baseline, dropout, interaction, and hard scenarios.
- **Algorithm-facing scaling result:** Solve stage is more expensive computationally due to increasing SNP density, however, the optimized configurations offers superior scalability when compared to the default configuration.

Based on the optimization sweeps, the recommended configuration for practical usage for this pipeline are:

- **Recommended `max_coverage`:** `8`
- **Recommended phasing thresholds:** `min_mapq = 20`, `min_baseq = 10`
- **Recommended calling thresholds:** `call_min_mapq = 20` (or leave at default, as the effect is negligible), `call_min_baseq = 10`

An alternative configuration which is more runtime-oriented, with `max_coverage` = 6 and other caller and phasing quality thresholds remain unchanged, maintains almost identical phasing accuracy. This alternative configuration maybe useful when efficiency is prioritized.

#### Final takeaway
The overall picture showed by the experiments is consistent. The recovery and preservation of useful heterozygous variant information is the primary limiting factor in end-to-end phasing performance, with good variant sites available, WhatsHap's phasing capability is near perfect. While coverage dropout is the strongest isolated realism stressor, the clearest compounded vulnerability is duplication combined with dropout. Instead of relying on tuning WhatsHap parameters, combining moderate variant-caller-side and phasing-side thresholds is the most effective practical optimization.