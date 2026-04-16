## 10. Experimental Evaluation

### 10.1 Experimental methodology and common setup

#### 10.1.1 Environment and toolchain

- OS / CPU: macOS / Apple M2 Max
- Python version: 3.12.12
- minimap2 version: 2.30-r1287
- samtools version: 1.23
- bcftools version: 1.23
- Vendored WhatsHap core path: `vendor/whatshap_core/whatshap/core.<platform-specific extension>`

#### 10.1.2 Pipeline regimes

- **Oracle regime:** phase using `*.oracle.vcf.gz` (phaser-limited)
- **Called regime:** phase using `*.called.vcf.gz` (end-to-end)
- **Indel policy:** when indels are enabled, use SNP-only phasing/evaluation (`phase_snps_only`, `eval_snps_only`)

For this project, a custom adaptor layer was used to construct the `ReadSet` consumed by the vendored WhatsHap core. This was necessary to support oracle-vs-called benchmarking, unified JSON provenance, and controlled SNP-only phasing/evaluation under indel-containing truth sets. The adaptor is therefore part of the project’s benchmarking infrastructure rather than a direct reproduction of the standard production-facing WhatsHap front-end.

#### 10.1.3 Fixed parameters (baseline constants)

Unless stated otherwise in a subsection, experiments use the baseline defaults defined by the experiment driver. These parameters cover the full pipeline from reference/truth generation to read simulation, alignment, calling, phasing, and evaluation.

**Reference generation**

- `ref_length`: length of the synthetic reference genome (bp). Larger references reduce boundary effects and allow realism stressors to be placed without overlap.
- `ref_preset`: reference complexity preset (e.g., plain vs realistic). Presets control the presence of hard regions and are recorded in `*.ref.meta.json`.
- `dup_segments`: number of duplicated segments inserted into the reference. Higher values increase mapping ambiguity and can reduce calling recall or increase fragmentation.
- `dup_len`: length (bp) of each duplicated segment.
- `dup_min_gap`: minimum separation (bp) between duplicated placements.

**Truth generation**

- `num_snps`: number of SNP sites inserted into the truth set; controls variant density.
- `het_rate`: fraction of truth SNPs that are heterozygous. Only heterozygous sites contribute to phasing difficulty.
- `phased_truth`: whether the truth VCF uses phased genotypes (e.g., `0|1`); required for switch-error-based evaluation.
- `random_phase`: whether truth phase orientation is randomized to avoid bias in haplotype assignment.
- Indel parameters (when enabled):
  - `num_indels`: number of truth indels
  - `indel_min_len`, `indel_max_len`: indel length range (bp)
  - `indel_het_rate`: fraction of indels that are heterozygous

**Read simulation**

- `platform`: read simulation platform (e.g., `ont`).
- `ont_profile`: ONT error/quality preset (e.g., `q20`).
- `num_reads`: number of simulated reads; used here as a proxy for sequencing depth.
- `min_len`, `max_len`: read length range (bp); affects overlap connectivity and phase-block size.
- `hap1_frac`: fraction of reads sampled from haplotype 1; typically 0.5 for balanced coverage.
- Length model parameters:
  - `len_model`: read length distribution (`uniform` or `lognormal`)
  - `ln_mean`, `ln_sigma`: lognormal parameters when `len_model=lognormal`
- Start / coverage model parameters:
  - `start_model`: read start-position model (`uniform` or `dropout`)
  - `dropout_fraction`: fraction of the reference affected by coverage dropout
  - `dropout_block_len`: typical length of dropout blocks (bp)
- Correlated error burst parameters:
  - `burst_prob`: probability that a read contains elevated-error bursts
  - `burst_count`: number of bursts per read
  - `burst_len`: burst length (bp)
  - `burst_mult`: multiplier applied to base error rates within bursts

**Alignment and calling thresholds**

- `map_preset`: minimap2 mapping preset (default: `map-ont`).
- `call_min_mapq`: minimum read mapping quality used in variant calling (`bcftools mpileup -q`).
- `call_min_baseq`: minimum base quality used in variant calling (`bcftools mpileup -Q`).

**Phasing parameters (WhatsHap integration)**

- `vcf_source`: which VCF to phase:
  - `oracle`: truth-derived variant set
  - `called`: called variant set
  - `both`: run both for attribution
- `max_coverage`: WhatsHap read-selection cap.
- `min_mapq`: minimum mapping quality for reads used in phasing.
- `min_baseq`: minimum base quality for allele observations used in phasing.
- Indel-safe evaluation policy flags:
  - `phase_snps_only`: restrict phasing input to biallelic SNPs
  - `eval_snps_only`: restrict evaluation truth to biallelic SNPs

These flags prevent indel representation differences from corrupting SNP phasing evaluation.

#### 10.1.4 Metrics reported (from `aggregate.csv`)

Metrics are computed per run from `*.pipeline.json` and `*.eval.json`, then aggregated across seeds into `aggregate.csv`. Both oracle (phaser-focused) and called (end-to-end) metrics are reported to attribute performance loss to calling vs phasing.

**Calling quality (called regime only)**

- `call_precision`: fraction of called SNPs that match truth SNP sites. Intuition: “How many called variants are correct?”
- `call_recall`: fraction of truth SNPs recovered by the caller. Intuition: “How many true variants were found?”
- `truth_snps`: number of truth SNPs after any SNP-only filtering.
- `called_snps`: number of called SNPs.
- `shared_snps`: number of SNPs present in both truth and called sets.

**Phasing coverage / completeness**

- `*_shared_het_recall` (called): fraction of truth heterozygous SNPs that survive into the phasing input overlap.
- `*_phasing_rate_shared_het` (called): fraction of shared heterozygous SNPs that are actually phased.
- `*_effective_phased_recall` (oracle/called): fraction of truth heterozygous SNPs that are both present in the evaluation set and correctly phased under best block flip. This is the main end-to-end usable-phasing metric.

**Phasing correctness**

- `*_phase_accuracy` (oracle/called): accuracy of phased heterozygous genotypes under an optimal global flip per phase set.
- `*_switch_error` (oracle/called): switch error rate, computed as `#switches / #adjacent_het_pairs_compared`. This is often the most informative correctness metric for long-read phasing.

**Phase-block fragmentation**

- `*_num_phase_sets` (oracle/called): number of distinct phase sets (PS) in the phased VCF.
- Phase-set size distribution (from eval JSON): per-phase-set sizes and best-flip matches, used to interpret whether fragmentation is driven by coverage gaps, ambiguous mapping, or filtering thresholds.

**Runtime / efficiency**

- `time_total_sec`: total pipeline runtime per run, when recorded.
- `called_time_*`: phasing-stage timing breakdown for called runs, used in the optimization discussion.

**Oracle vs called interpretation**

- Oracle metrics represent phaser-limited performance under correct variant sites.
- Called metrics represent end-to-end performance including variant-calling limitations.
- The gap between oracle and called effective phased recall indicates how much performance is lost due to variant overlap and callset quality rather than phasing itself.

---

### 10.2 Baseline: depth scaling (coverage proxy via num_reads)

#### Aim
Establish baseline performance curves as sequencing depth increases, and attribute performance limitations by comparing the oracle (phaser-limited) regime against the called (end-to-end) regime.

#### Setup

- Varied: `num_reads ∈ {50, 100, 200, 400}` (3 seeds each)
- Fixed: ONT-like simulation profile (`q20`), reference length 80 kb, 800 truth SNPs (`het_rate = 0.8`), read length range 2–6 kb, and the calling/phasing thresholds defined in §10.1.

**Figure 10.2.1 — Variant calling recall vs sequencing depth.**  
Calling recall is computed as `shared_snps / truth_snps`, where `shared_snps` are SNP records present in both the truth VCF and the called VCF. As read depth increases, recall rises sharply, indicating that variant discovery is the main limiter at low coverage. Points show mean over seeds; error bars show standard deviation.

**Figure 10.2.2 — Oracle effective phased recall vs sequencing depth (phasing-only upper bound).**  
Oracle effective phased recall measures phasing performance when the input variants are truth-derived SNP sites, removing variant calling errors from the pipeline. The curve rises with depth and approaches a high plateau, representing phasing-only capability under the simulated ONT Q20 model. Points show mean over seeds; error bars show standard deviation.

\begin{center}
\includegraphics[width=0.82\linewidth]{figures/10_experiments/fig_10_2_3_called_effective_phased_recall.png}
\end{center}

**Figure 10.2.3 — Called effective phased recall vs sequencing depth (end-to-end performance).**  
Called effective phased recall measures the fraction of truth heterozygous SNPs that end up both present in the called set and correctly phased (after best block flip). This curve is substantially lower than the oracle curve at low depth, showing that variant calling recall and callset composition dominate end-to-end phasing performance early on. Points show mean over seeds; error bars show standard deviation.

**Figure 10.2.4 — Oracle number of phase sets vs sequencing depth (block fragmentation, phasing-only).**  
The number of phase sets (PS blocks) reflects fragmentation: more phase sets means shorter blocks and weaker connectivity. Under oracle variants, increasing depth reduces fragmentation as reads bridge more heterozygous sites. Points show mean over seeds; error bars show standard deviation.

**Figure 10.2.5 — Called number of phase sets vs sequencing depth (block fragmentation, end-to-end).**  
The phase-set count using called variants reflects fragmentation under end-to-end conditions. Unlike the oracle curve, this may be non-monotonic at moderate depths because the called VCF itself changes with depth, affecting graph connectivity. Points show mean over seeds; error bars show standard deviation.

**Figure 10.2.6 — Called switch error rate vs sequencing depth.**  
Switch error rate is computed over adjacent phased heterozygous SNP pairs within phase sets after choosing the best global flip per phase set. The mean generally decreases as depth increases, but can show variance at intermediate depths due to small denominators and changing callsets. Error bars show standard deviation; values should be interpreted as bounded to [0, 1] even if `mean ± std` extends below 0.

#### Results summary (means over seeds)

Calling (called regime):

- `call_recall`: 0.224 (50) → 0.399 (100) → 0.662 (200) → 0.877 (400)
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

Most loss in called effective phased recall therefore arises from reduced shared heterozygous site recall rather than poor phase accuracy, indicating that variant recovery and connectivity are the dominant bottlenecks in the baseline setting.

#### Key observations

- O1 (Calling is depth-limited): Calling recall increases strongly with `num_reads` (0.224 → 0.877), indicating that the caller is mainly sensitivity-limited at low and moderate depth. Precision remains 1.0 in this controlled setup, so errors appear primarily as missed sites rather than false positives.
- O2 (WhatsHap is accurate given correct sites): In the oracle regime, switch error remains near zero and phase accuracy remains near 1 even at low depth. The main limitation at 50–100 reads is incompleteness: many heterozygous sites remain unphased because of connectivity gaps.
- O3 (Residual high-depth loss remains caller-limited): At 400 reads, called phasing is nearly perfect on shared heterozygous sites (`called_phasing_rate_shared_het ≈ 0.99`, accuracy ≈ 1.0). The remaining gap between oracle and called effective phased recall is therefore dominated by missing heterozygous sites in the called set (`called_shared_het_recall ≈ 0.85`).
- O4 (Mid-depth instability in the called regime): At 100–200 reads, the called regime shows occasional high switch error in individual seeds. This likely reflects both seed sensitivity and small switch-error denominators when the called VCF contains relatively few informative heterozygous sites.

#### Takeaway
The baseline depth sweep confirms that the pipeline behaves sensibly: calling recall increases with depth; oracle phasing becomes both accurate and contiguous; and end-to-end (called) phasing improves sharply but remains limited primarily by the overlap of correctly called heterozygous variants. At sufficiently high depth (400 reads in this setup), WhatsHap phasing is essentially correct on available het sites, implying that further end-to-end gains would require improved variant calling overlap rather than changes to the phasing algorithm itself.

---

### 10.3 Single-knob study: duplicated regions (reference stressor)

#### Aim
Measure how duplicated reference segments affect end-to-end phasing by increasing mapping ambiguity, and determine whether the main impact is on variant recovery, phasing continuity, or phasing correctness.

#### Setup
- Varied: `dup_segments ∈ {0, 1, 3, 5}` (3 seeds each)
- Fixed:
  - `dup_len = 3000`
  - `dup_min_gap = 500`
  - baseline depth fixed at `num_reads = 200`
  - ONT-like simulation profile (`ont_profile = q20`)
  - reference length `80 kb`
  - `800` truth SNPs (`het_rate = 0.8`)
  - read length range `2–6 kb`
  - `start_model = uniform`
  - no dropout (`dropout_fraction = 0.0`)
  - no correlated bursts (`burst_prob = 0.0`)
  - no truth indels (`num_indels = 0`)
  - phasing input source `vcf_source = both`

This experiment isolates duplicated sequence as a single realism stressor while keeping all other baseline conditions unchanged.

**Figure 10.3.1 — Variant calling recall vs duplicated-region count.**  
Calling recall remains nearly unchanged across duplication levels, indicating that the duplicated-region stressor does not substantially reduce SNP site recovery under the current alignment and calling settings.

**Figure 10.3.2 — Called effective phased recall vs duplicated-region count.**  
End-to-end phasing remains broadly stable at low-to-moderate duplication levels, but declines at the highest duplication setting. This indicates that duplicated regions become a substantial stressor in the called regime only once ambiguity is sufficiently strong.

**Figure 10.3.3 — Called switch error rate vs duplicated-region count.**  
Switch error remains modest at low-to-moderate duplication levels but increases markedly at the highest duplication setting, suggesting that duplicated regions can destabilize local phase relationships even when overall call recall remains similar.

**Figure 10.3.4 — Oracle effective phased recall vs duplicated-region count.**  
Oracle phasing performance remains high and stable across duplication levels, showing that WhatsHap itself is robust when correct variant sites are available and that the main duplication effect emerges in the end-to-end called regime.

**Figure 10.3.5 — Called number of phase sets vs duplicated-region count.**  
Phase-set count does not increase with duplication in this experiment and therefore should not be interpreted in isolation. At the highest duplication level, fewer phase sets coincide with worse phasing accuracy and higher switch error, indicating that phase-set count alone is not a sufficient measure of performance here.

#### Results summary (means over seeds)
- `call_recall`: 0.663 (0 dup) → 0.671 (1 dup) → 0.661 (3 dup) → 0.666 (5 dup)
- `call_precision`: 1.000 at all duplication levels in this synthetic setup
- `called_effective_phased_recall`: 0.436 → 0.473 → 0.479 → 0.372
- `called_shared_het_recall`: 0.590 → 0.599 → 0.586 → 0.593
- `called_phasing_rate_shared_het`: 0.819 → 0.836 → 0.848 → 0.778
- `called_phase_accuracy`: 0.890 → 0.920 → 0.960 → 0.802
- `called_switch_error`: 0.122 → 0.073 → 0.040 → 0.209
- `called_num_phase_sets`: 3.4 → 3.4 → 3.4 → 2.6
- `oracle_effective_phased_recall`: 0.925 → 0.932 → 0.924 → 0.918
- `oracle_switch_error`: near 0 across all duplication levels
- `oracle_num_phase_sets`: 1.8 → 1.8 → 1.2 → 1.8

Taken together, these results show that duplicated regions are a relatively mild stressor under the current parameterization. The strongest degradation appears only at the highest duplication level (`dup_segments = 5`) and is driven mainly by reduced called-regime phasing completeness and phase accuracy rather than by reduced callset overlap.

#### Key observations
- O1 (Calling recall is largely unaffected by duplication in this setup): Across `dup_segments ∈ {0,1,3,5}`, calling recall remains nearly flat (`0.663–0.671`) and call precision stays at 1.0. This indicates that the current duplication stressor does not substantially reduce SNP site recovery under the baseline alignment/calling settings.

- O2 (Oracle phasing remains stable): Oracle effective phased recall remains high and nearly unchanged (`0.918–0.932`), with near-zero oracle switch error throughout. This shows that when correct variant sites are supplied, WhatsHap itself remains robust to the duplicated-region stressor at these settings.

- O3 (Called-regime degradation appears mainly at the highest duplication level): At `dup_segments = 5`, called effective phased recall drops from `0.479` to `0.372`, while called switch error rises from `0.040` to `0.209`. However, `called_shared_het_recall` remains approximately constant (`~0.59`), indicating that the degradation is not primarily due to loss of shared heterozygous sites.

- O4 (Main loss is in phasing completeness/correctness on surviving sites): The highest-duplication condition reduces `called_phasing_rate_shared_het` (`0.848 → 0.778`) and `called_phase_accuracy` (`0.960 → 0.802`). This suggests that duplicated regions mainly weaken the consistency of phasing evidence in the called regime rather than simply reducing callset overlap.

- O5 (Phase-set count is not the main signal here): `called_num_phase_sets` does not increase with duplication, and is actually lower at `dup_segments = 5`. In this experiment, lower phase-set count does not imply better phasing, because it occurs alongside lower effective phased recall and higher switch error. Phase-set count therefore needs to be interpreted together with completeness and correctness metrics.

#### Takeaway
Under the current parameterization, duplicated regions are a relatively mild stressor for variant recovery and for oracle phasing. Their main impact appears only at the highest duplication level, where end-to-end (called) phasing becomes less complete and less accurate despite largely unchanged callset overlap. This suggests that duplicated regions primarily affect the quality and stability of phasing evidence in the called regime, rather than acting as a strong caller-limited bottleneck.

---

### 10.4 Single-knob study: coverage dropout (read start model)

#### Aim
Evaluate how localized coverage dropout affects phase-block continuity and effective phased recall, and determine whether the main impact is on variant recovery, phasing completeness, or both.

#### Setup

- Varied: `dropout_fraction ∈ {0.00, 0.05, 0.10, 0.20}`
- Fixed:
  - `dropout_block_len = 1000`
  - baseline depth fixed at `num_reads = 200`
  - ONT-like simulation profile (`ont_profile = q20`)
  - reference length `80 kb`
  - `800` truth SNPs (`het_rate = 0.8`)
  - read length range `2–6 kb`
  - `start_model = dropout`
  - no duplications (`dup_segments = 0`)
  - no correlated bursts (`burst_prob = 0.0`)
  - no truth indels (`num_indels = 0`)
  - phasing input source `vcf_source = both`

This experiment isolates coverage dropout as a single realism stressor while keeping the remaining baseline conditions unchanged.

**Figure 10.4.1 — Oracle number of phase sets vs dropout fraction.**  
The number of oracle phase sets increases strongly with dropout fraction, showing that localized coverage gaps fragment haplotype blocks even when correct variant sites are supplied to the phaser.

**Figure 10.4.2 — Called number of phase sets vs dropout fraction.**  
End-to-end phase fragmentation increases at low-to-moderate dropout levels, but phase-set count alone becomes harder to interpret at severe dropout because many heterozygous sites are no longer called or phased.

**Figure 10.4.3 — Oracle effective phased recall vs dropout fraction.**  
Oracle effective phased recall decreases as dropout increases, indicating that coverage gaps directly reduce phasing completeness by breaking read connectivity between heterozygous sites.

**Figure 10.4.4 — Called effective phased recall vs dropout fraction.**  
Called effective phased recall declines much more steeply than oracle effective phased recall, showing that dropout affects both phasing connectivity and, at stronger dropout levels, variant recovery in the called regime.

**Figure 10.4.5 — Variant calling recall vs dropout fraction.**  
Calling recall decreases gradually at moderate dropout and then sharply at the strongest dropout level, indicating that severe low-coverage regions substantially reduce the recovery of true SNP sites.

#### Results summary (means over seeds)

- `call_recall`: 0.663 (0.00) → 0.635 (0.05) → 0.605 (0.10) → 0.432 (0.20)
- `call_precision`: 1.000 at all dropout levels in this synthetic setup
- `oracle_effective_phased_recall`: 0.925 → 0.865 → 0.857 → 0.637
- `oracle_switch_error`: near 0 across all dropout levels
- `oracle_num_phase_sets`: 1.8 → 2.8 → 3.2 → 10.0
- `called_effective_phased_recall`: 0.436 → 0.359 → 0.269 → 0.081
- `called_shared_het_recall`: 0.590 → 0.565 → 0.529 → 0.362
- `called_phasing_rate_shared_het`: 0.819 → 0.704 → 0.514 → 0.239
- `called_phase_accuracy`: 0.890 → 0.894 → 0.981 → 0.922
- `called_switch_error`: 0.122 → 0.072 → 0.006 → 0.089
- `called_num_phase_sets`: 3.4 → 5.0 → 6.2 → 4.0

Decomposition of called effective phased recall shows that coverage dropout first increases loss through unphased shared heterozygous sites, and at the strongest dropout level also causes a substantial increase in loss from missing called sites.

#### Key observations

- O1 (Coverage dropout strongly increases fragmentation): In the oracle regime, the number of phase sets rises from 1.8 at `dropout_fraction = 0.00` to 10.0 at `0.20`, while oracle switch error remains near zero throughout. This shows that dropout primarily breaks phase-block continuity rather than causing incorrect phase assignments.

- O2 (Dropout directly reduces phasing completeness even with correct variant sites): Oracle effective phased recall declines from 0.925 to 0.637 as dropout increases. Because oracle shared-site recall remains fixed by construction, this reduction is attributable to weaker read connectivity and more unphased heterozygous sites.

- O3 (End-to-end performance degrades more strongly than oracle performance): Called effective phased recall falls from 0.436 to 0.081, a steeper decline than in the oracle regime. This indicates that dropout harms end-to-end phasing through both connectivity loss and reduced overlap of correctly called heterozygous variants.

- O4 (Moderate dropout is mainly a connectivity / completeness problem): From `dropout_fraction = 0.00` to `0.10`, the main worsening in the called regime is a drop in `called_phasing_rate_shared_het` (`0.819 → 0.514`) together with increased fragmentation. This indicates that moderate dropout primarily prevents surviving heterozygous sites from being phased into long connected blocks.

- O5 (Severe dropout becomes a mixed calling + phasing bottleneck): At `dropout_fraction = 0.20`, `call_recall` drops sharply to 0.432 and `called_shared_het_recall` falls to 0.362, while `called_phasing_rate_shared_het` also falls to 0.239. This shows that strong dropout not only fragments phase blocks but also substantially reduces variant recovery in the called regime.

- O6 (Switch error is not the main signal in this study): Although called switch error does not increase monotonically, this should not be interpreted as improved phasing at moderate dropout. As dropout increases, fewer heterozygous sites remain phased, making switch error less informative than completeness and fragmentation metrics.

#### Takeaway
Coverage dropout is a strong fragmentation-dominated stressor in this pipeline. Its primary effect is to create low-coverage gaps that prevent reads from connecting heterozygous sites into long phase blocks, which reduces phasing completeness even in the oracle regime. In the called regime, this effect is amplified, and at severe dropout levels the stressor becomes mixed: both phase connectivity and variant recovery deteriorate substantially. Overall, coverage dropout is one of the clearest realism knobs in this study for directly disrupting haplotype continuity.

---

### 10.5 Single-knob study: correlated error bursts

#### Aim
Assess how localized high-error segments within reads affect variant calling and phasing stability, and determine whether the main impact is on callset recovery, phasing correctness, or phasing completeness.

#### Setup
- Varied: `burst_prob ∈ {0.0, 0.3, 0.6}`
- Fixed:
  - `burst_count = 2`
  - `burst_len = 300`
  - `burst_mult = 8.0`
  - baseline depth fixed at `num_reads = 200`
  - ONT-like simulation profile (`ont_profile = q20`)
  - reference length `80 kb`
  - `800` truth SNPs (`het_rate = 0.8`)
  - read length range `2–6 kb`
  - `start_model = uniform`
  - no duplications (`dup_segments = 0`)
  - no coverage dropout (`dropout_fraction = 0.0`)
  - no truth indels (`num_indels = 0`)
  - phasing input source `vcf_source = both`

This experiment isolates correlated error bursts as a single realism stressor while keeping the remaining baseline conditions unchanged.

#### Figures to include
- **Figure 10.5.1 — Called switch error rate vs burst probability.**
- **Figure 10.5.2 — Called effective phased recall vs burst probability.**
- **Figure 10.5.3 — Variant calling recall vs burst probability.**
- **Figure 10.5.4 — Called number of phase sets vs burst probability.**
- **Figure 10.5.5 — Oracle effective phased recall vs burst probability.**

**Figure 10.5.1 — Called switch error rate vs burst probability.**  
Called switch error shows a modest increase at the intermediate burst setting, suggesting that correlated high-error segments can destabilize local phase relationships in the end-to-end called regime. However, the effect is not monotonic across the burst sweep.

**Figure 10.5.2 — Called effective phased recall vs burst probability.**  
End-to-end phased recall decreases at `burst_prob = 0.3` but rebounds at `0.6`, indicating that the current burst stressor produces a weak and somewhat noisy effect rather than a strong monotonic degradation trend.

**Figure 10.5.3 — Variant calling recall vs burst probability.**  
Calling recall remains broadly stable across burst settings, indicating that correlated error bursts do not substantially reduce SNP site recovery under the current parameterization.

**Figure 10.5.4 — Called number of phase sets vs burst probability.**  
Phase-set count does not increase with burst probability and therefore is not the main signal in this study. Under the current settings, correlated error bursts appear to affect phasing correctness more than fragmentation.

**Figure 10.5.5 — Oracle effective phased recall vs burst probability.**  
Oracle phasing performance remains high and stable across burst settings, showing that WhatsHap itself is robust when correct variant sites are supplied and that the observed burst effect arises mainly in the end-to-end called regime.

#### Results summary (means over seeds)
- `call_recall`: 0.663 (`burst_prob = 0.0`) → 0.681 (`0.3`) → 0.676 (`0.6`)
- `call_precision`: 1.000 at all burst levels in this synthetic setup
- `oracle_effective_phased_recall`: 0.925 → 0.926 → 0.921
- `oracle_switch_error`: near 0 across all burst levels
- `oracle_num_phase_sets`: 1.8 → 1.0 → 1.6
- `called_effective_phased_recall`: 0.436 → 0.369 → 0.455
- `called_shared_het_recall`: 0.590 → 0.611 → 0.608
- `called_phasing_rate_shared_het`: 0.819 → 0.790 → 0.837
- `called_phase_accuracy`: 0.890 → 0.766 → 0.887
- `called_switch_error`: 0.122 → 0.142 → 0.108
- `called_num_phase_sets`: 3.4 → 2.0 → 2.6

Decomposition of called effective phased recall shows that the main additional loss at `burst_prob = 0.3` comes from switch / phase error rather than from reduced shared-site recall. However, this effect is not monotonic, because the `burst_prob = 0.6` condition rebounds close to baseline on most metrics.

#### Key observations
- O1 (Calling recall is largely unaffected by bursty errors in this setup): Across `burst_prob ∈ {0.0, 0.3, 0.6}`, calling recall remains nearly flat (`0.663–0.681`) and call precision stays at 1.0. This indicates that the current burst stressor does not substantially reduce SNP site recovery under the baseline alignment and calling settings.

- O2 (Oracle phasing remains stable): Oracle effective phased recall stays near `0.92–0.93`, oracle switch error remains near zero, and oracle phase-set counts stay low. This shows that WhatsHap itself is robust to the current burst stressor when correct variant sites are supplied.

- O3 (The main visible effect is a dip in called-regime phasing correctness at intermediate burst probability): At `burst_prob = 0.3`, called effective phased recall drops from `0.436` to `0.369`, while called phase accuracy decreases from `0.890` to `0.766` and called switch error increases from `0.122` to `0.142`. This indicates that bursty errors can reduce the stability of phasing decisions in the end-to-end called regime.

- O4 (The burst effect is not monotonic under the current parameterization): At `burst_prob = 0.6`, called effective phased recall rebounds to `0.455`, while called phase accuracy and switch error return close to baseline. This suggests that the current burst sweep produces a weak and somewhat noisy signal rather than a strong monotonic degradation trend.

- O5 (The main additional loss at `burst_prob = 0.3` comes from correctness rather than site recovery): Decomposition shows that the largest extra loss term at `burst_prob = 0.3` is switch / phase error, whereas `called_shared_het_recall` remains stable (`~0.61`). This suggests that correlated error bursts mainly affect the quality of phasing evidence on surviving called sites rather than substantially reducing callset overlap.

- O6 (Fragmentation is not the dominant signal in this study): Called phase-set counts do not increase with burst probability and are lower than baseline in both burst settings. Therefore, correlated error bursts in the current setup are better characterized as a weak phasing-correctness stressor than as a fragmentation stressor.

#### Takeaway
Under the current parameterization, correlated error bursts are a relatively weak and noisy stressor in this pipeline. They do not substantially reduce calling recall or oracle phasing performance, indicating that neither variant recovery nor phasing with correct sites is strongly disrupted. The main visible effect is a reduction in called-regime phasing correctness at intermediate burst probability, suggesting that structured local errors can sometimes destabilize end-to-end phasing decisions. However, because this effect is not monotonic across the burst sweep, the present results should be interpreted as preliminary rather than as evidence of a strong burst-driven failure mode.

---

### 10.6 Single-knob study: read length model

#### Aim
Compare uniform and lognormal read length distributions to assess how read-length variability affects phasing connectivity, fragmentation, and end-to-end phased recall.

#### Setup
- Varied: `len_model ∈ {uniform, lognormal}`
- Fixed:
  - baseline depth fixed at `num_reads = 200`
  - ONT-like simulation profile (`ont_profile = q20`)
  - reference length `80 kb`
  - `800` truth SNPs (`het_rate = 0.8`)
  - read length range `2–6 kb` for the uniform baseline
  - lognormal parameters fixed at `ln_mean = 8.3`, `ln_sigma = 0.6`
  - `start_model = uniform`
  - no duplications (`dup_segments = 0`)
  - no coverage dropout (`dropout_fraction = 0.0`)
  - no correlated bursts (`burst_prob = 0.0`)
  - no truth indels (`num_indels = 0`)
  - phasing input source `vcf_source = both`

This experiment isolates the read length distribution as a single realism stressor while keeping the remaining baseline conditions unchanged.

#### Figures to include
- **Figure 10.6.1 — Oracle number of phase sets by read length model.**
- **Figure 10.6.2 — Called number of phase sets by read length model.**
- **Figure 10.6.3 — Oracle effective phased recall by read length model.**
- **Figure 10.6.4 — Called effective phased recall by read length model.**
- **Figure 10.6.5 — Variant calling recall by read length model.**

**Figure 10.6.1 — Oracle number of phase sets by read length model.**  
The lognormal read length model reduces the number of oracle phase sets relative to the uniform baseline, indicating improved long-range connectivity between heterozygous sites when correct variant sites are supplied.

**Figure 10.6.2 — Called number of phase sets by read length model.**  
The lognormal read length model substantially reduces called phase-set count, showing that it improves phase-block continuity under end-to-end conditions. However, this continuity gain does not translate into a large increase in effective phased recall.

**Figure 10.6.3 — Oracle effective phased recall by read length model.**  
Oracle effective phased recall is slightly higher under the lognormal model, suggesting that the alternative read length distribution modestly improves phasing completeness when correct variant sites are available.

**Figure 10.6.4 — Called effective phased recall by read length model.**  
Called effective phased recall changes only slightly between the uniform and lognormal models, indicating that improved connectivity and overlap are largely offset by somewhat noisier called-regime phasing correctness.

**Figure 10.6.5 — Variant calling recall by read length model.**  
Calling recall is modestly higher under the lognormal read length model, indicating that the alternative read length distribution slightly improves recovery of true SNP sites in the called regime.

#### Results summary (means over seeds)
- `call_recall`: 0.663 (uniform) → 0.701 (lognormal)
- `call_precision`: 1.000 under both read length models in this synthetic setup
- `oracle_effective_phased_recall`: 0.925 → 0.931
- `oracle_switch_error`: near 0 under both models
- `oracle_num_phase_sets`: 1.8 → 1.4
- `called_effective_phased_recall`: 0.436 → 0.443
- `called_shared_het_recall`: 0.590 → 0.633
- `called_phasing_rate_shared_het`: 0.819 → 0.810
- `called_phase_accuracy`: 0.890 → 0.855
- `called_switch_error`: 0.122 → 0.165
- `called_num_phase_sets`: 3.4 → 1.6

Taken together, these results show that the lognormal read length model improves callset overlap and phase-block continuity, but these gains translate into only a very small increase in called effective phased recall because called-regime phasing correctness becomes slightly worse and more variable across seeds.

#### Key observations
- O1 (Lognormal read lengths improve callset overlap): Compared with the uniform baseline, the lognormal model increases `call_recall` (`0.663 → 0.701`) and `called_shared_het_recall` (`0.590 → 0.633`). This suggests that the alternative read length distribution provides slightly stronger variant-support coverage under the current settings.

- O2 (Lognormal read lengths improve phase-block continuity): Oracle phase-set count decreases from `1.8` to `1.4`, and called phase-set count decreases more strongly from `3.4` to `1.6`. Oracle effective phased recall also increases slightly (`0.925 → 0.931`). Together, these results indicate improved long-range connectivity between heterozygous sites under the lognormal read length model.

- O3 (End-to-end improvement is small despite better connectivity): Although overlap and continuity improve, called effective phased recall changes only slightly (`0.436 → 0.443`). This indicates that better block connectivity alone is not sufficient to produce a large end-to-end gain in this setting.

- O4 (Completeness gains are offset by noisier called-regime correctness): Under the lognormal model, `called_phasing_rate_shared_het` decreases slightly (`0.819 → 0.810`), `called_phase_accuracy` decreases (`0.890 → 0.855`), and `called_switch_error` increases (`0.122 → 0.165`). This suggests that while the lognormal distribution helps connect more sites, the resulting called-regime phasing is somewhat less stable across seeds.

- O5 (The read length model is therefore a mixed, not clearly beneficial, knob): The current comparison shows that read-length distribution can improve continuity-related metrics without producing a correspondingly strong improvement in end-to-end usable phasing output. This makes the read length model informative for mechanism analysis, but not yet a clear optimization lever under the present parameterization.

#### Takeaway
Under the current parameterization, the read length model is a mixed stressor rather than a clear optimization knob. The lognormal model improves callset overlap and phase-block continuity, indicating better long-range connectivity, but these benefits are largely offset by slightly worse and noisier called-regime phasing correctness. As a result, end-to-end effective phased recall changes only marginally. This suggests that changing the read length distribution alone is not sufficient to deliver a strong practical improvement in this pipeline, although it does reveal an important trade-off between connectivity and phasing stability.

---

### 10.7 Truth indels and SNP-only policy (evaluation robustness)

#### Aim
Show that introducing indels into the truth set does not invalidate SNP phasing evaluation when SNP-only phasing and evaluation are enforced. This section verifies that the benchmarking framework remains robust to indel-containing truth without conflating SNP phasing performance with indel representation mismatches.

#### Setup
- Varied: `num_indels ∈ {0, 80, 200}`
- Fixed:
  - `indel_min_len = 1`
  - `indel_max_len = 5`
  - `indel_het_rate = 0.5`
  - baseline depth fixed at `num_reads = 200`
  - ONT-like simulation profile (`ont_profile = q20`)
  - reference length `80 kb`
  - `800` truth SNPs (`het_rate = 0.8`)
  - read length range `2–6 kb`
  - `start_model = uniform`
  - no duplications (`dup_segments = 0`)
  - no coverage dropout (`dropout_fraction = 0.0`)
  - no correlated bursts (`burst_prob = 0.0`)
  - phasing input source `vcf_source = both`
  - SNP-only safeguards enabled automatically when indels are present:
    - `phase_snps_only = true`
    - `eval_snps_only = true`

This experiment isolates the effect of adding indels to truth generation while preserving SNP-only phasing and evaluation comparability.

#### Figures to include
- **Figure 10.7.1 — Oracle effective phased recall vs number of truth indels.**
- **Figure 10.7.2 — Called effective phased recall vs number of truth indels.**
- **Figure 10.7.3 — Variant calling recall vs number of truth indels.**
- **Figure 10.7.4 — Called number of phase sets vs number of truth indels.**

**Figure 10.7.1 — Oracle effective phased recall vs number of truth indels.**  
Oracle SNP-phasing performance remains essentially unchanged as truth indels are introduced, showing that SNP-only phasing and evaluation preserve a stable phasing-only benchmark in indel-containing truth sets.

**Figure 10.7.2 — Called effective phased recall vs number of truth indels.**  
Called SNP-phasing performance remains broadly stable across indel settings under SNP-only evaluation, indicating that truth indels do not invalidate end-to-end SNP phasing assessment when the appropriate safeguards are enabled.

**Figure 10.7.3 — Variant calling recall vs number of truth indels.**  
Calling recall changes only modestly across indel settings, suggesting that the presence of truth indels has limited indirect impact on SNP recovery under the current pipeline and SNP-only evaluation policy.

**Figure 10.7.4 — Called number of phase sets vs number of truth indels.**  
Called phase-set count varies across indel settings but does not show evidence of catastrophic fragmentation or evaluation instability, supporting the conclusion that SNP-only safeguards preserve interpretable SNP phasing metrics in indel-containing truth sets.

#### Results summary (means over seeds)
- `call_recall`: 0.663 (`num_indels = 0`) → 0.691 (`80`) → 0.654 (`200`)
- `call_precision`: 1.000 → 0.999 → 0.997
- `oracle_effective_phased_recall`: 0.925 → 0.925 → 0.921
- `oracle_switch_error`: near 0 across all indel settings
- `oracle_num_phase_sets`: 1.8 → 1.2 → 1.8
- `called_effective_phased_recall`: 0.436 → 0.573 → 0.493
- `called_shared_het_recall`: 0.590 → 0.623 → 0.577
- `called_phasing_rate_shared_het`: 0.819 → 0.919 → 0.856
- `called_phase_accuracy`: 0.890 → 0.999 → 1.000
- `called_switch_error`: 0.122 → 0.001 → 0.000
- `called_num_phase_sets`: 3.4 → 2.4 → 4.0

Taken together, these results show that introducing truth indels does not destabilize SNP phasing evaluation when `phase_snps_only = true` and `eval_snps_only = true` are enabled. Oracle SNP-phasing metrics remain essentially unchanged, and called SNP metrics remain within a comparable range across indel settings.

#### Key observations
- O1 (Oracle SNP-phasing remains stable under truth indels): Oracle effective phased recall stays essentially unchanged (`0.925 → 0.925 → 0.921`), oracle switch error remains near zero, and oracle phase-set counts remain low. This shows that adding indels to truth generation does not corrupt SNP-only phasing evaluation when SNP-only safeguards are enabled.

- O2 (Called SNP-phasing does not collapse in the presence of indels): Called effective phased recall remains within a similar overall range across `num_indels ∈ {0, 80, 200}`, and calling recall stays relatively stable (`0.663–0.691`, then `0.654` at the strongest indel setting). This indicates that the benchmarking pipeline continues to produce interpretable SNP-level results even when truth indels are present.

- O3 (The main result is evaluation robustness, not a strong indel-driven degradation trend): The `80`-indel condition performs slightly better than the no-indel baseline on several called metrics, including called effective phased recall and switch error. This should not be interpreted as evidence that indels improve SNP phasing; rather, it suggests that under SNP-only filtering the remaining differences are dominated by ordinary stochastic variation and indirect alignment/calling effects rather than by evaluation artifacts.

- O4 (Even stronger indel conditions remain manageable under SNP-only policy): At `num_indels = 200`, oracle metrics remain near baseline and called effective phased recall remains above the no-indel baseline, while called switch error is effectively zero. This supports the claim that the SNP-only policy successfully prevents indel representation differences from invalidating SNP phasing evaluation.

- O5 (SNP-only safeguards are therefore justified methodologically): Overall, the indel suite supports the use of `phase_snps_only` and `eval_snps_only` as robust safeguards whenever truth indels are present, because they preserve stable and interpretable SNP phasing metrics without obvious representation-induced distortion.

#### Takeaway
The SNP-only policy successfully preserves interpretable SNP phasing evaluation in the presence of truth indels. Under the current settings, introducing indels into truth generation does not materially destabilize oracle or called SNP-level phasing metrics, indicating that indel representation mismatch has been effectively controlled. This validates the use of `phase_snps_only = true` and `eval_snps_only = true` as an evaluation safeguard for indel-containing truth sets in the remainder of the study.

---

### 10.8 Interaction study: duplications × dropout

#### Aim
Test whether duplicated regions and coverage dropout compound non-linearly, and determine whether the resulting weakness is mainly caller-limited, connectivity-limited, or mixed.

#### Setup
Small 2×2 interaction grid:
- `dup_segments ∈ {0, 5}`
- `dropout_fraction ∈ {0.0, 0.1}`

Fixed:
- `dropout_block_len = 1000`
- `dup_len = 3000`
- `dup_min_gap = 500`
- baseline depth fixed at `num_reads = 200`
- ONT-like simulation profile (`ont_profile = q20`)
- reference length `80 kb`
- `800` truth SNPs (`het_rate = 0.8`)
- read length range `2–6 kb`
- `start_model = dropout`
- no correlated bursts (`burst_prob = 0.0`)
- no truth indels (`num_indels = 0`)
- phasing input source `vcf_source = both`

This experiment tests whether two realistic stressors that were individually interpretable in earlier sections combine into a more severe end-to-end phasing failure mode.

#### Figures to include
- **Figure 10.8.1 — Variant calling recall across duplication × dropout conditions.**
- **Figure 10.8.2 — Called effective phased recall across duplication × dropout conditions.**
- **Figure 10.8.3 — Called number of phase sets across duplication × dropout conditions.**
- **Figure 10.8.4 — Called switch error across duplication × dropout conditions.** (secondary; interpret with caution)

**Figure 10.8.1 — Variant calling recall across duplication × dropout conditions.**  
Calling recall is only mildly affected by duplication alone, but decreases further when duplication is combined with coverage dropout. This shows that ambiguous mapping compounds the loss of callable sites under low-coverage conditions.

\begin{center}
\includegraphics[width=0.82\linewidth]{figures/10_experiments/fig_10_8_2_called_effective_phased_recall.png}
\end{center}

**Figure 10.8.2 — Called effective phased recall across duplication × dropout conditions.**
The combined duplication + dropout condition produces the lowest end-to-end phased recall, showing that the two stressors interact to create a substantially more difficult phasing problem than either one alone.

**Figure 10.8.3 — Called number of phase sets across duplication × dropout conditions.**  
Phase fragmentation is driven primarily by dropout, but remains high in the combined condition. This indicates that low-coverage gaps are the main cause of broken block continuity, while duplication compounds the end-to-end difficulty in other ways.

**Figure 10.8.4 — Called switch error across duplication × dropout conditions.**  
Switch error should be interpreted cautiously in the interaction study. Although duplication alone produces the highest switch error, the combined condition is much more fragmented and leaves fewer phased adjacencies to compare, making completeness and fragmentation metrics more informative than switch error alone.

#### Results summary (means over seeds)
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
- `called_phase_accuracy`:
  - 0.890 → 0.981 → 0.802 → 0.940
- `called_switch_error`:
  - 0.122 → 0.006 → 0.207 → 0.062
- `called_num_phase_sets`:
  - 3.4 → 6.2 → 2.6 → 6.2

Taken together, these results show that the combined duplication + dropout condition is substantially worse than either single stressor alone, and that the strongest compounded loss occurs in the called regime rather than in the oracle regime.

#### Key observations
- O1 (The combined condition is the worst end-to-end setting): Called effective phased recall decreases from `0.436` in the baseline to `0.269` with dropout alone, `0.373` with duplication alone, and `0.197` when both stressors are combined. This shows that duplication and dropout together create a materially more difficult end-to-end phasing condition than either stressor in isolation.

- O2 (Dropout remains the main driver of phasing fragmentation): In the oracle regime, duplication alone has little effect (`oracle_num_phase_sets = 1.8`), whereas dropout increases fragmentation (`1.8 → 3.2`) and the combined condition increases it further (`4.0`). This indicates that coverage gaps are the primary cause of broken phase-block connectivity.

- O3 (The interaction is strongest in the called regime): Relative to dropout alone, the combined condition further reduces `call_recall` (`0.605 → 0.582`), `called_shared_het_recall` (`0.529 → 0.501`), and `called_phasing_rate_shared_het` (`0.514 → 0.410`). This shows that duplication compounds dropout by worsening both variant recovery and phasing completeness on surviving shared heterozygous sites.

- O4 (The compounded weakness is mixed rather than purely phaser-limited): Oracle effective phased recall decreases only modestly from `0.857` to `0.847` when duplication is added on top of dropout, but called effective phased recall decreases more substantially from `0.269` to `0.197`. Therefore, the main interaction effect is not purely within phasing itself; it arises more strongly from the end-to-end combination of reduced callable overlap and weaker phasing connectivity.

- O5 (Switch error is not the main interaction signal): Although duplication alone produces the highest called switch error, the combined condition does not. This should not be interpreted as an improvement, because the combined condition is much more fragmented and leaves fewer adjacent phased heterozygous pairs to compare. In this interaction study, completeness and fragmentation metrics are more informative than switch error alone.

#### Optimization implication
The duplication × dropout interaction suggests that the most important optimization opportunities under compounded stress are not limited to WhatsHap read selection. Preserving phasing evidence remains important (e.g., avoiding overly strict `min_baseq` filtering), but the interaction also indicates that upstream variant recovery in ambiguous and low-coverage regions is a major limiting factor. This motivates future tuning of calling-side thresholds and, potentially, more region-aware handling of ambiguous reads near coverage gaps.

#### Takeaway
Duplicated regions and coverage dropout interact to produce a more severe end-to-end phasing weakness than either stressor alone. Dropout drives fragmentation, while duplication further reduces callable overlap and phasing completeness in the called regime. This makes the interaction a useful diagnostic for optimization: it shows that practical improvements will likely require both careful preservation of phasing evidence and better upstream handling of difficult regions, rather than relying on a single WhatsHap parameter alone.

---

### 10.9 Optimization sweeps (WhatsHap / pipeline parameters)

#### Aim
Identify parameter settings that improve phasing performance under a composite hard scenario, and quantify trade-offs between phased recall, switch error, fragmentation, and robustness.

#### Hard scenario
Optimization sweeps were conducted on a fixed composite stress condition designed to combine several realism knobs:
- duplicated regions (`dup_segments = 5`)
- coverage dropout (`dropout_fraction = 0.1`, `dropout_block_len = 1000`)
- correlated error bursts (`burst_prob = 0.6`, `burst_count = 2`, `burst_len = 300`, `burst_mult = 8.0`)
- truth indels (`num_indels = 120`, SNP-only phasing/evaluation enabled)
- `ref_length = 120 kb`
- `num_snps = 1200`
- `num_reads = 300`

This hard scenario was chosen to test whether WhatsHap-side tuning can recover performance under combined realistic stresses rather than under isolated perturbations.

#### 10.9.1 Sweep A: `max_coverage`

##### Setup
- Varied: `max_coverage ∈ {10, 15, 25, 40}`
- Fixed: all other hard-scenario parameters

##### Figures to include
- **Figure 10.9.1.1 — Called effective phased recall vs `max_coverage`.**
- **Figure 10.9.1.2 — Called switch error vs `max_coverage`.**
- **Figure 10.9.1.3 — Called number of phase sets vs `max_coverage`.**
- **Figure 10.9.1.4 — Total pipeline runtime vs `max_coverage`.**

**Figure 10.9.1.1 — Called effective phased recall vs `max_coverage`.**  
Called effective phased recall remains essentially flat across the `max_coverage` sweep, indicating that increasing the WhatsHap read-selection cap does not provide a meaningful gain under the current hard scenario.

**Figure 10.9.1.2 — Called switch error vs `max_coverage`.**  
Called switch error shows no meaningful dependence on `max_coverage`, further indicating that read-selection depth is not the limiting factor for phasing correctness in the current hard scenario.

**Figure 10.9.1.3 — Called number of phase sets vs `max_coverage`.**  
Phase-set fragmentation remains unchanged across the `max_coverage` sweep, showing that retaining more reads does not improve haplotype continuity in this stress condition.

**Figure 10.9.1.4 — Total pipeline runtime vs `max_coverage`.**  
Runtime increases substantially as `max_coverage` increases, but phasing performance remains nearly unchanged. This indicates that higher `max_coverage` mainly adds compute cost rather than practical benefit.

##### Results summary
- `call_recall`: 0.5895 at all `max_coverage` settings
- `called_effective_phased_recall`: 0.2917 → 0.2919 → 0.2919 → 0.2919
- `called_switch_error`: 0.0056 at all `max_coverage` settings
- `called_num_phase_sets`: 9.6 at all `max_coverage` settings
- `called_shared_het_recall`: 0.5098 at all `max_coverage` settings
- `called_phasing_rate_shared_het`: 0.5860 → 0.5864 → 0.5864 → 0.5864
- `called_phase_accuracy`: 0.9748 at all `max_coverage` settings
- `oracle_effective_phased_recall`: 0.8256 → 0.8298 → 0.8300 → 0.8300
- `oracle_num_phase_sets`: 6.4 → 6.2 → 6.2 → 6.2
- `time_total_sec`: 3.76 → 4.37 → 6.30 → 6.39

##### Key observations
- O1 (Changing `max_coverage` has almost no effect on phasing performance): Called effective phased recall, switch error, phase-set count, and phasing-rate-on-shared-hets remain essentially unchanged across the full `max_coverage` sweep.

- O2 (Runtime increases without measurable accuracy gain): Total runtime increases substantially as `max_coverage` increases, but this additional compute does not translate into improved end-to-end phasing metrics.

- O3 (The read-selection cap is not the active bottleneck in this hard scenario): Although higher `max_coverage` retains more reads internally, the resulting phasing output is almost identical. This indicates that the limiting factor is not insufficient retained read coverage at the phasing stage.

- O4 (Optimization headroom from `max_coverage` is negligible): Under the current hard scenario, increasing the WhatsHap read-selection cap does not provide a practically useful optimization lever.

##### Takeaway
Adjusting `max_coverage` does not provide a meaningful optimization benefit under the composite hard scenario. Since phasing metrics remain flat while runtime increases, the default or lower `max_coverage` setting is preferable on efficiency grounds.

#### 10.9.2 Sweep B: phasing quality thresholds (`min_mapq × min_baseq`)

##### Setup
- Varied:
  - `min_mapq ∈ {0, 10, 20}`
  - `min_baseq ∈ {0, 10, 20}`
- Fixed: all other hard-scenario parameters

##### Figures to include
- **Figure 10.9.2.1 — Called effective phased recall across the `min_mapq × min_baseq` grid.**
- **Figure 10.9.2.2 — Called switch error across the `min_mapq × min_baseq` grid.**
- **Figure 10.9.2.3 — Called number of phase sets across the phasing quality-threshold grid (heatmap).**
- **Figure 10.9.2.4 — Oracle effective phased recall across the phasing quality-threshold grid (heatmap).**

**Figure 10.9.2.1 — Called effective phased recall across the `min_mapq × min_baseq` grid.**  
Called effective phased recall is consistently higher when `min_baseq = 0` or `10` than when `min_baseq = 20`, while `min_mapq` has little visible effect. This shows that overly strict base-quality filtering reduces usable phasing evidence.

**Figure 10.9.2.2 — Called switch error across the `min_mapq × min_baseq` grid.**  
Called switch error is lowest when `min_baseq = 0` or `10` and higher when `min_baseq = 20`, indicating that stricter base-quality filtering does not improve phasing correctness in this hard scenario and instead weakens phasing evidence.

**Figure 10.9.2.3 — Called number of phase sets across the `min_mapq × min_baseq` grid.**  
Phase fragmentation is substantially lower when `min_baseq = 0` or `10` than when `min_baseq = 20`, showing that strict base-quality filtering breaks block continuity by removing informative allele observations.

**Figure 10.9.2.4 — Oracle effective phased recall across the `min_mapq × min_baseq` grid.**  
Oracle phasing performance is markedly worse when `min_baseq = 20`, confirming that the main threshold effect arises within the phasing stage itself rather than through changes in the upstream callset.

##### Results summary
- For `min_baseq = 0 or 10`:
  - `called_effective_phased_recall`: 0.2964
  - `called_switch_error`: 0.0007
  - `called_num_phase_sets`: 6.6
  - `called_phasing_rate_shared_het`: 0.6027
  - `oracle_effective_phased_recall`: 0.9402–0.9415
  - `oracle_num_phase_sets`: 3.8
- For `min_baseq = 20`:
  - `called_effective_phased_recall`: 0.2919–0.2923
  - `called_switch_error`: 0.0056
  - `called_num_phase_sets`: 9.6
  - `called_phasing_rate_shared_het`: 0.5864–0.5872
  - `oracle_effective_phased_recall`: 0.8298–0.8319
  - `oracle_num_phase_sets`: 6.2–6.4
- `min_mapq` has negligible effect across the grid.
- `call_recall` remains constant at 0.5895 because variant-calling thresholds are fixed and only phasing-side filters are varied here.

##### Key observations
- O1 (Base-quality filtering is the main useful phasing-side knob): The grid separates into two clear regimes. Using `min_baseq = 0 or 10` improves called effective phased recall, lowers switch error, and reduces phase fragmentation relative to `min_baseq = 20`.

- O2 (Overly strict `min_baseq` harms both oracle and called phasing): Raising `min_baseq` to 20 reduces oracle effective phased recall from about `0.94` to about `0.83`, and increases oracle phase-set count from `3.8` to `6.2–6.4`. This shows that strict base-quality filtering discards too much informative phasing evidence.

- O3 (`min_mapq` contributes little under the current hard scenario): Changing `min_mapq` from 0 to 20 produces almost no measurable difference in called or oracle metrics, indicating that mapping-quality filtering is not a strong optimization lever here.

- O4 (The gain from threshold tuning is modest but real): Lowering `min_baseq` improves phasing completeness and continuity, but the absolute gain in called effective phased recall is still small because the dominant loss term remains missing shared heterozygous sites.

##### Takeaway
Phasing quality thresholds provide a limited but real optimization lever under the composite hard scenario. The most useful change is to avoid overly strict base-quality filtering: `min_baseq = 10` (or 0) clearly outperforms `min_baseq = 20`, while `min_mapq` has little practical impact.

#### 10.9.3 Sweep C: calling thresholds (`call_min_mapq × call_min_baseq`)

##### Setup
- Varied:
  - `call_min_mapq ∈ {0, 10, 20}`
  - `call_min_baseq ∈ {5, 10, 15, 20}`
- Fixed:
  - hard-scenario stressors unchanged
  - phasing parameters fixed to the best practical settings identified earlier:
    - `max_coverage = 10`
    - `min_mapq = 20`
    - `min_baseq = 10`

##### Figures to include
- **Figure 10.9.3.1 — Variant calling recall across the caller-threshold grid (heatmap).**
- **Figure 10.9.3.2 — Called effective phased recall across the `call_min_mapq × call_min_baseq` grid.**
- **Figure 10.9.3.3 — Called shared heterozygous recall across the `call_min_mapq × call_min_baseq` grid.**
- **Figure 10.9.3.4 — Call precision across the `call_min_mapq × call_min_baseq` grid.**

**Figure 10.9.3.1 — Calling recall across the `call_min_mapq × call_min_baseq` grid.**  
Calling recall depends strongly on `call_min_baseq` but only weakly on `call_min_mapq`. Moderate base-quality thresholds retain more true SNP evidence, whereas `call_min_baseq = 20` is overly strict and sharply reduces variant recovery.

**Figure 10.9.3.2 — Called effective phased recall across the `call_min_mapq × call_min_baseq` grid.**  
End-to-end phased recall is highest when `call_min_baseq = 5` or `10`, showing that caller-side base-quality filtering has a direct downstream impact on usable phasing performance.

**Figure 10.9.3.3 — Called shared heterozygous recall across the `call_min_mapq × call_min_baseq` grid.**  
Called shared heterozygous recall follows the same pattern as calling recall, indicating that the main benefit of moderate caller thresholds is better recovery of true heterozygous SNPs that can subsequently be phased.

**Figure 10.9.3.4 — Call precision across the `call_min_mapq × call_min_baseq` grid.**  
Call precision remains very high across the full threshold grid, showing that the gains from moderate caller base-quality thresholds are achieved without a substantial precision penalty.

##### Results summary
- For `call_min_baseq = 5 or 10`:
  - `call_precision`: ~0.9995
  - `call_recall`: ~0.596–0.597
  - `called_shared_het_recall`: ~0.517–0.518
  - `called_effective_phased_recall`: ~0.3155–0.3157
  - `called_phasing_rate_shared_het`: ~0.609–0.610
  - `called_switch_error`: 0.0000
  - `called_num_phase_sets`: 6.6
- For `call_min_baseq = 15`:
  - `call_recall`: ~0.589–0.590
  - `called_shared_het_recall`: ~0.509–0.510
  - `called_effective_phased_recall`: ~0.307–0.3074
- For `call_min_baseq = 20`:
  - `call_precision`: ~0.9991
  - `call_recall`: ~0.3645–0.3647
  - `called_shared_het_recall`: ~0.287
  - `called_effective_phased_recall`: ~0.1144–0.1146
  - `called_phasing_rate_shared_het`: ~0.399
  - `called_num_phase_sets`: 7.0
- `call_min_mapq` has negligible effect across the grid.

##### Key observations
- O1 (Caller base-quality filtering is a strong optimization lever): Lowering `call_min_baseq` from 20 to 10 or 5 substantially improves call recall, called shared-heterozygous recall, and called effective phased recall. This shows that overly strict caller base-quality filtering removes too much useful variant evidence under the composite hard scenario.

- O2 (Caller mapping-quality filtering has little practical effect): Changing `call_min_mapq` from 0 to 20 produces almost no measurable difference in calling or phasing metrics. Under the current hard scenario, base-quality filtering is therefore the dominant caller-side threshold knob.

- O3 (The gain arises through overlap recovery rather than phasing correctness): Oracle effective phased recall remains fixed (~0.9402), while called phase accuracy and switch error remain near-perfect in the best settings. This indicates that the improvement comes mainly from recovering more true shared heterozygous sites rather than from changing phasing correctness itself.

- O4 (Moderate caller thresholds outperform strict thresholds without sacrificing precision): The best-performing settings (`call_min_baseq = 5 or 10`) retain very high call precision (~0.9995) while increasing end-to-end phased recall. This makes moderate caller base-quality thresholds a practical optimization opportunity in this pipeline.

##### Takeaway
Calling thresholds provide a more effective optimization lever than the WhatsHap-only sweeps tested earlier. In particular, a moderate caller base-quality threshold (`call_min_baseq = 10`, or 5) improves variant recovery and end-to-end phased recall without materially harming precision, while `call_min_baseq = 20` is clearly too strict under the current hard scenario.

#### 10.9.4 Sweep D: runtime-focused lower `max_coverage`

##### Setup
- Varied: `max_coverage ∈ {4, 6, 8, 10, 15}`
- Fixed:
  - hard-scenario stressors unchanged
  - caller thresholds fixed to a strong practical setting:
    - `call_min_mapq = 20`
    - `call_min_baseq = 15`
  - phasing thresholds fixed to the best practical settings identified earlier:
    - `min_mapq = 20`
    - `min_baseq = 10`

##### Figures to include
- **Figure 10.9.4.1 — Called effective phased recall vs lower `max_coverage`.**
- **Figure 10.9.4.2 — Called number of phase sets vs lower `max_coverage`.**
- **Figure 10.9.4.3 — Oracle effective phased recall vs lower `max_coverage`.**
- **Figure 10.9.4.4 — Total pipeline runtime vs lower `max_coverage`.**

**Figure 10.9.4.1 — Called effective phased recall vs lower `max_coverage`.**  
Called effective phased recall is essentially unchanged for `max_coverage = 8, 10, 15`, indicating that retained phasing coverage above 8 provides no practical benefit in the current hard scenario.

**Figure 10.9.4.2 — Called number of phase sets vs lower `max_coverage`.**  
Phase fragmentation is stable from `max_coverage = 6` upward, but becomes slightly worse at `max_coverage = 4`, suggesting that very aggressive down-selection can begin to remove useful connectivity.

**Figure 10.9.4.3 — Oracle effective phased recall vs lower `max_coverage`.**  
Oracle phasing performance improves slightly from `max_coverage = 4` to `8`, then plateaus. This indicates that a small amount of retained coverage is necessary for full phasing continuity, but higher retained coverage beyond 8 is unnecessary.

**Figure 10.9.4.4 — Total pipeline runtime vs lower `max_coverage`.**  
Runtime is similar for `max_coverage = 6–10` but increases noticeably at `15`, showing that a lower retained-coverage cap can reduce compute cost without sacrificing phasing performance.

##### Results summary
- `call_recall`: 0.5895 at all `max_coverage` settings
- `call_precision`: 0.9994 at all `max_coverage` settings
- `called_effective_phased_recall`:
  - 0.2958 (`max_coverage = 4`)
  - 0.2960 (`6`)
  - 0.2964 (`8`)
  - 0.2964 (`10`)
  - 0.2964 (`15`)
- `called_switch_error`:
  - 0.0051 → 0.0038 → 0.0007 → 0.0007 → 0.0007
- `called_num_phase_sets`:
  - 6.8 → 6.6 → 6.6 → 6.6 → 6.6
- `called_phase_accuracy`:
  - 0.9604 → 0.9610 → 0.9625 → 0.9625 → 0.9625
- `oracle_effective_phased_recall`:
  - 0.9228 → 0.9348 → 0.9402 → 0.9402 → 0.9402
- `oracle_num_phase_sets`:
  - 5.0 → 3.8 → 3.8 → 3.8 → 3.8
- `time_total_sec`:
  - 3.6877 → 3.6720 → 3.6769 → 3.6774 → 4.2029

##### Key observations
- O1 (Performance plateaus by `max_coverage = 8`): Called effective phased recall, switch error, phase-set count, and phase accuracy are essentially identical for `max_coverage = 8, 10, 15`. This shows that retained phasing coverage above 8 does not provide additional practical benefit under the current hard scenario.

- O2 (Very low `max_coverage` begins to hurt modestly): At `max_coverage = 4`, oracle effective phased recall and oracle continuity are slightly worse, and called metrics are also marginally lower. This suggests that reducing retained coverage too aggressively can begin to remove useful bridging evidence.

- O3 (The best efficiency region is around `max_coverage = 8–10`): Runtime remains around 3.67 seconds for `max_coverage = 6–10`, but increases to about 4.20 seconds at `15` without any measurable gain in phasing performance.

- O4 (`max_coverage = 15` is unnecessarily expensive): Compared with `max_coverage = 8–10`, the setting `15` increases runtime by roughly 14% while leaving the main phasing metrics unchanged. This makes higher retained coverage unattractive on efficiency grounds.

##### Takeaway
A lower `max_coverage` setting provides a practical runtime optimization under the composite hard scenario. The best trade-off is around `max_coverage = 8` or `10`, where performance is indistinguishable from `15` but runtime is lower. In contrast, `max_coverage = 4` is slightly too aggressive and begins to degrade continuity-related metrics.

#### 10.9.5 Sweep E: fine phasing `min_baseq`

##### Setup
- Varied: `min_baseq ∈ {0, 5, 10, 15, 20}`
- Fixed:
  - hard-scenario stressors unchanged
  - caller thresholds fixed to a strong practical setting:
    - `call_min_mapq = 20`
    - `call_min_baseq = 15`
  - phasing parameters otherwise fixed:
    - `max_coverage = 10`
    - `min_mapq = 20`

##### Figures to include
- **Figure 10.9.5.1 — Called effective phased recall vs fine phasing `min_baseq`.**
- **Figure 10.9.5.2 — Called number of phase sets vs fine phasing `min_baseq`.**
- **Figure 10.9.5.3 — Oracle effective phased recall vs fine phasing `min_baseq`.**
- **Figure 10.9.5.4 — Called switch error vs fine phasing `min_baseq`.**

**Figure 10.9.5.1 — Called effective phased recall vs fine phasing `min_baseq`.**  
Called effective phased recall remains flat from `min_baseq = 0` to `15`, then decreases at `20`, indicating that only overly strict phasing base-quality filtering harms usable end-to-end phasing output.

**Figure 10.9.5.2 — Called number of phase sets vs fine phasing `min_baseq`.**  
Phase fragmentation is unchanged across `min_baseq = 0–15` but increases substantially at `20`, showing that strict phasing base-quality filtering breaks block continuity by discarding informative allele observations.

**Figure 10.9.5.3 — Oracle effective phased recall vs fine phasing `min_baseq`.**  
Oracle phasing performance is stable up to `min_baseq = 15` and then drops sharply at `20`, confirming that the harmful threshold effect arises within the phasing stage itself rather than through upstream callset changes.

**Figure 10.9.5.4 — Called switch error vs fine phasing `min_baseq`.**  
Called switch error remains very low and stable from `min_baseq = 0` to `15`, but increases at `20`, indicating that overly strict phasing base-quality filtering weakens the consistency of retained phasing evidence.

##### Results summary
- `call_recall`: 0.5895 at all `min_baseq` settings
- `call_precision`: 0.9994 at all `min_baseq` settings
- `called_effective_phased_recall`:
  - 0.2964 (`min_baseq = 0`)
  - 0.2964 (`5`)
  - 0.2964 (`10`)
  - 0.2964 (`15`)
  - 0.2917 (`20`)
- `called_switch_error`:
  - 0.0007 → 0.0007 → 0.0007 → 0.0007 → 0.0056
- `called_num_phase_sets`:
  - 6.6 → 6.6 → 6.6 → 6.6 → 9.6
- `called_phase_accuracy`:
  - 0.9625 → 0.9625 → 0.9625 → 0.9625 → 0.9748
- `oracle_effective_phased_recall`:
  - 0.9402 → 0.9402 → 0.9402 → 0.9375 → 0.8256
- `oracle_num_phase_sets`:
  - 3.8 → 3.8 → 3.8 → 3.8 → 6.4
- `time_total_sec`:
  - 3.6758 → 3.6505 → 3.6501 → 3.6172 → 3.6158

##### Key observations
- O1 (Phasing performance is flat from `min_baseq = 0` to `15`): Called effective phased recall, switch error, phase-set count, and phase accuracy are effectively identical across `min_baseq ∈ {0, 5, 10, 15}`. This indicates that moderate variation in phasing base-quality filtering has no measurable effect within this range.

- O2 (The harmful threshold is `min_baseq = 20`): Raising `min_baseq` to 20 reduces called effective phased recall, increases switch error, and increases phase fragmentation. Oracle effective phased recall also drops sharply (`0.9402 → 0.8256`), confirming that the loss is caused within the phasing stage itself.

- O3 (There is no practical benefit to using values below 10): Since `0`, `5`, and `10` perform identically, more permissive phasing base-quality thresholds do not yield additional measurable gains under the current hard scenario.

- O4 (`min_baseq = 10` is the best practical recommendation): Because `10` lies on the same performance plateau as `0`, `5`, and `15`, but is easier to justify as a moderate threshold, it provides the cleanest practical phasing setting for the remainder of the study.

##### Takeaway
The fine phasing `min_baseq` sweep confirms that performance is stable across `min_baseq = 0–15` and degrades only when the threshold becomes too strict (`20`). A moderate phasing base-quality threshold such as `min_baseq = 10` is therefore the most practical recommendation: it retains full performance while avoiding unnecessary permissiveness and clearly outperforming the overly strict setting.

#### 10.9.6 Default vs optimized comparison

##### Setup
A final confirmation run compared the default hard-scenario configuration against the best practical tuned configuration identified from the preceding optimization sweeps.

- **Default configuration**
  - `call_min_mapq = 20`
  - `call_min_baseq = 15`
  - `max_coverage = 15`
  - `min_mapq = 20`
  - `min_baseq = 20`

- **Optimized configuration**
  - `call_min_mapq = 20`
  - `call_min_baseq = 10`
  - `max_coverage = 8`
  - `min_mapq = 20`
  - `min_baseq = 10`

All other hard-scenario stressors were kept fixed.

##### Figures to include
- **Figure 10.9.6.1 — Called effective phased recall: default vs optimized.**
- **Figure 10.9.6.2 — Called shared heterozygous recall: default vs optimized.**
- **Figure 10.9.6.3 — Called number of phase sets: default vs optimized.**
- **Figure 10.9.6.4 — Total pipeline runtime: default vs optimized.**

**Figure 10.9.6.1 — Called effective phased recall: default vs optimized.**  
The optimized configuration achieves higher end-to-end phased recall than the default hard-scenario setting, confirming that coordinated parameter tuning yields a measurable improvement under compounded realistic stress.

**Figure 10.9.6.2 — Called shared heterozygous recall: default vs optimized.**  
The optimized configuration increases called shared heterozygous recall, indicating that part of the improvement comes from recovering more true heterozygous sites that can subsequently be phased.

**Figure 10.9.6.3 — Called number of phase sets: default vs optimized.**  
The optimized configuration substantially reduces the number of phase sets, showing that it preserves stronger phase-block continuity than the default setting under the hard scenario.

**Figure 10.9.6.4 — Total pipeline runtime: default vs optimized.**  
The optimized configuration is slightly faster than the default while also improving phasing performance, showing that the observed gains are not achieved at the cost of higher compute.

##### Results summary
Means over 5 seeds:

- `time_total_sec`: 4.0614 → 3.9783
- `call_precision`: 0.9994 → 0.9995
- `call_recall`: 0.5895 → 0.5970
- `oracle_effective_phased_recall`: 0.8298 → 0.9402
- `oracle_switch_error`: 0.0083 → 0.0031
- `oracle_num_phase_sets`: 6.2 → 3.8
- `called_effective_phased_recall`: 0.2919 → 0.3041
- `called_shared_het_recall`: 0.5098 → 0.5181
- `called_phasing_rate_shared_het`: 0.5864 → 0.6090
- `called_phase_accuracy`: 0.9748 → 0.9622
- `called_switch_error`: 0.0056 → 0.0007
- `called_num_phase_sets`: 9.6 → 6.6

##### Key observations
- O1 (The optimized configuration improves end-to-end phased recall): Called effective phased recall increases from `0.2919` to `0.3041`, confirming that the tuned configuration provides a real, if moderate, end-to-end gain under the composite hard scenario.

- O2 (The improvement comes mainly from better overlap and phasing completeness): The optimized configuration increases `call_recall` (`0.5895 → 0.5970`), `called_shared_het_recall` (`0.5098 → 0.5181`), and `called_phasing_rate_shared_het` (`0.5864 → 0.6090`). This shows that the gain arises primarily from recovering and retaining more useful phasing evidence rather than from a large change in phase accuracy itself.

- O3 (Fragmentation is substantially reduced): Called phase-set count decreases from `9.6` to `6.6`, while oracle phase-set count decreases from `6.2` to `3.8`. This indicates that the optimized configuration preserves much stronger block continuity under the hard scenario.

- O4 (The optimized configuration is also slightly faster): Total runtime decreases from `4.0614 s` to `3.9783 s`, showing that the optimized setting improves performance without increasing compute cost.

##### Takeaway
The final confirmation run shows that the optimized configuration is practically preferable to the default hard-scenario setting. It improves end-to-end phased recall, increases shared-site overlap and phasing completeness, reduces fragmentation, lowers switch error, and slightly reduces runtime. Although the absolute gain in called effective phased recall is moderate, this experiment confirms that coordinated tuning of caller-side and phasing-side parameters can produce a measurable improvement under compounded realistic stress.

#### 10.9.7 Sweep F: local joint search around recommended thresholds

##### Setup
A focused local search was performed around the current best practical threshold settings to test whether the recommended caller/phasing base-quality pair is locally optimal or whether a nearby combination performs better.

- Varied:
  - `call_min_baseq ∈ {5, 10, 15}`
  - `min_baseq ∈ {5, 10, 15}`
- Fixed:
  - `call_min_mapq = 20`
  - `min_mapq = 20`
  - `max_coverage = 8`
  - hard-scenario stressors unchanged

##### Figures to include
- **Figure 10.9.7.1 — Called effective phased recall across the local joint base-quality grid (heatmap).**
- **Figure 10.9.7.2 — Called shared heterozygous recall across the local caller/phasing base-quality grid.**
- **Figure 10.9.7.3 — Calling recall across the local caller/phasing base-quality grid.**
- **Figure 10.9.7.4 — Total pipeline runtime across the local caller/phasing base-quality grid.**

**Figure 10.9.7.1 — Called effective phased recall across the local caller/phasing base-quality grid.**  
Called effective phased recall is highest whenever `call_min_baseq = 5` or `10`, while varying phasing `min_baseq` between 5 and 15 has no measurable effect. This shows that local sensitivity is driven mainly by the caller-side threshold.

**Figure 10.9.7.2 — Called shared heterozygous recall across the local caller/phasing base-quality grid.**  
Called shared heterozygous recall follows the same pattern as called effective phased recall, indicating that the local optimization gain arises from recovering more true heterozygous sites rather than from changing phasing correctness.

**Figure 10.9.7.3 — Calling recall across the local caller/phasing base-quality grid.**  
Calling recall improves when `call_min_baseq` is reduced from 15 to 10 or 5, while phasing `min_baseq` has no effect as expected. This confirms that the decisive local optimization lever is the caller base-quality threshold.

**Figure 10.9.7.4 — Total pipeline runtime across the local caller/phasing base-quality grid.**  
Runtime varies only slightly across the local search grid, indicating that the local threshold refinement primarily affects accuracy-related metrics rather than compute cost.

##### Results summary
Grouped means over 5 seeds:

- For `call_min_baseq = 5 or 10`:
  - `call_recall`: 0.5970
  - `called_shared_het_recall`: 0.5181
  - `called_effective_phased_recall`: 0.3041
  - `called_phasing_rate_shared_het`: 0.6090
  - `called_switch_error`: 0.0007
  - `called_num_phase_sets`: 6.6
- For `call_min_baseq = 15`:
  - `call_recall`: 0.5895
  - `called_shared_het_recall`: 0.5098
  - `called_effective_phased_recall`: 0.2964
  - `called_phasing_rate_shared_het`: 0.6027
  - `called_switch_error`: 0.0007
  - `called_num_phase_sets`: 6.6
- Across `min_baseq ∈ {5, 10, 15}`, the main called metrics are effectively unchanged for a fixed `call_min_baseq`.

##### Key observations
- O1 (The recommended setting is locally robust): No combination in the 3×3 local grid outperforms the existing recommended region. The best results are obtained whenever `call_min_baseq = 5 or 10`, regardless of whether phasing `min_baseq` is 5, 10, or 15.

- O2 (Caller base-quality remains the decisive local knob): Lowering `call_min_baseq` from 15 to 10 or 5 improves `call_recall`, `called_shared_het_recall`, and `called_effective_phased_recall`. This confirms that the main optimization gain still comes from recovering more useful variant evidence upstream.

- O3 (Phasing base-quality is insensitive within the moderate range): For fixed caller thresholds, varying phasing `min_baseq` between 5 and 15 has no measurable effect on the main called metrics. This indicates that, once overly strict filtering has been avoided, phasing base-quality is not a sensitive local optimization lever.

- O4 (`call_min_baseq = 10`, `min_baseq = 10` remains the best practical recommendation): Although `call_min_baseq = 5` performs equivalently to 10, the value 10 is easier to justify as a moderate threshold. Likewise, phasing `min_baseq = 10` lies on the same performance plateau as 5 and 15, making it the cleanest practical choice.

##### Takeaway
The local search confirms that the recommended optimized threshold pair is robust. The main local sensitivity lies in the caller base-quality threshold, where 10 (or 5) is clearly better than 15, while phasing base-quality is effectively flat within the moderate range 5–15. This strengthens confidence that the chosen optimized configuration is not a fragile one-off result.

#### 10.9.8 Sweep G: accuracy/runtime frontier comparison

##### Setup
Representative configurations were compared under the same hard scenario to summarize the practical trade-off between end-to-end phased recall, block continuity, and runtime.

Configurations compared:
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
- **`balanced`**
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

All other hard-scenario stressors were kept fixed.

##### Figures to include
- **Figure 10.9.8.1 — Called effective phased recall across representative configurations.**
- **Figure 10.9.8.2 — Called shared heterozygous recall across representative configurations.**
- **Figure 10.9.8.3 — Called number of phase sets across representative configurations.**
- **Figure 10.9.8.4 — Total pipeline runtime across representative configurations.**

\begin{center}
\includegraphics[width=0.82\linewidth]{figures/10_experiments/fig_10_9_8_1_called_effective_phased_recall.png}
\end{center}

**Figure 10.9.8.1 — Called effective phased recall across representative configurations.**
The balanced configuration achieves the highest end-to-end phased recall, while the runtime-biased configuration performs almost identically. This shows that combining caller-side and phasing-side tuning yields the strongest practical performance.

**Figure 10.9.8.2 — Called shared heterozygous recall across representative configurations.**
Called shared heterozygous recall is highest for the configurations that use the improved caller threshold, confirming that part of the overall gain comes from recovering more true heterozygous sites before phasing.

**Figure 10.9.8.3 — Called number of phase sets across representative configurations.**  
Phase fragmentation is much lower in the phasing-tuned, balanced, and runtime configurations than in the default or caller-only configurations. This shows that phasing-side tuning is the main driver of improved block continuity.

**Figure 10.9.8.4 — Total pipeline runtime across representative configurations.**  
Runtime differences between the representative configurations are modest, but the balanced and runtime-biased settings both remain competitive while outperforming the default configuration on the main phasing metrics.

##### Results summary
Means over 5 seeds:

- **`default`**
  - `time_total_sec`: 3.9338
  - `call_precision`: 0.9994
  - `call_recall`: 0.5895
  - `oracle_effective_phased_recall`: 0.8298
  - `oracle_num_phase_sets`: 6.2
  - `called_effective_phased_recall`: 0.2919
  - `called_shared_het_recall`: 0.5098
  - `called_phasing_rate_shared_het`: 0.5864
  - `called_switch_error`: 0.0056
  - `called_num_phase_sets`: 9.6

- **`caller_only`**
  - `time_total_sec`: 4.0540
  - `call_precision`: 0.9995
  - `call_recall`: 0.5970
  - `oracle_effective_phased_recall`: 0.8298
  - `oracle_num_phase_sets`: 6.2
  - `called_effective_phased_recall`: 0.2994
  - `called_shared_het_recall`: 0.5181
  - `called_phasing_rate_shared_het`: 0.5934
  - `called_switch_error`: 0.0077
  - `called_num_phase_sets`: 9.0

- **`phasing_only`**
  - `time_total_sec`: 3.7970
  - `call_precision`: 0.9994
  - `call_recall`: 0.5895
  - `oracle_effective_phased_recall`: 0.9402
  - `oracle_num_phase_sets`: 3.8
  - `called_effective_phased_recall`: 0.2964
  - `called_shared_het_recall`: 0.5098
  - `called_phasing_rate_shared_het`: 0.6027
  - `called_switch_error`: 0.0007
  - `called_num_phase_sets`: 6.6

- **`balanced`**
  - `time_total_sec`: 3.8682
  - `call_precision`: 0.9995
  - `call_recall`: 0.5970
  - `oracle_effective_phased_recall`: 0.9402
  - `oracle_num_phase_sets`: 3.8
  - `called_effective_phased_recall`: 0.3041
  - `called_shared_het_recall`: 0.5181
  - `called_phasing_rate_shared_het`: 0.6090
  - `called_switch_error`: 0.0007
  - `called_num_phase_sets`: 6.6

- **`runtime`**
  - `time_total_sec`: 3.8504
  - `call_precision`: 0.9995
  - `call_recall`: 0.5970
  - `oracle_effective_phased_recall`: 0.9348
  - `oracle_num_phase_sets`: 3.8
  - `called_effective_phased_recall`: 0.3039
  - `called_shared_het_recall`: 0.5181
  - `called_phasing_rate_shared_het`: 0.6090
  - `called_switch_error`: 0.0022
  - `called_num_phase_sets`: 6.6

##### Key observations
- O1 (The balanced configuration provides the best overall trade-off): Among the representative configurations, `balanced` achieves the highest called effective phased recall (`0.3041`) while also maintaining low fragmentation (`6.6` phase sets), very low switch error (`0.0007`), and lower runtime than the default configuration.

- O2 (Default is dominated by the tuned configurations): The default configuration is worse than `balanced` and `runtime` on called effective phased recall, called shared heterozygous recall, switch error, and fragmentation, while also not being the fastest option. This indicates that the untuned hard-scenario settings are not practically competitive.

- O3 (Caller-only and phasing-only tuning recover complementary parts of the gain): `caller_only` improves overlap-related metrics (`call_recall`, `called_shared_het_recall`) but leaves fragmentation high, whereas `phasing_only` strongly improves continuity and oracle phasing metrics but does not recover additional called-site overlap. The balanced configuration combines both benefits.

- O4 (The runtime-biased configuration is a viable alternative): `runtime` performs almost identically to `balanced` on the main called metrics (`0.3039` vs `0.3041` called effective phased recall) while using a slightly lower `max_coverage`. This makes it a reasonable alternative when efficiency is prioritized, although `balanced` remains the cleaner general recommendation.

##### Takeaway
The frontier comparison shows that the best practical performance comes from combining moderate caller-side and phasing-side tuning rather than optimizing either stage in isolation. The balanced configuration is the preferred general setting because it delivers the highest end-to-end phased recall together with strong continuity and low switch error. A slightly more runtime-biased configuration performs almost identically and may be acceptable when efficiency is prioritized, but the default configuration is clearly dominated.

#### 10.9.9 Optimization robustness across scenarios

##### Setup
To test whether the recommended optimized settings generalize beyond the hard tuning scenario, three representative configurations were compared across four scenarios:

Configurations compared:
- **`default`**
- **`balanced`**
- **`runtime`**

Scenarios compared:
- **baseline**
- **dropout**
- **interaction**
- **hard**

The aim was to determine whether the recommended tuning is broadly useful or only beneficial in the hard scenario used for parameter selection.

##### Figures to include
- **Figure 10.9.9.1 — Called effective phased recall across scenarios and representative configurations (heatmap).**
- **Figure 10.9.9.2 — Called shared heterozygous recall across scenarios and representative configurations.**
- **Figure 10.9.9.3 — Called number of phase sets across scenarios and representative configurations.**
- **Figure 10.9.9.4 — Total pipeline runtime across scenarios and representative configurations.**

\begin{center}
\includegraphics[width=0.88\linewidth]{figures/10_experiments/fig_10_9_9_1_called_effective_phased_recall_heatmap.png}
\end{center}

**Figure 10.9.9.1 — Called effective phased recall across scenarios and representative configurations.**
Across all tested scenarios, the balanced and runtime configurations achieve higher end-to-end phased recall than the default configuration, showing that the selected tuning generalizes beyond the hard optimization scenario.

**Figure 10.9.9.2 — Called shared heterozygous recall across scenarios and representative configurations.**  
The tuned configurations consistently improve called shared heterozygous recall, indicating that part of the transferable gain comes from recovering more true heterozygous sites across multiple scenarios.

**Figure 10.9.9.3 — Called number of phase sets across scenarios and representative configurations.**  
The balanced and runtime configurations reduce phase fragmentation relative to the default configuration in all tested scenarios, showing that the tuning also improves block continuity rather than only site recovery.

**Figure 10.9.9.4 — Total pipeline runtime across scenarios and representative configurations.**  
Runtime remains competitive for the tuned configurations across scenarios, and the runtime-biased variant offers nearly identical accuracy with similar or slightly improved efficiency.

##### Results summary
Means over 5 seeds:

- **Baseline**
  - `default`: `called_effective_phased_recall = 0.4364`, `called_shared_het_recall = 0.5902`, `called_num_phase_sets = 3.4`, `time_total_sec = 2.9850`
  - `balanced`: `0.5448`, `0.6097`, `1.4`, `3.0223`
  - `runtime`: `0.5448`, `0.6097`, `1.4`, `2.9596`

- **Dropout**
  - `default`: `0.2694`, `0.5292`, `6.2`, `2.9768`
  - `balanced`: `0.2908`, `0.5420`, `4.2`, `2.9362`
  - `runtime`: `0.2908`, `0.5420`, `4.2`, `2.9391`

- **Interaction**
  - `default`: `0.1972`, `0.5014`, `6.2`, `3.0043`
  - `balanced`: `0.2118`, `0.5094`, `5.4`, `2.8789`
  - `runtime`: `0.2118`, `0.5094`, `5.4`, `2.9750`

- **Hard**
  - `default`: `0.2919`, `0.5098`, `9.6`, `4.3644`
  - `balanced`: `0.3041`, `0.5181`, `6.6`, `4.2497`
  - `runtime`: `0.3039`, `0.5181`, `6.6`, `4.3050`

##### Key observations
- O1 (The tuned configurations generalize across scenarios): In all four scenarios, both `balanced` and `runtime` outperform `default` on called effective phased recall and called shared heterozygous recall. This shows that the tuned settings are not overfit to the hard scenario alone.

- O2 (The balanced configuration is the safest general recommendation): `balanced` is never worse than `default` and is either best or tied for best on the main called metrics across all scenarios. This supports its use as the default recommended configuration.

- O3 (The runtime-biased configuration is highly competitive): `runtime` is effectively identical to `balanced` on the main called metrics in this sweep and remains competitive on runtime. This makes it a reasonable alternative when efficiency is prioritized.

- O4 (The gains transfer to both easy and hard cases): The tuned configurations improve performance not only in the hard scenario, but also in the baseline, dropout, and interaction scenarios. This indicates that the optimization is robust across qualitatively different data conditions.

##### Takeaway
The recommended tuned settings are robust across the representative scenarios tested in this study. The balanced configuration is the cleanest general recommendation because it consistently improves end-to-end phased recall and reduces fragmentation relative to the default. A runtime-biased alternative performs almost identically and may be preferred when efficiency is prioritized.

#### 10.9.10 Sweep H: caller `call_min_mapq` rule-out

##### Setup
A final small rule-out sweep was conducted to test whether caller mapping-quality filtering provides meaningful optimization headroom once the other recommended settings are fixed.

- Varied:
  - `call_min_mapq ∈ {0, 10, 20}`
- Fixed:
  - `call_min_baseq = 10`
  - `max_coverage = 8`
  - `min_mapq = 20`
  - `min_baseq = 10`
  - hard-scenario stressors unchanged

##### Figures to include
- **Figure 10.9.10.1 — Called effective phased recall vs `call_min_mapq`.**
- **Figure 10.9.10.2 — Calling recall vs `call_min_mapq`.**
- **Figure 10.9.10.3 — Called shared heterozygous recall vs `call_min_mapq`.**
- **Figure 10.9.10.4 — Call precision vs `call_min_mapq`.**

**Figure 10.9.10.1 — Called effective phased recall vs `call_min_mapq`.**  
Called effective phased recall changes only slightly across the caller mapping-quality sweep, indicating that `call_min_mapq` is not a strong optimization lever under the current hard scenario.

**Figure 10.9.10.2 — Calling recall vs `call_min_mapq`.**  
Calling recall remains essentially unchanged across the caller mapping-quality sweep, showing that mapQ filtering has little influence on variant recovery in this setup.

**Figure 10.9.10.3 — Called shared heterozygous recall vs `call_min_mapq`.**  
Called shared heterozygous recall is nearly flat across the caller mapping-quality sweep, confirming that caller mapQ does not materially alter the overlap of useful true heterozygous sites.

**Figure 10.9.10.4 — Call precision vs `call_min_mapq`.**  
Call precision remains extremely high across all caller mapping-quality settings, indicating that stricter mapQ filtering does not provide a meaningful precision advantage in this pipeline.

##### Results summary
Means over 5 seeds:

- `call_min_mapq = 0`
  - `call_precision = 0.9995`
  - `call_recall = 0.5962`
  - `called_shared_het_recall = 0.5171`
  - `called_effective_phased_recall = 0.3155`
  - `called_phasing_rate_shared_het = 0.6102`
  - `called_phase_accuracy = 1.0000`
  - `called_switch_error = 0.0000`
  - `time_total_sec = 4.1980`

- `call_min_mapq = 10`
  - `call_precision = 0.9995`
  - `call_recall = 0.5970`
  - `called_shared_het_recall = 0.5181`
  - `called_effective_phased_recall = 0.3157`
  - `called_phasing_rate_shared_het = 0.6094`
  - `called_phase_accuracy = 1.0000`
  - `called_switch_error = 0.0000`
  - `time_total_sec = 4.1432`

- `call_min_mapq = 20`
  - `call_precision = 0.9995`
  - `call_recall = 0.5970`
  - `called_shared_het_recall = 0.5181`
  - `called_effective_phased_recall = 0.3041`
  - `called_phasing_rate_shared_het = 0.6090`
  - `called_phase_accuracy = 0.9622`
  - `called_switch_error = 0.0007`
  - `time_total_sec = 4.1246`

##### Key observations
- O1 (Caller mapping-quality filtering has little effect on overlap recovery): Calling recall and called shared heterozygous recall remain essentially unchanged across `call_min_mapq ∈ {0, 10, 20}`. This indicates that caller mapping-quality filtering does not materially alter the amount of useful variant evidence reaching the phasing stage.

- O2 (The optimization signal is much weaker than for caller base-quality): In contrast to the earlier caller base-quality sweep, changing `call_min_mapq` produces only very small differences in the main called metrics. This shows that caller mapping-quality is not a strong optimization lever in the current pipeline.

- O3 (No compelling benefit is obtained from stricter caller mapQ filtering): Although `call_min_mapq = 20` shows a slight reduction in called effective phased recall, the effect is small and not accompanied by meaningful changes in call precision or shared-site recall. This suggests that aggressive caller mapQ filtering provides no practical advantage here.

##### Takeaway
Caller mapping-quality filtering does not provide meaningful optimization headroom under the current hard scenario. Unlike caller base-quality, which showed a clear and useful tuning signal, `call_min_mapq` has only weak effects on the main end-to-end phasing metrics. This supports treating caller mapQ as a low-priority or effectively negligible optimization knob in this pipeline.

#### 10.9.11 DP-scaling under hard conditions

##### Aim
Assess how the downstream phasing problem scales with increasing SNP density under the composite hard scenario, and determine whether the tuned configuration changes the runtime growth of the solve stage.

##### Setup
A small scaling study was conducted under the fixed hard scenario, varying only the number of SNPs:

- `num_snps ∈ {400, 800, 1200, 1600}`

Two configurations were compared:
- **`default`**
  - `call_min_baseq = 15`
  - `max_coverage = 15`
  - `min_baseq = 20`
- **`balanced`**
  - `call_min_baseq = 10`
  - `max_coverage = 8`
  - `min_baseq = 10`

All other hard-scenario stressors were kept fixed.

##### Figures to include
- **Figure 10.9.11.1 — Called solve time vs SNP count under default and optimized configurations.**
- **Figure 10.9.11.2 — Oracle solve time vs SNP count under default and optimized configurations.**
- **Figure 10.9.11.3 — Called effective phased recall vs SNP count under default and optimized configurations.**
- **Figure 10.9.11.4 — Called number of phase sets vs SNP count under default and optimized configurations.**

\begin{center}
\includegraphics[width=0.82\linewidth]{figures/10_experiments/fig_10_9_11_1_called_time_solve_sec.png}
\end{center}

**Figure 10.9.11.1 — Called solve time vs SNP count under default and optimized configurations.**
Called solve time increases with SNP density under both configurations, but grows much more steeply under the default setting. This shows that the tuned configuration makes the residual phasing problem easier to solve as variant density increases.

**Figure 10.9.11.2 — Oracle solve time vs SNP count under default and optimized configurations.**  
Oracle solve time also increases with SNP density, confirming that the solver becomes more expensive as the phasing-only problem grows. The optimized configuration remains substantially faster, indicating that its benefit is not solely due to changes in the called callset.

**Figure 10.9.11.3 — Called effective phased recall vs SNP count under default and optimized configurations.**  
Across the scaling range, the optimized configuration maintains slightly higher end-to-end phased recall than the default configuration, showing that its runtime benefit is not obtained by sacrificing phasing quality.

**Figure 10.9.11.4 — Called number of phase sets vs SNP count under default and optimized configurations.**  
The optimized configuration generally produces fewer phase sets than the default configuration as SNP density increases, indicating better phase-block continuity and a less fragmented residual phasing problem.

##### Results summary
Means over 3 seeds:

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

##### Key observations
- O1 (Solve time increases with SNP density): In both oracle and called regimes, solve time increases as `num_snps` increases from 400 to 1600, showing that the downstream phasing problem becomes more expensive as variant density rises.

- O2 (The optimized configuration scales much better than default): The growth in solve time is much steeper under the default configuration than under the optimized configuration. At `num_snps = 1600`, default solve time is substantially higher than optimized in both oracle and called regimes, indicating that the tuned settings make the residual phasing problem easier to solve.

- O3 (optimized remains better on phasing continuity as the problem grows): Across the scaling range, the optimized configuration generally maintains fewer phase sets and slightly better phased recall than the default configuration. This suggests that the tuning improves not only final phasing output but also the structure of the instance presented to the solver.

- O4 (Solve time still does not dominate total runtime end-to-end): Although solve time grows with SNP density, total phasing runtime remains dominated by preprocessing/readset construction in the current research pipeline. Therefore, the main value of this result is algorithmic interpretation rather than end-to-end runtime engineering.

##### Takeaway
The DP-scaling study shows that increasing SNP density makes the phasing solve stage progressively more expensive, but also that the optimized tuned configuration scales substantially better than the default configuration. This suggests that practical tuning can improve not only phased output quality but also the difficulty of the residual wMEC-style problem faced by the solver.

#### Overall takeaway
Under the current hard scenario, optimization headroom exists but is uneven across knobs. WhatsHap-side tuning provides only limited gains: increasing `max_coverage` does not help, and the main useful phasing-side adjustment is to avoid overly strict `min_baseq` filtering. In contrast, caller-side tuning is more impactful: lowering `call_min_baseq` from 20 to a moderate value such as 10 improves variant recovery, shared heterozygous overlap, and end-to-end phased recall without materially reducing precision. Taken together, these results indicate that practical optimization in this pipeline is not purely a WhatsHap parameter problem; meaningful gains are more likely to come from a combination of permissive phasing evidence retention and improved upstream variant recovery.

---

### 10.10 Summary of results (cross-experiment synthesis)

#### 10.10.1 Attribution summary (oracle vs called)

Taken together, the oracle-vs-called comparison provides a consistent attribution picture across the stress studies.

- **Duplications:** primarily a **called-regime phasing/evidence-quality stressor**, rather than a strong caller-limited bottleneck. Calling recall remains nearly unchanged across duplication levels, and oracle effective phased recall remains high and stable. The main degradation appears only at the strongest duplication level, where called phasing completeness and phase accuracy worsen without a corresponding drop in shared-site recall. This indicates that duplicated regions mainly weaken the consistency of phasing evidence on surviving called sites.

- **Coverage dropout:** primarily a **fragmentation-dominated stressor**, becoming **mixed** at high severity. Oracle phase-set counts increase strongly with dropout fraction and oracle effective phased recall declines, showing that dropout directly breaks read connectivity even when correct sites are supplied. In the called regime, performance degrades more strongly still, and at severe dropout both shared-site recall and phasing completeness collapse. This makes dropout the clearest example of a connectivity bottleneck that later becomes both a calling and phasing bottleneck.

- **Correlated error bursts:** a **weak and somewhat noisy called-regime correctness stressor** under the current parameterization. Oracle metrics remain stable, calling recall is broadly unchanged, and the main visible effect is a dip in called-regime phase accuracy and switch stability at intermediate burst probability. Because the response is non-monotonic, this stressor should be interpreted cautiously and does not currently provide a strong optimization signal.

- **Read length model:** a **mixed continuity/overlap trade-off** rather than a clear stressor or optimization knob. The lognormal model improves calling recall, shared-site recall, and phase-block continuity, but these gains are largely offset by slightly worse called-regime phasing correctness. This shows that better connectivity alone does not guarantee a strong end-to-end gain.

- **Truth indels with SNP-only policy:** primarily a **methodology-validity result**, not a performance stressor. Oracle SNP-phasing metrics remain stable and called SNP-level metrics stay within a comparable range across indel settings, indicating that `phase_snps_only` and `eval_snps_only` successfully preserve interpretable SNP phasing evaluation even when truth indels are present.

- **Duplication × dropout interaction:** a **mixed compounded weakness**, strongest in the called regime. Dropout remains the main source of fragmentation, but duplication further worsens both called-site overlap and phasing completeness on surviving sites. This interaction shows that realistic end-to-end failure modes arise most clearly when ambiguous mapping and low-coverage gaps are combined.

Overall, the attribution pattern is consistent: whenever oracle performance remains strong but called performance degrades, the dominant bottleneck lies in the overlap and quality of called heterozygous sites. When oracle performance also degrades, the stressor is directly harming phase-block connectivity or phasing completeness.

#### 10.10.2 Most impactful stressors

Across the single-knob and interaction studies, the realism knobs differ substantially in practical impact.

- **Most impactful overall:** **coverage dropout**. It produces the clearest and strongest reduction in called effective phased recall, strongly increases phase fragmentation, and eventually also reduces calling recall. It is the most important isolated stressor in this study.

- **Most impactful compounded weakness:** **duplication × dropout interaction**. While duplication alone is relatively mild, its combination with dropout creates a substantially worse end-to-end condition than either stressor alone, especially in the called regime. This reveals a realistic mixed failure mode involving both reduced callable overlap and weaker phasing connectivity.

- **Moderate impact:** **read length model** and **strong duplication**. The read length model changes continuity and overlap in opposite directions, while strong duplication mainly affects called-regime phasing quality at the highest severity.

- **Weak / noisy impact:** **correlated error bursts**. Under the current parameterization, bursts do not produce a strong monotonic degradation trend and therefore are less informative as a practical optimization target.

- **Methodological rather than stress-inducing:** **truth indels under SNP-only policy**. Their main value is to validate that the evaluation remains interpretable under indel-containing truth sets.

If stressors are ranked by effect on the main end-to-end metric (`called_effective_phased_recall`), the most important findings are:  
1. **dropout**,  
2. **duplication × dropout interaction**,  
3. **strong duplication / read-length trade-off effects**,  
4. **bursts (weak/noisy)**,  
with the indel study serving primarily as an evaluation-robustness check rather than a degradation study.

#### 10.10.3 Optimization summary

Optimization sweeps under the composite hard scenario showed that optimization headroom exists, but is uneven across knobs.

A small runtime breakdown of the called phasing stage was also examined for representative configurations under the hard scenario:

| Configuration | Build ReadSet (s) | Read Selection (s) | Solve (s) | Called Phasing Total (s) |
|---|---:|---:|---:|---:|
| Optimized | 0.3679 | 0.0061 | 0.0056 | 0.3854 |
| Default | 0.3734 | 0.0023 | 0.0610 | 0.4433 |
| Speed-oriented | 0.3685 | 0.0050 | 0.0030 | 0.3827 |

This runtime breakdown shows that, in the current research pipeline, the largest share of called phasing runtime is spent in readset construction, while read selection and the downstream solve stage contribute smaller but configuration-dependent fractions. In particular, the default configuration has a noticeably larger solve-stage cost than the tuned configurations, whereas the optimized and speed-oriented settings are dominated even more strongly by readset construction. This result should be interpreted cautiously, because the present study uses a custom adaptor to support oracle/called benchmarking and unified reporting rather than the standard production-facing WhatsHap front-end. Nevertheless, it suggests that evidence extraction can be a significant practical cost in WhatsHap-based phasing workflows, and that the tuned configurations improve runtime mainly through a better overall phasing setup rather than through a large reduction in the core solve stage alone.

The main optimization findings are:

- **Useful caller-side knob:** `call_min_baseq`. Lowering the caller base-quality threshold from 20 to a moderate value such as 10 substantially improves calling recall, shared heterozygous overlap, and end-to-end phased recall without materially harming precision. This is the strongest optimization signal found in the study.

- **Useful phasing-side knob:** `min_baseq`, but mainly as a **rule-out of overly strict filtering**. Phasing performance is stable for `min_baseq = 0–15`, but degrades at `20`, indicating that strict phasing base-quality filtering discards too much informative evidence. A moderate setting such as `min_baseq = 10` is the cleanest practical recommendation.

- **Useful runtime knob:** `max_coverage`. Increasing `max_coverage` beyond 8–10 does not improve phasing performance, while lower values such as 8 preserve full performance with lower runtime. This makes `max_coverage = 8` a practical efficiency setting.

- **Weak / negligible knobs:** caller `call_min_mapq`, phasing `min_mapq`, and high `max_coverage` values above the plateau. These do not provide meaningful practical gains under the tested hard scenario.

- **Robustness of tuning:** the tuned configurations generalize across baseline, dropout, interaction, and hard scenarios. In all tested cases, the balanced tuned configuration improves called effective phased recall and reduces fragmentation relative to default, while a runtime-biased configuration performs almost identically and remains a viable alternative when efficiency is prioritized.

- **Algorithm-facing scaling result:** increasing SNP density makes the solve stage more expensive, but the balanced tuned configuration scales much better than the default configuration. This suggests that practical tuning not only improves phased output quality, but also reduces the difficulty of the residual wMEC-style problem faced by the solver.

Based on the optimization sweeps, the recommended practical settings for this pipeline are:

- **Recommended `max_coverage`:** `8`  
- **Recommended phasing thresholds:** `min_mapq = 20`, `min_baseq = 10`  
- **Recommended calling thresholds:** `call_min_mapq = 20` (or leave at default, as the effect is negligible), `call_min_baseq = 10`

A slightly more runtime-biased alternative (`max_coverage = 6`, with the same caller/phasing quality thresholds) performs almost identically on the main accuracy metrics and may be acceptable when efficiency is prioritized.

Future work could compare the custom adaptor against the standard WhatsHap front-end directly, or optimize the adaptor separately if pipeline-level runtime becomes a primary objective.