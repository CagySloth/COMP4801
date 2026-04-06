## 10. Experimental Evaluation

### 10.1 Experimental methodology and common setup

#### 10.1.1 Environment and toolchain
- OS/CPU: MacOS/M2 Max
- Python version: 3.12.12
- minimap2 versions: 2.30-r1287
- samtools version: 1.23
- bcftools version: 1.23
- Vendored WhatsHap core path (from `*.ws*.summary.json`): <TODO>

#### 10.1.2 Pipeline regimes
- **Oracle regime:** phase using `*.oracle.vcf.gz` (phaser-limited)
- **Called regime:** phase using `*.called.vcf.gz` (end-to-end)
- **Indel policy:** when indels are enabled, use SNP-only phasing/evaluation (`phase_snps_only`, `eval_snps_only`)

#### 10.1.3 Fixed parameters (baseline constants)

Unless stated otherwise in a subsection, experiments use the baseline defaults defined by the experiment driver. These parameters control each stage of the end-to-end pipeline (reference/truth generation → read simulation → alignment → calling → phasing → evaluation).

**Reference generation**
- **`ref_length`**: Length of the synthetic reference genome (bp). Larger references reduce boundary effects and allow realism stressors (e.g., duplications) to be placed without overlap.
- **`ref_preset`**: Reference “complexity preset” (e.g., `plain` vs `realistic`). Presets control the presence of hard regions (e.g., repetitive or biased regions) and are recorded in `*.ref.meta.json`.
- **`dup_segments`**: Number of duplicated segments introduced into the reference. Higher values create more ambiguous mapping, typically reducing calling recall and increasing phasing fragmentation.
- **`dup_len`**: Length (bp) of each duplicated segment copied from one locus to another.
- **`dup_min_gap`**: Minimum separation (bp) enforced between duplicated segments/placements to avoid trivial overlaps.

**Truth generation**
- **`num_snps`**: Number of SNP sites inserted into the truth set. Controls variant density and the potential number of heterozygous loci for phasing.
- **`het_rate`**: Fraction of truth SNPs that are heterozygous (e.g., 0/1). WhatsHap only needs to phase heterozygous sites; homozygous sites do not contribute phasing difficulty.
- **`phased_truth`**: Whether the truth VCF contains phased GT (e.g., `0|1`). This is required for switch-error-based evaluation.
- **`random_phase`**: Whether truth phase orientation is randomly assigned per site (used to avoid biasing truth haplotypes).
- **Indel parameters (only when enabled)**:
  - **`num_indels`**: Number of indel variants introduced into truth.
  - **`indel_min_len`, `indel_max_len`**: Indel length range (bp). Larger indels can increase representation differences between truth/called sets.
  - **`indel_het_rate`**: Fraction of indels that are heterozygous.

**Read simulation**
- **`platform`**: Read simulation platform (e.g., `ont`). Determines error profile assumptions.
- **`ont_profile`**: Preset for substitution/insertion/deletion rates and base-quality characteristics (e.g., `q20`). Used to approximate practical ONT error rates.
- **`num_reads`**: Number of reads simulated. This acts as a proxy for sequencing depth in this synthetic setting.
- **`min_len`, `max_len`**: Bounds on read lengths (bp). Affects overlap connectivity (and therefore phasing block sizes).
- **`hap1_frac`**: Fraction of reads sampled from haplotype 1 (vs haplotype 2). Typically 0.5 for balanced coverage.
- **Length model parameters**:
  - **`len_model`**: Distribution used to draw read lengths (e.g., `uniform` or `lognormal`).
  - **`ln_mean`, `ln_sigma`**: Parameters of the lognormal model (used when `len_model=lognormal`).
- **Start/coverage model parameters**:
  - **`start_model`**: How read start positions are sampled (e.g., `uniform` vs `dropout`).
  - **`dropout_fraction`**: Fraction of the reference affected by coverage dropout (low/near-zero coverage windows).
  - **`dropout_block_len`**: Typical length of dropout blocks (bp), controlling the size of connectivity gaps.
- **Correlated error burst parameters**:
  - **`burst_prob`**: Probability a read contains one or more “bursts” of elevated error.
  - **`burst_count`**: Number of bursts injected into a read (often fixed or bounded).
  - **`burst_len`**: Length of each burst (bp).
  - **`burst_mult`**: Multiplier applied to base error rates within bursts (higher → more severe).

**Alignment and calling thresholds**
- **`map_preset`**: Minimap2 mapping preset (e.g., `map-ont`). Chosen to match ONT alignment behavior.
- **`call_min_mapq`**: Minimum mapping quality for reads to be included in variant calling (`bcftools mpileup -q`). Higher thresholds reduce ambiguous mappings but may reduce coverage.
- **`call_min_baseq`**: Minimum base quality for bases to contribute to variant calling (`bcftools mpileup -Q`). Higher thresholds reduce noisy evidence but may reduce sensitivity.

**Phasing parameters (WhatsHap integration)**
- **`vcf_source`**: Which VCF to use for phasing:
  - `oracle`: truth-derived variant set (isolates phasing),
  - `called`: called variant set (end-to-end),
  - `both`: run both for attribution.
- **`max_coverage`**: WhatsHap read-selection cap (maximum effective coverage retained). Larger values retain more reads (potentially better connectivity) but increase compute.
- **`min_mapq`**: Minimum mapping quality for reads to be used in phasing (driver filter). Controls how aggressively ambiguous alignments are excluded.
- **`min_baseq`**: Minimum base quality for allele observations used in phasing (driver filter). Controls how aggressively noisy allele calls are excluded.
- **Indel-safe evaluation policy flags**:
  - **`phase_snps_only`**: Filter *input to phasing* to biallelic SNPs only.
  - **`eval_snps_only`**: Filter *truth used for evaluation* to biallelic SNPs only.
  These flags prevent indel representation differences from corrupting SNP phasing evaluation.

#### 10.1.4 Metrics reported (from `aggregate.csv`)

Metrics are computed per run from `*.pipeline.json` and `*.eval.json`, then aggregated across seeds into `aggregate.csv`. We report both **oracle** (phaser-focused) and **called** (end-to-end) metrics to attribute performance loss to variant calling vs phasing.

**Calling quality (called regime only)**
- **`call_precision`**: Fraction of called SNPs that match truth SNP sites (low precision indicates false positives).  
  Intuition: “How many called variants are correct?”
- **`call_recall`**: Fraction of truth SNPs recovered by the caller (low recall indicates missed variants).  
  Intuition: “How many true variants were found?”
- **`truth_snps`**: Number of SNP records in the truth set (after any SNP-only filtering).
- **`called_snps`**: Number of SNP records produced by the caller.
- **`shared_snps`**: Number of SNP records present in both truth and called sets (the overlap the phaser can actually work with in called regime).

**Phasing coverage / completeness**
- **`*_shared_het_recall`** (called): Fraction of truth heterozygous SNPs that appear in the phasing input overlap (i.e., how many true hets survive calling and exact matching).  
  This captures the “caller bottleneck” for phasing.
- **`*_phasing_rate_shared_het`** (called): Fraction of shared heterozygous SNPs that end up phased (i.e., GT becomes `0|1` or `1|0`).  
  This captures the “phaser completeness” *given* available het sites.
- **`*_effective_phased_recall`** (oracle/called): Fraction of truth heterozygous SNPs that are both (a) present in the evaluation set and (b) phased correctly (under best block flip).  
  This is the main “end-to-end usable phasing” metric.

**Phasing correctness**
- **`*_phase_accuracy`** (oracle/called): Accuracy of phased heterozygous genotypes under an optimal global flip per phase set.  
  Reason: haplotypes are label-symmetric (swap hap1/hap2 yields equivalent solution), so evaluation must be flip-invariant.
- **`*_switch_error`** (oracle/called): Switch error rate, computed as `#switches / #adjacent_het_pairs_compared`.  
  Interpretation: how often the predicted phase relationship between adjacent heterozygous variants is inconsistent with truth.  
  Switch error is often the most informative correctness metric in long-read phasing.

**Phase block fragmentation**
- **`*_num_phase_sets`** (oracle/called): Number of distinct phase sets (`PS`) in the phased VCF.  
  Lower is better (fewer blocks means longer haplotype continuity), but extremely low values can also occur if very few het sites are phased.
- **Phase set size distribution (from eval JSON)**: Per-phase-set sizes and best-flip matches.  
  Used to interpret whether fragmentation is driven by coverage gaps (dropout), ambiguous mapping (duplications), or filtering thresholds.

**Runtime / efficiency**
- **`time_total_sec`**: Total pipeline runtime per run (if recorded).  
  Used in optimization sweeps to quantify trade-offs between accuracy/continuity and compute.

**Oracle vs called interpretation**
- Oracle metrics represent “phaser-limited” performance (variants are correct).
- Called metrics represent “end-to-end” performance (caller + phaser).
- The gap between oracle and called effective phased recall indicates how much performance is lost due to variant calling overlap rather than phasing itself.

---

### 10.2 Baseline: depth scaling (coverage proxy via num_reads)

#### Aim
Establish baseline performance curves as sequencing depth increases, and attribute performance limitations by comparing the oracle (phaser-limited) regime against the called (end-to-end) regime.

#### Setup
- Varied: `num_reads ∈ {50, 100, 200, 400}` (3 seeds each)
- Fixed: ONT-like simulation profile (q20), reference length 80 kb, 800 truth SNPs (het rate 0.8), read length range 2–6 kb, calling thresholds (mapQ/baseQ), and phasing thresholds (max_coverage, min_mapq, min_baseq) as specified in §10.1.

#### Figures to include
- 

**Figure 10.2.1 — Variant calling recall vs sequencing depth.**  
Calling recall is computed as `shared_snps / truth_snps`, where `shared_snps` are SNP records present in both the truth VCF and the called VCF (exact-position matching in the pipeline’s callset summary). As read depth increases (num_reads), recall increases sharply, indicating that *variant discovery* is the main limiter at low coverage. Points show mean over seeds; error bars show standard deviation.

**Figure 10.2.2 — Oracle effective phased recall vs sequencing depth (phasing-only upper bound).**  
Oracle effective phased recall measures WhatsHap’s phasing performance when the input variants are “oracle” (truth-derived SNP sites), removing variant calling errors from the pipeline. The curve rises with depth and approaches a high plateau, representing the phasing-only capability under the simulated ONT Q20 error model. Points show mean over seeds; error bars show standard deviation.

**Figure 10.2.3 — Called effective phased recall vs sequencing depth (end-to-end performance).**  
Called effective phased recall measures the fraction of truth heterozygous SNPs that end up *both present in the called set and correctly phased* (after best block flip). This curve is substantially lower than the oracle curve at low depth, demonstrating that variant calling recall (and callset composition) dominates end-to-end phasing performance early on. Points show mean over seeds; error bars show standard deviation.

**Figure 10.2.4 — Oracle number of phase sets vs sequencing depth (block fragmentation, phasing-only).**  
The number of phase sets (PS blocks) reflects *fragmentation*: more phase sets means shorter blocks / less connectivity. Under oracle variants, increasing depth reduces fragmentation as reads bridge more heterozygous sites. Points show mean over seeds; error bars show standard deviation.

**Figure 10.2.5 — Called number of phase sets vs sequencing depth (block fragmentation, end-to-end).**  
The phase-set count using called variants reflects fragmentation under end-to-end conditions. Unlike the oracle curve, this may be non-monotonic at moderate depths because the called VCF changes with depth (different sites called and filtered), affecting graph connectivity. Points show mean over seeds; error bars show standard deviation.

**Figure 10.2.6 — Called switch error rate vs sequencing depth.**  
Switch error rate is computed over adjacent phased heterozygous SNP pairs within phase sets (after choosing the best global flip per phase set). The mean generally decreases as depth increases, but can show variance at intermediate depths due to small denominators and changing callsets. Error bars show standard deviation; interpret values as bounded to [0, 1] even if `mean ± std` extends below 0.

#### Results summary (means over seeds)
Calling (called regime):
- `call_recall`: 0.224 (50) → 0.399 (100) → 0.662 (200) → 0.877 (400)
- `call_precision`: 1.000 at all depths in this synthetic setup (no SNP false positives observed)

Oracle phasing (phaser-limited):
- `oracle_effective_phased_recall`: 0.513 → 0.750 → 0.921 → 0.969
- `oracle_switch_error`: near 0 across all depths
- `oracle_num_phase_sets`: 10.7 → 5.7 → 1.7 → 1.0

Called phasing (end-to-end):
- `called_effective_phased_recall`: 0.040 → 0.140 → 0.426 → 0.845
- `called_shared_het_recall`: 0.129 → 0.281 → 0.590 → 0.854
- `called_phasing_rate_shared_het`: 0.355 → 0.517 → 0.808 → 0.991
- `called_switch_error`: high variance at 100–200 reads (outlier seeds), near 0 at 400 reads
- `called_num_phase_sets`: 5.0 → 7.7 → 3.0 → 1.0 (non-monotonic due to changing number of phased het sites)

Taken together, these decomposition terms show that most loss in called effective phased recall arises from reduced shared heterozygous site recall rather than from poor phase accuracy, indicating that variant recovery and connectivity are the dominant bottlenecks in the baseline setting.

#### Key observations
- O1 (Calling is depth-limited): Calling recall increases strongly with `num_reads` (0.224 → 0.877), indicating the caller is primarily sensitivity-limited at low/moderate depth. Precision remains 1.0 in this controlled setup, so errors manifest mainly as missed sites rather than spurious sites.
- O2 (WhatsHap is accurate given correct sites): In the oracle regime, switch error is ~0 and phase accuracy is ~1 even at low depth; the main limitation at 50–100 reads is incompleteness (many het sites remain unphased due to connectivity gaps).
- O3 (Residual high-depth loss remains caller-limited): At 400 reads, called phasing is nearly perfect on shared heterozygous sites (phasing rate ≈ 0.99, accuracy ≈ 1.0). The remaining gap between oracle and called effective phased recall is therefore dominated by missing heterozygous sites in the called set (`called_shared_het_recall` ≈ 0.85).
- O4 (Mid-depth instability in called regime): At 100–200 reads, the called regime shows occasional high switch error in individual seeds. This likely reflects both seed sensitivity and small denominators in switch-error calculation when the called VCF contains relatively few informative heterozygous sites. This motivates using more seeds in the final pass.

#### Takeaway
The baseline depth sweep confirms that the pipeline behaves sensibly: calling recall increases with depth; oracle phasing becomes both accurate and contiguous with increasing depth; and end-to-end (called) phasing improves sharply but remains limited primarily by the overlap of correctly called heterozygous variants. At sufficiently high depth (400 reads in this setup), WhatsHap phasing is essentially correct on available het sites, implying that further end-to-end gains would require improvements in variant calling overlap rather than the phasing algorithm itself.

---

### 10.3 Single-knob study: duplicated regions (reference stressor)

#### Aim
Measure how duplicated reference segments affect end-to-end phasing by increasing mapping ambiguity, and determine whether the main impact is on variant recovery, phasing continuity, or phasing correctness.

#### Setup
- Varied: `dup_segments ∈ {0, 1, 3, 5}` (3 seeds each by default; 5 seeds recommended if runtime allows)
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

#### Figures to include
- `call_recall` vs `dup_segments`
- `called_effective_phased_recall` vs `dup_segments`
- `called_switch_error` vs `dup_segments`
- `oracle_effective_phased_recall` vs `dup_segments`
- `called_num_phase_sets` vs `dup_segments` (secondary / interpret with caution)

**Figure 10.3.1 — Variant calling recall vs duplicated-region count.**  
Calling recall remains nearly unchanged across duplication levels, indicating that the duplicated-region stressor does not substantially reduce SNP site recovery under the current alignment and calling settings.

**Figure 10.3.2 — Called effective phased recall vs duplicated-region count.**  
End-to-end phasing remains broadly stable at low-to-moderate duplication levels, but declines at the highest duplication setting. This indicates that duplicated regions only become a substantial stressor in the called regime once ambiguity is sufficiently strong.

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

Taken together, these results show that duplicated regions are a relatively mild stressor under the current parameterization. The strongest degradation appears only at the highest duplication level (`dup_segments = 5`), and is driven mainly by reduced called-regime phasing completeness and phase accuracy rather than by reduced callset overlap.

#### Key observations
- O1 (Calling recall is largely unaffected by duplication in this setup): Across `dup_segments ∈ {0,1,3,5}`, calling recall remains nearly flat (`0.663–0.671`) and call precision stays at 1.0. This indicates that the current duplication stressor does not substantially reduce SNP site recovery under the baseline alignment/calling settings.

- O2 (Oracle phasing remains stable): Oracle effective phased recall remains high and nearly unchanged (`0.918–0.932`), with near-zero oracle switch error throughout. This shows that when correct variant sites are supplied, WhatsHap itself remains robust to the duplicated-region stressor at these settings.

- O3 (Called-regime degradation appears mainly at the highest duplication level): At `dup_segments = 5`, called effective phased recall drops from `0.479` to `0.372`, while called switch error rises from `0.040` to `0.209`. However, `called_shared_het_recall` remains approximately constant (`~0.59`), indicating that the degradation is not primarily due to loss of shared heterozygous sites.

- O4 (Main loss is in phasing completeness/correctness on surviving sites): The highest-duplication condition reduces `called_phasing_rate_shared_het` (`0.848 → 0.778`) and `called_phase_accuracy` (`0.960 → 0.802`). This suggests that duplicated regions mainly weaken the consistency of phasing evidence in the called regime rather than simply reducing callset overlap.

- O5 (Phase-set count is not the main signal here): `called_num_phase_sets` does not increase with duplication, and is actually lower at `dup_segments = 5`. In this experiment, lower phase-set count does not imply better phasing, because it occurs alongside lower effective phased recall and higher switch error. This indicates that phase-set count should be interpreted together with completeness and correctness metrics rather than in isolation.

#### Takeaway
Under the current parameterization, duplicated regions are a relatively mild stressor for variant recovery and for oracle phasing. Their main impact appears only at the highest duplication level, where end-to-end (called) phasing becomes less complete and less accurate despite largely unchanged callset overlap. This suggests that duplicated regions primarily affect the quality and stability of phasing evidence in the called regime, rather than acting as a strong caller-limited bottleneck. In other words, the duplication stressor in this setup is better characterized as a called-regime phasing/evidence-quality stressor than as a pure site-recovery stressor.

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

#### Figures to include
- `oracle_num_phase_sets` vs `dropout_fraction`
- `called_num_phase_sets` vs `dropout_fraction` (secondary; interpret with completeness metrics)
- `oracle_effective_phased_recall` vs `dropout_fraction`
- `called_effective_phased_recall` vs `dropout_fraction`
- `call_recall` vs `dropout_fraction` (supporting evidence for severe-dropout calling loss)

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

Decomposition of called effective phased recall shows that coverage dropout first increases loss through unphased shared heterozygous sites (fragmentation / reduced connectivity), and at the strongest dropout level also causes a substantial increase in loss from missing called sites.

#### Key observations
- O1 (Coverage dropout strongly increases fragmentation): In the oracle regime, the number of phase sets rises from `1.8` at `dropout_fraction = 0.00` to `10.0` at `0.20`, while oracle switch error remains near zero throughout. This shows that dropout primarily breaks phase-block continuity rather than causing incorrect phase assignments.

- O2 (Dropout directly reduces phasing completeness even with correct variant sites): Oracle effective phased recall declines from `0.925` to `0.637` as dropout increases. Because oracle shared-site recall remains fixed by construction, this reduction is attributable to weaker read connectivity and more unphased heterozygous sites.

- O3 (End-to-end performance degrades more strongly than oracle performance): Called effective phased recall falls from `0.436` to `0.081`, a steeper decline than in the oracle regime. This indicates that dropout harms end-to-end phasing through both connectivity loss and reduced overlap of correctly called heterozygous variants.

- O4 (Moderate dropout is mainly a connectivity / completeness problem): From `dropout_fraction = 0.00` to `0.10`, the main worsening in the called regime is a drop in `called_phasing_rate_shared_het` (`0.819 → 0.514`) together with increased fragmentation. This indicates that moderate dropout primarily prevents surviving heterozygous sites from being phased into long connected blocks.

- O5 (Severe dropout becomes a mixed calling + phasing bottleneck): At `dropout_fraction = 0.20`, `call_recall` drops sharply to `0.432` and `called_shared_het_recall` falls to `0.362`, while `called_phasing_rate_shared_het` also falls to `0.239`. This shows that strong dropout not only fragments phase blocks but also substantially reduces variant recovery in the called regime.

- O6 (Switch error is not the main signal in this study): Although called switch error does not increase monotonically, this should not be interpreted as improved phasing at moderate dropout. As dropout increases, fewer heterozygous sites remain phased, making switch error less informative than completeness and fragmentation metrics.

#### Takeaway
Coverage dropout is a strong fragmentation-dominated stressor in this pipeline. Its primary effect is to create low-coverage gaps that prevent reads from connecting heterozygous sites into long phase blocks, which reduces phasing completeness even in the oracle regime. In the called regime, this effect is amplified, and at severe dropout levels the stressor becomes mixed: both phase connectivity and variant recovery deteriorate substantially. Overall, coverage dropout is one of the clearest examples in this study of a realism knob that directly disrupts haplotype continuity.

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
- **Figure 10.X — Calling recall vs burst probability.**  
  Shows whether structured high-error segments reduce SNP recovery in the called regime.

- **Figure 10.X — Called effective phased recall vs burst probability.**  
  Measures the overall end-to-end phasing impact of error bursts.

- **Figure 10.X — Called switch error rate vs burst probability.**  
  Measures whether bursty errors destabilize local phase relationships between adjacent heterozygous sites.

- **Figure 10.X — Called number of phase sets vs burst probability.**  
  Measures whether bursts also increase fragmentation, in addition to correctness loss.

- **Figure 10.X — Oracle effective phased recall vs burst probability.**  
  Provides attribution by showing whether burst-induced degradation persists even when correct variant sites are supplied.

#### Results summary (means over seeds)
- `call_recall`: <TODO>
- `called_effective_phased_recall`: <TODO>
- `called_switch_error`: <TODO>
- `called_num_phase_sets`: <TODO>
- `called_shared_het_recall`: <TODO, if informative>
- `called_phasing_rate_shared_het`: <TODO, if informative>
- `called_phase_accuracy`: <TODO, if informative>
- `oracle_effective_phased_recall`: <TODO>

Taken together, these metrics show whether correlated error bursts primarily act as a calling-quality stressor, a phasing-correctness stressor, a fragmentation stressor, or a mixture of these effects.

#### Key observations
- O1: <Does `call_recall` decrease as `burst_prob` increases?>
- O2: <Does `called_switch_error` increase, indicating reduced phasing stability?>
- O3: <Does oracle performance remain stable, suggesting the main effect is in the called regime, or does oracle also degrade?>
- O4: <Does decomposition show that the main loss is in shared-site recall, phasing rate on shared sites, or phase accuracy?>

#### Takeaway
Correlated error bursts primarily act as a <calling-quality / phasing-correctness / mixed> stressor in this pipeline. Their main effect is <reducing reliable variant evidence / destabilizing phase relationships / both>, as shown by the relative changes in called versus oracle metrics.

---

### 10.6 Single-knob study: read length model

#### Aim
Compare uniform vs lognormal read length distributions on connectivity and phasing.

#### Setup
- Vary: `len_model ∈ {uniform, lognormal}` (and `ln_mean/ln_sigma` if lognormal)
- Keep fixed: baseline constants

#### Commands
- `python -m benchmark.experiment_driver --outdir <...> --seeds ... --only lenmodel`

#### Outputs
- `<outdir>/<lenmodel_section>/aggregate.csv`
- `<outdir>/<lenmodel_section>/plots/*.png`

#### Figures to include
- `called_num_phase_sets` (uniform vs lognormal)
- oracle/called effective phased recall (uniform vs lognormal)

#### Key observations
- O1: <TODO>
- O2: <TODO>

#### Takeaway
- <TODO>

---

### 10.7 Truth indels and SNP-only policy (evaluation robustness)

#### Aim
Show that indel-heavy truth does not invalidate SNP phasing evaluation when SNP-only policy is applied.

#### Setup
- Enable indels in truth generation
- Ensure SNP-only flags are enabled:
  - `phase_snps_only = true`
  - `eval_snps_only = true`
- Vary: indel severity (count and/or length distribution)

#### Commands
- `python -m benchmark.experiment_driver --outdir <...> --seeds ... --only indels`

#### Outputs
- `<outdir>/<indels_section>/aggregate.csv`
- `<outdir>/<indels_section>/plots/*.png`

#### Figures to include
- oracle phasing metrics (SNP-only) vs indel severity:
  - `oracle_switch_error`, `oracle_effective_phased_recall`
- called metrics vs indel severity:
  - `called_effective_phased_recall`, `call_recall`

#### Key observations
- O1: <TODO>
- O2: <TODO>

#### Takeaway
- <TODO>

---

### 10.8 Interaction study: duplications × dropout

#### Aim
Test whether stressors compound non-linearly.

#### Setup
Small grid (example 2×2):
- `dup_segments ∈ {0, high}`
- dropout ∈ {off, on}

#### Commands
- `python -m benchmark.experiment_driver --outdir <...> --seeds ... --only interaction`

#### Outputs
- `<outdir>/<interaction_section>/aggregate.csv`
- `<outdir>/<interaction_section>/plots/*.png`

#### Figures to include
- `called_effective_phased_recall` across conditions (grouped bars or heatmap)
- `called_num_phase_sets` across conditions (grouped bars or heatmap)

#### Key observations
- O1: <TODO>
- O2: <TODO>

#### Takeaway
- <TODO>

---

### 10.9 Optimization sweeps (WhatsHap / pipeline parameters)

#### Aim
Identify parameter settings that improve phasing metrics under hard conditions and quantify trade-offs.

#### Setup
Run on a fixed “hard scenario” (document the chosen stressors), then:
- Sweep A: `max_coverage`
- Sweep B: `min_mapq × min_baseq`

#### Commands
- `python -m benchmark.experiment_driver --outdir <...> --seeds ... --only optimize`

#### Outputs
- `<outdir>/<optimize_section>/.../aggregate.csv`
- `<outdir>/<optimize_section>/plots/*.png`

#### Figures to include
- effective phased recall vs `max_coverage` (+ runtime if available)
- switch error vs quality thresholds (grid)
- phase sets vs quality thresholds (grid)

#### Key observations
- O1: <TODO>
- O2: <TODO>

#### Takeaway (lead-in to Section 11 recommendations)
- <TODO>

---

### 10.10 Summary of results (cross-experiment synthesis)

#### 10.10.1 Attribution summary (oracle vs called)
For each stressor, classify the main bottleneck using oracle vs called gap and decomposition columns:
- Duplications: <caller-limited / phaser-limited / mixed> (evidence: <TODO>)
- Dropout: <...>
- Bursts: <...>
- Indels (SNP-only): <...>

#### 10.10.2 Most impactful stressors
Rank stressors by impact on:
- `called_effective_phased_recall`
- `called_switch_error`
- `called_num_phase_sets`

#### 10.10.3 Optimization summary
Summarize best trade-off settings found in §10.9 (to be expanded in Section 11):
- Recommended `max_coverage`: <TODO>
- Recommended `min_mapq/min_baseq`: <TODO>
- Notes on runtime and trade-offs: <TODO>