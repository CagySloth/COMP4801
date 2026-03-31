## 4. System Design (SDD-style)

This section describes the design of the WhatsHap benchmarking platform implemented, focusing on architecture, component responsibilities, interfaces, and key design decisions that support reproducible end-to-end experiments.

### 4.1 Design objectives

The system is designed to:

1. **Benchmark WhatsHap in practical workflows**  
   Evaluate phasing in an end-to-end long-read pipeline (FASTA → FASTQ → BAM → VCF → phased VCF), not only on idealized matrices.

2. **Separate error sources**  
   Support both:
   - **Oracle VCF** (ground-truth variants; isolates phasing behavior), and
   - **Called VCF** (from bcftools; measures end-to-end caller + phaser performance).

3. **Be reproducible and traceable**  
   Every run is seed-controlled and produces a machine-readable report (`*.pipeline.json`) that records parameters, file paths, and metrics.

4. **Be modular and scriptable**  
   Each pipeline stage has a dedicated `python -m ...` entrypoint so it can be run independently or orchestrated by runners/drivers.

5. **Support realism knobs and systematic sweeps**  
   “Realism knobs” (duplications, dropout, bursts, indels) are parameterized and can be swept via an experiment driver.

---

### 4.2 High-level pipeline structure and tool mapping

The system is organized as a sequential pipeline with clearly defined stage boundaries. Each stage has (1) a single responsibility, (2) a primary implementation entrypoint in the repository, and (3) a well-defined input/output artifact contract.

#### 4.2.1 Pipeline stages (conceptual view)

1. **Reference synthesis**  
   Create a reference FASTA with optional complexity/duplications.

2. **Truth synthesis (diploid ground truth)**  
   Create phased truth variants (truth VCF), derived haplotypes (hap1/hap2 FASTA), and an oracle VCF.

3. **Read simulation (ONT-like)**  
   Simulate reads from haplotypes, with configurable length models and error models.

4. **Alignment + preprocessing**  
   Align reads to the reference, producing a sorted/indexed BAM suitable for calling and phasing.

5. **Variant calling**  
   Produce a called VCF from the BAM and reference.

6. **Phasing (WhatsHap)**  
   Phase variants given BAM+VCF using the vendored WhatsHap core, producing phased VCF and summaries.

7. **Evaluation + reporting**  
   Compare predicted phasing to truth and produce JSON metrics; record run-level provenance in a pipeline report.

8. **Experiment orchestration (optional layer)**  
   Run multiple pipeline instances across seeds/knobs, aggregate reports into CSV, and generate plots.

This pipeline structure supports two evaluation regimes:
- **Oracle regime:** stages 1–4–6–7 using the oracle VCF to isolate phasing behavior.
- **Called regime:** stages 1–5–6–7 using the called VCF to measure end-to-end behavior.

#### 4.2.2 Implementation mapping (codebase view)

| Stage | Primary entrypoint in COMP4801_20 | External tools | Key outputs |
|------:|-----------------------------------|----------------|-------------|
| 1 | `python -m dataset.longread.reference` | — | `*.ref.fasta`, `*.ref.meta.json` |
| 2 | `python -m dataset.longread.truth` | `bcftools view/index` (bgzip+index) | `*.truth.vcf.gz`, `*.oracle.vcf.gz`, `*.hap1.fasta`, `*.hap2.fasta` |
| 3 | `python -m dataset.longread.readsim` | — | `*.reads.fastq`, `*.reads.truth.tsv`, optional `*.reads.meta.json` |
| 4 | orchestrated by `benchmark.longread_pipeline_runner` | `minimap2`, `samtools sort/index` | `*.bam`, `*.bam.bai` |
| 5 | orchestrated by `benchmark.longread_pipeline_runner` | `bcftools mpileup/call`, `bcftools index` | `*.called.vcf.gz`, `*.called.vcf.gz.tbi` |
| 6 | `python -m algorithms.cli.phase diploid-whats-bam` → `algorithms/diploid/whatshap_bam_driver.py` | (vendored) `whatshap_core` | `*.phased.vcf`, `*.summary.json` |
| 7 | `python -m benchmark.vcf_phase_eval` + report writer in runner | — | `*.eval.json`, `*.pipeline.json` |
| 8 | `python -m benchmark.experiment_driver` + `benchmark.aggregate_pipeline_reports_full` | — | `aggregate.csv`, `plots/*.png` |

The single-run orchestrator `python -m benchmark.longread_pipeline_runner` is the canonical “glue” that runs stages 1–7 and writes `*.pipeline.json`. The experiment driver is a wrapper that repeatedly invokes the single-run orchestrator and then aggregates outputs.

#### 4.2.3 Stage boundaries and responsibilities (quality view)

- **Data generation stages (1–3)** must be deterministic given a seed and must export metadata describing realism knobs applied.
- **Alignment/calling stages (4–5)** are treated as “tool-driven” steps; parameters are surfaced at the pipeline runner layer and recorded into `pipeline.json`.
- **Phasing stage (6)** must use the vendored WhatsHap core to avoid version drift; run summaries record module provenance.
- **Evaluation stage (7)** must be robust to caller differences; when indels are present, the pipeline supports SNP-only phasing/evaluation to avoid representation mismatches.
- **Experiment orchestration stage (8)** must not re-implement logic from stages 1–7; it only schedules runs and aggregates artifacts.

---

### 4.3 Component decomposition and responsibilities

#### 4.3.1 Reference generator (`dataset.longread.reference`)

**Responsibility:** Generate a synthetic reference FASTA with configurable complexity and duplicated segments.

**Outputs:**
- `*.ref.fasta`
- `*.ref.meta.json` containing:
  - `regions` (complex windows such as homopolymers/STR/GC windows)
  - `duplications` (source/destination coordinates for copied segments)
  - `coord_note`: **0-based half-open** intervals `[start0, end0)`

**Key extension points:**
- `--preset {plain,toy,realistic}`
- `--dup-segments`, `--dup-len`, `--dup-min-gap`

---

#### 4.3.2 Truth + haplotype generator (`dataset.longread.truth`)

**Responsibility:** Generate diploid ground truth variants and haplotypes.

**Outputs:**
- `*.truth.vcf` (+ `*.truth.vcf.gz` + index)
- `*.hap1.fasta`, `*.hap2.fasta`
- `*.truth.meta.json`
- Oracle VCF template (written by pipeline runner; see §4.4)

**Notes:**
- The truth VCF is **phased** when `--phased-truth` is used.
- Indels can be introduced via `--num-indels` and indel length parameters.
- If `--avoid-regions` is enabled, truth generation can avoid complex reference regions using `*.ref.meta.json`.

---

#### 4.3.3 Long-read simulator (`dataset.longread.readsim`)

**Responsibility:** Simulate long reads from haplotypes and emit FASTQ with qualities and read-truth metadata.

**Outputs:**
- `*.reads.fastq` (sequence + ASCII qualities)
- `*.reads.truth.tsv` (read-level truth such as haplotype origin; used for debugging/validation)
- optional `*.reads.meta.json` (read simulation parameters/statistics)

**Realism knobs currently supported by the runner:**
- Platform: `--platform {ont,perfect}`
- ONT profile: `--ont-profile {classic,q20}`
- Length distribution: `--len-model {uniform,lognormal}`, `--ln-mean`, `--ln-sigma`
- Start/coverage model: `--start-model {uniform,dropout}`, plus dropout params
- Correlated error bursts: `--burst-prob`, `--burst-count`, `--burst-len`, `--burst-mult`

---

#### 4.3.4 Orchestration runner (`benchmark.longread_pipeline_runner`)

**Responsibility:** Execute the full pipeline for a single run and write a complete run report.

**Workflow steps (as implemented):**
1. Run `dataset.longread.reference`
2. Run `dataset.longread.truth`
3. Create and index truth VCF (`bcftools view/index`)
4. Create oracle VCF from truth (`_write_oracle_vcf`) and bgzip/index (optional)
5. Run `dataset.longread.readsim`
6. Align reads → BAM using `minimap2 | samtools sort`, then `samtools index`
7. Call variants using `bcftools mpileup | bcftools call` (optional)
8. Phase using `algorithms.cli.phase diploid-whats-bam` for oracle/called (optional)
9. Evaluate phasing using `benchmark.vcf_phase_eval`

**Design details:**
- Uses `_which_or_fail()` to ensure external tools exist on `PATH`.
- Uses `_piped()` to reliably run piped commands with error propagation.
- Uses prefix-based naming and writes a single report:
  - `*.pipeline.json` containing:
    - `steps` (file paths)
    - `params` (all run parameters)
    - `callset` (calling precision/recall when applicable)
    - `phasing_runs` (oracle/called subreports with derived metrics)
    - `counts_raw` (record counts for sanity checks)

**Indel-aware evaluation policy:**
- When indels are present, the runner supports:
  - `--phase-snps-only`: filter phasing input VCF to biallelic SNPs
  - `--eval-snps-only`: filter truth VCF to biallelic SNPs for evaluation  
  This prevents evaluation artifacts due to indel representation mismatches.

---

#### 4.3.5 Unified phasing CLI (`algorithms.cli.phase`)

**Responsibility:** Provide a single CLI namespace for all phasing algorithms.

**Relevant command for this project:**
- `diploid-whats-bam`: WhatsHap-core phasing from BAM+VCF using vendored core.

This is the phasing entrypoint called by the pipeline runner.

---

#### 4.3.6 WhatsHap BAM driver (`algorithms.diploid.whatshap_bam_driver`)

**Responsibility:** Implement WhatsHap-style phasing on BAM+VCF inputs using the **vendored** WhatsHap core.

**Pipeline inside the driver:**
1. Parse VCF (minimal parser) and keep **biallelic SNP records**; identify heterozygous sites.
2. Build WhatsHap `ReadSet` from BAM pileups at heterozygous SNP sites:
   - filters by `--min-mapq` and `--min-baseq`
   - only keeps reads covering **≥2 variants** (required by WhatsHap read selection)
3. Perform read selection using `whatshap.readselect.readselection` with `--max-coverage`.
4. Solve phasing using one of:
   - `core.PedigreeDPTable` (WhatsHap default-like solver), or
   - `core.HapChatCore` (optional)
5. Derive `PS` phase-set labels using connectivity of selected reads.
6. Write phased VCF by replacing GT and adding PS.
7. Emit auxiliary outputs:
   - `*.haplotypes.tsv`
   - `*.assignments.tsv` (placeholder assignments; consistent with project output format)
   - `*.summary.json` (counts + module provenance)

**Vendored WhatsHap enforcement:**
- `algorithms.vendor.whatshap_vendor.import_whatshap_vendor()` ensures `import whatshap` resolves to `vendor/whatshap_core`, preventing version drift and keeping benchmarking controlled.

---

#### 4.3.7 Evaluation (`benchmark.vcf_phase_eval`)

**Responsibility:** Compare phased VCF to truth VCF and compute phasing metrics.

**Outputs:**
- `*.eval.json` containing:
  - record overlap counts (shared SNPs, shared hets, phased hets)
  - `num_phase_sets`
  - `phase_accuracy_blockflip`
  - switch error numerator/denominator and rate
  - per-phase-set breakdown (`phase_sets`)

**Assumptions:**
- Evaluation focuses on **biallelic SNPs**, and SNP-only filtering is used when indels are enabled.

---

#### 4.3.8 Aggregation (`benchmark.aggregate_pipeline_reports_full`)

**Responsibility:** Convert many `*.pipeline.json` reports into a single `aggregate.csv` with standardized columns for plotting and analysis.

This is the “data contract” between experiment execution and plotting.

---

#### 4.3.9 Experiment driver (`benchmark.experiment_driver`)

**Responsibility:** Run predefined experiment suites across multiple seeds, then aggregate and plot results.

**Key design features:**
- Runs the pipeline runner repeatedly with controlled argument sets.
- Skips existing runs unless `--force` is used.
- Automatically enables SNP-only phasing/evaluation when indels are configured.
- Produces:
  - `aggregate.csv` per experiment directory
  - `plots/*.png` (error-bar plots across seeds)
  - optional decomposition CSVs (caller vs phaser bottlenecks)

---

### 4.4 Data contracts and file-level interfaces

The design relies on stable file formats exchanged between stages:

- **FASTA**: reference and haplotypes (`*.ref.fasta`, `*.hap1.fasta`, `*.hap2.fasta`)
- **FASTQ**: reads with qualities (`*.reads.fastq`)
- **BAM**: sorted + indexed alignments (`*.bam`, `*.bam.bai`)
- **VCF**:
  - truth (`*.truth.vcf.gz`), oracle (`*.oracle.vcf.gz`), called (`*.called.vcf.gz`)
  - phased outputs (`*.ws*.phased.vcf`)
- **JSON**:
  - `*.ref.meta.json` (regions + duplications; 0-based half-open coords)
  - `*.truth.meta.json` (truth generation parameters/summary)
  - `*.ws*.summary.json` (phasing run summary)
  - `*.ws*.eval.json` (phasing evaluation)
  - `*.pipeline.json` (single-run full report)

Prefix-based naming ensures that all artifacts from a run can be discovered from a single identifier.

---

### 4.5 Realism knobs as configurable stressors

The platform treats “realism knobs” as **configurable stressors**: each knob is a controlled perturbation that (a) approximates a real-world difficulty mode, and (b) is expected to measurably affect either calling metrics, phasing metrics, or both.

These knobs are intentionally:
- **parameterized** (CLI flags),
- **stage-local** (implemented in the stage where the phenomenon arises), and
- **fully recorded** in `*.pipeline.json` (so aggregated results are traceable).

#### 4.5.1 Stressor catalogue and where it applies

| Stressor (knob) | Pipeline stage | Why it matters | CLI / recorded fields |
|---|---|---|---|
| Reference complexity preset | Reference | Creates hard contexts (repeats, GC bias) affecting mapping and calling | `--ref-preset`; `P.ref.meta.json` regions; `pipeline.json.params.ref_preset` |
| Duplicated segments | Reference | Introduces ambiguous mapping → lower call recall / fragmented phase sets | `--dup-segments`, `--dup-len`, `--dup-min-gap`; `P.ref.meta.json.duplications` |
| Avoid complex regions for truth | Truth | Controls whether variants appear inside hard regions | `--avoid-regions` (uses `P.ref.meta.json`) |
| Truth indels | Truth | Indel representation mismatch affects both calling and evaluation if not controlled | `--num-indels`, `--indel-min-len`, `--indel-max-len`, `--indel-het-rate` |
| Read length distribution | Readsim | Affects effective coverage uniformity and phase connectivity | `--len-model`, `--ln-mean`, `--ln-sigma`; `P.reads.meta.json` |
| Coverage dropout model | Readsim | Creates local low-coverage gaps → phase block fragmentation | `--start-model dropout`, `--dropout-fraction`, `--dropout-block-len` |
| Correlated error bursts | Readsim | Models read segments with elevated errors → call errors / read selection changes | `--burst-prob`, `--burst-count`, `--burst-len`, `--burst-mult` |
| SNP-only phasing/eval | Runner / Eval | Prevents indel encoding from dominating evaluation | `--phase-snps-only`, `--eval-snps-only` |

#### 4.5.2 Measurement expectations (what changes when knobs increase)

- **Caller-focused effects** (primarily in called regime):
  - `call_recall` decreases under duplications / dropout / bursts
  - `precision` may vary depending on caller thresholds and error severity

- **Phaser-focused effects** (oracle regime isolates these):
  - `num_phase_sets` increases under dropout and duplications (fragmentation)
  - switch error can increase when reads contain bursty errors and fewer reliable alleles survive filtering

- **End-to-end effects** (called regime):
  - “effective phased recall” typically decreases the most because it is limited by both calling overlap and phasing completeness

These expectations justify treating realism knobs as stressors that can be validated independently before running optimization experiments.

---

### 4.6 Key design decisions (rationale)

1. **Prefix-based artifact naming**  
   Simplifies orchestration, cleanup, and aggregation; avoids implicit state.

2. **Oracle vs called regimes**  
   Enables attribution of performance loss to calling vs phasing.

3. **Vendored WhatsHap core**  
   Prevents dependency drift and ensures controlled benchmarking; provenance is recorded in summary JSON.

4. **SNP-only phasing/evaluation when indels are present**  
   Avoids false failures due to indel representation differences across truth/called/phased VCFs.

5. **Read filtering for WhatsHap compatibility**  
   Reads covering fewer than 2 variants are excluded because WhatsHap read selection requires ≥2 covered variants.

6. **Single-run JSON report as the source of truth**  
   Aggregation and plotting rely on `*.pipeline.json`, not on parsing logs.

---

### 4.7 Extensibility

The design supports incremental realism and experimentation:

- **New realism knobs** can be added in `dataset.longread.reference/truth/readsim`, then surfaced via:
  - `benchmark.longread_pipeline_runner` (single-run support)
  - `benchmark.experiment_driver` (sweep support)
  - `benchmark.aggregate_pipeline_reports_full` (CSV column support)

- **Alternative phasers** can be added by:
  - implementing a new driver under `algorithms/`
  - registering a new subcommand in `algorithms.cli.phase`
  - adding calls and metrics extraction in the pipeline runner/driver
```

**Brief explanation (what this section covered):**
This design section explains *how the system is structured* (architecture + components) and *how data flows* through your end-to-end pipeline. It also documents key design decisions (oracle vs called, vendored WhatsHap, SNP-only policy for indels, prefix-based artifacts) and describes how experiments are automated and aggregated.

**How it aligns with your codebase (`COMP4801_20.zip`):**
Every component described is mapped to an actual module and CLI in your repo: `dataset.longread.reference/truth/readsim`, the orchestrator `benchmark.longread_pipeline_runner`, the experiment automation `benchmark.experiment_driver`, the phasing entrypoint `algorithms.cli.phase diploid-whats-bam` and its implementation `algorithms/diploid/whatshap_bam_driver.py`, plus evaluation/aggregation scripts under `benchmark/`. The file contracts listed match the prefix-based artifacts your runner writes (e.g., `*.pipeline.json`, `*.ws*.eval.json`, `*.ref.meta.json`).
