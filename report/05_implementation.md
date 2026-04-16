## 5. Implementation

This section describes how the benchmarking platform is implemented, focusing on repository organization, the canonical runner, WhatsHap integration, evaluation logic, and experiment automation.

### 5.1 Repository structure

The repository is organized around three main concerns: data generation, phasing algorithms, and benchmarking/orchestration.

- `dataset/longread/`  
  Implements long-read data generation:
  - `reference.py`: synthetic reference generation with presets and duplicated-region support
  - `truth.py`: diploid truth variants, haplotypes, oracle VCF, and optional indels
  - `readsim.py`: ONT-like read simulation with configurable length, dropout, and burst-error models

- `algorithms/`  
  Implements phasing CLIs and backend drivers:
  - `cli/phase.py`: unified phasing entrypoint
  - `diploid/whatshap_bam_driver.py`: WhatsHap phasing from BAM + VCF using the vendored core
  - vendoring helper(s): ensure imports resolve to the vendored WhatsHap package

- `benchmark/`  
  Implements orchestration, evaluation, aggregation, and experiment control:
  - `longread_pipeline_runner.py`: canonical single-run orchestrator producing `*.pipeline.json`
  - `vcf_phase_eval.py`: phased VCF evaluation against truth
  - `aggregate_pipeline_reports.py`: converts many run reports into `aggregate.csv`
  - `experiment_driver.py`: runs experiment suites across seeds and parameter settings

- `tests/`  
  Unit and integration tests for core pipeline behaviour, evaluation, aggregation, and realism-specific logic

- `vendor/whatshap_core/`  
  Vendored WhatsHap core used for controlled benchmarking and stable integration behaviour

This structure keeps data generation, phasing, and benchmarking concerns separate while preserving a single orchestration path for full runs. :contentReference[oaicite:2]{index=2}

### 5.2 End-to-end orchestration: `benchmark.longread_pipeline_runner`

The canonical workflow is implemented as a single-run pipeline runner:

- Entry point: `python -m benchmark.longread_pipeline_runner`
- Primary argument: `--prefix <path-prefix>`

The runner enforces two main contracts:

1. **Prefix-based artifact naming**  
   All run outputs share a common prefix, making related artifacts easy to discover and aggregate.

2. **Machine-readable run report**  
   Every run produces a `*.pipeline.json` file recording parameters, file outputs, and evaluation metrics.

#### 5.2.1 Pipeline execution flow

For a run with prefix `P`, the runner performs the following high-level steps:

1. generate the synthetic reference and metadata  
2. generate truth variants, haplotypes, and oracle VCF  
3. simulate ONT-like reads  
4. align reads and produce a sorted/indexed BAM  
5. call variants to produce `P.called.vcf.gz`  
6. phase oracle and/or called variants using the configured phasing backend  
7. evaluate phased output against truth  
8. write a unified `P.pipeline.json` report

The runner is therefore the single source of truth for how a complete pipeline instance is executed. :contentReference[oaicite:3]{index=3}

#### 5.2.2 Toolchain validation and robustness

Before running stages that depend on external tools, the runner validates that the required executables are available on `PATH`. This includes:

- `minimap2`
- `samtools`
- `bcftools`

This early validation prevents partially completed runs caused by missing dependencies and improves reproducibility in scripted experiment execution. :contentReference[oaicite:4]{index=4}

### 5.3 WhatsHap integration: `diploid-whats-bam`

WhatsHap-based phasing is exposed through the unified CLI:

- `python -m algorithms.cli.phase diploid-whats-bam`

This entrypoint delegates to `algorithms.diploid.whatshap_bam_driver`, which implements phasing from BAM + VCF using the vendored WhatsHap core. The main phasing parameters surfaced through the driver are:

- `--bam`
- `--vcf`
- `--output-vcf`
- `--output-prefix`
- `--max-coverage`
- `--min-mapq`
- `--min-baseq`

#### 5.3.1 Vendored WhatsHap core

The phasing backend uses the vendored WhatsHap core under `vendor/whatshap_core/` rather than a system-installed copy. This design keeps benchmarking behaviour controlled and reproducible, and avoids dependency drift across environments. It also makes backend behaviour easier to inspect and reason about during controlled experiments. :contentReference[oaicite:5]{index=5}

#### 5.3.2 Read filtering and selection

The driver constructs the `ReadSet` consumed by the vendored core, applies phasing-side filters such as mapping-quality and base-quality thresholds, and enforces compatibility with WhatsHap’s read-selection assumptions. In particular, reads covering fewer than two informative variants are excluded before read selection, since they cannot contribute to haplotype linkage.

This stage is important because the study evaluates not only final phasing output, but also how upstream evidence retention affects end-to-end performance. :contentReference[oaicite:6]{index=6}

#### 5.3.3 Phased VCF output

The driver writes phased VCF output together with a summary JSON containing instrumentation, counts, and phasing-related metadata. These outputs are then consumed by the evaluation stage and folded into the final `*.pipeline.json` report. :contentReference[oaicite:7]{index=7}

### 5.4 Evaluation implementation: `benchmark.vcf_phase_eval`

Phasing evaluation is implemented by `benchmark.vcf_phase_eval`, which compares phased VCF output against truth and writes `*.eval.json`.

The evaluation focuses on three main aspects:

- **correctness**
  - phase accuracy under best block flip
  - switch error rate

- **completeness**
  - shared heterozygous recall
  - phasing rate on shared heterozygous sites
  - effective phased recall

- **fragmentation**
  - number of phase sets
  - phase-set structure implied by PS tags

This evaluator is used for both oracle and called regimes, allowing direct attribution of losses to phasing vs upstream variant recovery. :contentReference[oaicite:8]{index=8}

### 5.5 Indel-aware policy: SNP-only phasing/evaluation flags

When truth indels are enabled, the pipeline supports two important safeguards:

- `--phase-snps-only`
- `--eval-snps-only`

These options restrict phasing input and/or evaluation truth to biallelic SNPs so that SNP phasing metrics remain interpretable even when indel representations differ across truth, called, and predicted VCFs. This policy is a key implementation choice because it prevents representation mismatch from corrupting SNP-focused phasing evaluation. :contentReference[oaicite:9]{index=9}

### 5.6 Experiment automation and aggregation

#### 5.6.1 Experiment driver

Systematic studies are executed through:

- `python -m benchmark.experiment_driver --outdir ... --seeds ...`

The driver repeatedly invokes the single-run pipeline under predefined settings, supports selective execution through `--only`, and produces per-experiment directories containing per-seed outputs, `aggregate.csv`, and plots. This keeps experiment logic separate from single-run orchestration while preserving one canonical execution path. :contentReference[oaicite:10]{index=10}

#### 5.6.2 Aggregation contract

Aggregation is implemented by `benchmark.aggregate_pipeline_reports`, which converts many run-level `*.pipeline.json` files into a single `aggregate.csv`. This CSV acts as the standard interface between experiment execution and downstream plotting/reporting. By aggregating from machine-readable reports rather than logs, the design improves traceability and reproducibility. :contentReference[oaicite:11]{index=11}

### 5.7 Realism knobs: where they are implemented and how they are surfaced

Realism knobs are implemented close to the stage where they conceptually belong, and then surfaced through the pipeline runner and experiment driver.

- **Reference-level knobs**
  - complexity preset
  - duplicated segments

- **Truth-level knobs**
  - optional indels
  - heterozygosity and SNP density

- **Read-level knobs**
  - read length model
  - coverage dropout model
  - correlated error bursts

- **Evaluation-policy knobs**
  - oracle vs called VCF choice
  - SNP-only phasing/evaluation flags

This staged design is important because it allows each realism factor to be varied independently in controlled sweeps while still being recorded consistently in the run report and downstream aggregates. :contentReference[oaicite:12]{index=12}