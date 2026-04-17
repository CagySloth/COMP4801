## 6. Testing and Validation (STD-style)

The testing strategy and validation evidence for the platform, with the aim of achieving correctness and reproducibility of data generation, WhatsHap integration, evaluation, and end-to-end pipeline execution, will be covered by this section.

### 6.1 Testing objectives

The testing plan aims to validate these follow characteristics:

1. **Correctness of core representations and I/O**
   - matrix representation and TSV/NPZ conversions
   - correctness of adapter(s) for converting data representations between project-specific format to WhatsHap-compatible inputs

2. **Correctness of phasing outputs**
   - phased output formatting
   - vendored WhatsHap core compatiblity
   - correct propagation of phase set and genotype (PS/GT) information

3. **Correctness of evaluation metrics**
   - block-flip-invariant phase accuracy
   - switch error rate
   - phase-set counting and shared-site metrics

4. **End-to-end reproducibility**
   - seed-controlled execution
   - deterministic prefix-based artifact naming and file generation
   - standardized reporting from `*.pipeline.json` to `aggregate.csv`

### 6.2 Test levels and scope

A layered testing strategy is implemented:

- **Unit tests**
  Small deterministic components such as parsers, metric calculations, adapters, and aggregation logic are covered.

- **Integration tests**
  CLIs and drivers with minimal synthetic inputs and basic sanity checks are tested to verify output contracts (naming-system and file generation).

- **System / end-to-end smoke tests**
  The full long-read workflow is validated through the pipeline runner and experiment driver when external tools are available.

This separation is necessary because the practical long-read pipeline depends on external binaries (`minimap2`, `samtools`, `bcftools`) that may not be present in every test environment.

### 6.3 Test environment and dependencies

#### 6.3.1 Python environment

- Python 3.11
- `pytest` for test execution

#### 6.3.2 Vendored WhatsHap core

Tests that invoke the WhatsHap driver require the vendored `vendor/whatshap_core` package to be importable as `whatshap`. Where appropriate, such tests are skipped when `whatshap.core` is unavailable.

#### 6.3.3 External toolchain

System-level long-read runs require:

- `minimap2`
- `samtools`
- `bcftools`

These are treated as system dependencies and are validated at runtime by the pipeline runner.{index=8}

### 6.4 Implemented tests

The repository includes unit- and integration-level tests covering the main correctness-critical parts of the platform.

#### 6.4.1 Adapter and WhatsHap-compatibility tests

`tests/test_whatshap_adapter.py` validates compatibility between the project’s internal read representation and WhatsHap `ReadSet` requirements. In particular, it checks that variant indexing and allele encoding match the expected semantics used by the vendored phasing backend.

#### 6.4.2 Evaluation and aggregation tests

Core reporting and aggregation logic is validated through tests that ensure run-level reports can be converted into standardized `aggregate.csv` outputs with the expected columns and derived metrics. This is important because downstream plots and tables depend on these machine-readable contracts rather than on log parsing.

#### 6.4.3 Driver and phasing-output sanity tests

Integration-level checks validate that phasing drivers and CLIs produce the required output artifacts and that phased results satisfy basic structural expectations, including presence of phase-set information and well-formed summary/evaluation JSON outputs. These tests are intended to detect interface breakage early without requiring full experiment runs.

### 6.5 Smoke validation of the long-read pipeline

Because the full long-read workflow is more expensive and depends on external tools, it is validated as a smoke test rather than as part of the smallest unit-test layer.

Two orchestrators are used:

- **Single-run pipeline runner**
  `python -m benchmark.longread_pipeline_runner --prefix ...`

- **Experiment driver**
  `python -m benchmark.experiment_driver --outdir ... --seeds ...`

These smoke runs validate the practical path:

**reference → truth → read simulation → BAM → VCF → WhatsHap phasing → evaluation → `pipeline.json`**

Pass criteria include:

- expected artifacts exist for the selected regime
- `*.pipeline.json` is produced and parseable
- `*.eval.json` contains required keys such as switch error and phase-set count
- aggregation produces a non-empty `aggregate.csv` with expected columns

### 6.6 Validation of realism knobs

The realism knobs are validated using both metadata checks and metric-based sanity checks.

- **Duplications**
  - `*.ref.meta.json` must contain valid duplication coordinates
  - expected effect: reduced calling recall and/or increased fragmentation as duplication severity increases

- **Coverage dropout**
  - dropout parameters must be recorded in read metadata
  - expected effect: more phase sets and reduced shared-heterozygous recall in the called regime

- **Correlated error bursts**
  - burst parameters must be recorded in read metadata
  - expected effect: reduced calling quality and/or weaker phasing evidence

- **Truth indels**
  - indel introduction must be recorded in truth metadata
  - when SNP-only policy is enabled, SNP phasing metrics should remain broadly stable across indel sweeps, indicating that indel representation mismatch is not corrupting evaluation

This validation layer is important because the realism knobs are central to the benchmarking contribution of the project, not merely optional simulator features.

### 6.7 Coverage limitations and recommended additional tests

The current automated suite focuses mainly on matrix-mode, adapter, driver, and aggregation correctness, while the long-read toolchain workflow is validated through smoke runs that depend on external binaries.

Recommended additions include:

1. a very small long-read runner smoke test that is skipped automatically when external tools are unavailable
2. schema/contract tests for `pipeline.json` and `aggregate.csv`
3. deterministic reference-meta integrity tests for duplication coordinates
4. explicit SNP-only policy regression tests under indel-enabled truth generation

These additions would strengthen automated coverage further, but the current testing strategy is already sufficient to support the project’s main claims about correctness, reproducibility, and controlled experiment execution.