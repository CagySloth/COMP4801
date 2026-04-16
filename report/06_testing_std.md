## 6. Testing and Validation (STD-style)

This section describes the testing strategy and validation evidence for the benchmarking platform. The goal is to establish correctness and reproducibility of core data handling, WhatsHap integration, evaluation metrics, and the end-to-end long-read pipeline. :contentReference[oaicite:2]{index=2}

### 6.1 Testing objectives

The testing plan is designed to verify four properties:

1. **Correctness of core representations and I/O**
   - matrix and TSV/NPZ conversions
   - adapter correctness when bridging project data structures to WhatsHap-compatible inputs

2. **Correctness of phasing outputs**
   - phased output formatting
   - compatibility with the vendored WhatsHap core
   - correct propagation of PS/GT information

3. **Correctness of evaluation metrics**
   - block-flip-invariant phase accuracy
   - switch error rate
   - phase-set counting and shared-site statistics

4. **End-to-end reproducibility**
   - seed-controlled execution
   - deterministic prefix-based artifact contracts
   - stable reporting from `*.pipeline.json` to `aggregate.csv` :contentReference[oaicite:3]{index=3}

### 6.2 Test levels and scope

The project uses a layered testing strategy:

- **Unit tests**  
  Cover small deterministic components such as parsers, metric calculations, adapters, and aggregation logic.

- **Integration tests**  
  Exercise CLIs and drivers with minimal synthetic inputs to verify output contracts and sanity constraints.

- **System / end-to-end smoke tests**  
  Validate the full long-read workflow through the pipeline runner and experiment driver when external tools are available. :contentReference[oaicite:4]{index=4}

This separation is necessary because the practical long-read pipeline depends on external binaries (`minimap2`, `samtools`, `bcftools`) that may not be present in every test environment. :contentReference[oaicite:5]{index=5}

### 6.3 Test environment and dependencies

#### 6.3.1 Python environment

- Python 3.11
- `pytest` for test execution :contentReference[oaicite:6]{index=6}

#### 6.3.2 Vendored WhatsHap core

Tests that invoke the WhatsHap driver require the vendored `vendor/whatshap_core` package to be importable as `whatshap`. Where appropriate, such tests are skipped when `whatshap.core` is unavailable. :contentReference[oaicite:7]{index=7}

#### 6.3.3 External toolchain

System-level long-read runs require:

- `minimap2`
- `samtools`
- `bcftools`

These are treated as system dependencies and are validated at runtime by the pipeline runner. :contentReference[oaicite:8]{index=8}

### 6.4 Implemented tests

The repository includes unit- and integration-level tests covering the main correctness-critical parts of the platform.

#### 6.4.1 Adapter and WhatsHap-compatibility tests

`tests/test_whatshap_adapter.py` validates compatibility between the project’s internal read representation and WhatsHap `ReadSet` requirements. In particular, it checks that variant indexing and allele encoding match the expected semantics used by the vendored phasing backend. :contentReference[oaicite:9]{index=9}

#### 6.4.2 Evaluation and aggregation tests

Core reporting and aggregation logic is validated through tests that ensure run-level reports can be converted into standardized `aggregate.csv` outputs with the expected columns and derived metrics. This is important because downstream plots and tables depend on these machine-readable contracts rather than on log parsing. :contentReference[oaicite:10]{index=10}

#### 6.4.3 Driver and phasing-output sanity tests

Integration-level checks validate that phasing drivers and CLIs produce the required output artifacts and that phased results satisfy basic structural expectations, including presence of phase-set information and well-formed summary/evaluation JSON outputs. These tests are intended to detect interface breakage early without requiring full experiment runs. :contentReference[oaicite:11]{index=11}

### 6.5 Smoke validation of the long-read pipeline

Because the full long-read workflow is more expensive and depends on external tools, it is validated as a smoke test rather than as part of the smallest unit-test layer.

Two orchestrators are used:

- **Single-run pipeline runner**  
  `python -m benchmark.longread_pipeline_runner --prefix ...`

- **Experiment driver**  
  `python -m benchmark.experiment_driver --outdir ... --seeds ...` :contentReference[oaicite:12]{index=12}

These smoke runs validate the practical path:

**reference → truth → read simulation → BAM → VCF → WhatsHap phasing → evaluation → `pipeline.json`**

Pass criteria include:

- expected artifacts exist for the selected regime
- `*.pipeline.json` is produced and parseable
- `*.eval.json` contains required keys such as switch error and phase-set count
- aggregation produces a non-empty `aggregate.csv` with expected columns :contentReference[oaicite:13]{index=13}

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
  - when SNP-only policy is enabled, SNP phasing metrics should remain broadly stable across indel sweeps, indicating that indel representation mismatch is not corrupting evaluation :contentReference[oaicite:14]{index=14}

This validation layer is important because the realism knobs are central to the benchmarking contribution of the project, not merely optional simulator features. :contentReference[oaicite:15]{index=15}

### 6.7 Coverage limitations and recommended additional tests

The current automated suite focuses mainly on matrix-mode, adapter, driver, and aggregation correctness, while the long-read toolchain workflow is validated through smoke runs that depend on external binaries. :contentReference[oaicite:16]{index=16}

Recommended additions include:

1. a very small long-read runner smoke test that is skipped automatically when external tools are unavailable
2. schema/contract tests for `pipeline.json` and `aggregate.csv`
3. deterministic reference-meta integrity tests for duplication coordinates
4. explicit SNP-only policy regression tests under indel-enabled truth generation :contentReference[oaicite:17]{index=17}

These additions would strengthen automated coverage further, but the current testing strategy is already sufficient to support the project’s main claims about correctness, reproducibility, and controlled experiment execution.