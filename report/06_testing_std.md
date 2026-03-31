## 6. Testing and Validation (STD-style)

This section describes the testing strategy, test environment assumptions, implemented test cases, and validation evidence for the benchmarking platform. The goal is to ensure correctness and reproducibility of (i) core data structures and metrics, (ii) WhatsHap integration, and (iii) the end-to-end pipeline orchestration and reporting.

### 6.1 Testing objectives

The test plan is designed to verify:

1. **Correctness of core representations and I/O**
   - reads matrices and conversions (TSV/NPZ)
   - adapter correctness when bridging project data to WhatsHap data structures

2. **Correctness of phasing outputs**
   - phased haplotype TSV format
   - phased VCF output formatting (`GT` and `PS`)
   - use of **vendored WhatsHap core** (provenance checks)

3. **Correctness of evaluation metrics**
   - block-flip–invariant accuracy
   - switch error rate
   - phase set counting

4. **End-to-end reproducibility**
   - seed-controlled generation
   - deterministic output contracts (prefix-based artifacts)
   - stable reporting (`*.pipeline.json` → `aggregate.csv`)

### 6.2 Test levels and scope

The project uses a layered testing approach:

- **Unit tests**
  - Focus on small, deterministic components (I/O parsing, metric computations, adapters).

- **Integration tests**
  - Exercise CLIs and drivers with minimal synthetic inputs (NPZ/VCF mode).
  - Verify required outputs and sanity constraints.

- **System / end-to-end tests (smoke)**
  - Use the experiment/pipeline runner to validate the real long-read workflow when the external toolchain is installed.

This separation is necessary because the long-read pipeline depends on external binaries (`minimap2`, `samtools`, `bcftools`) that may not be available in all test environments.

### 6.3 Test environment and dependencies

#### 6.3.1 Python environment
- Python 3.11
- `pytest` for test execution

#### 6.3.2 Vendored WhatsHap core
- Tests that invoke WhatsHap drivers require the vendored `vendor/whatshap_core` to be importable as `whatshap`.
- Where appropriate, tests skip when `whatshap.core` is not importable.

#### 6.3.3 External toolchain (system tests)
System/end-to-end long-read runs require:
- `minimap2`
- `samtools`
- `bcftools`

These are treated as “system dependencies” and are validated at runtime by the pipeline runner.

### 6.4 Implemented unit tests (current repository)

The following unit-level tests are included under `tests/`:

- **`tests/test_whatshap_adapter.py`**
  - Validates compatibility between the project’s read representation and WhatsHap `ReadSet` requirements.
  - Ensures variant indexing and allele encoding match expected semantics.

- **`tests/test_diploid_whatshap.py`**
  - Minimal correctness test for the diploid WhatsHap NPZ-mode driver on a toy matrix.
  - Verifies output artifact creation (`haplotypes.tsv`, `assignments.tsv`, `summary.json`) and basic allele sanity.

- **`tests/test_pedigreedptable_vcf_mode.py`**
  - Regression test for VCF-mode phasing using the WhatsHap DP solver path.
  - Asserts that heterozygous sites become phased (`0|1` / `1|0`) and share a consistent `PS` for a connected component.

- **`tests/test_simulation_whatshap_compat.py`, `tests/test_data_simulation_io.py`**
  - Legacy / matrix-mode tests validating sparse/dense semantics and older simulation I/O.
  - These provide regression coverage for earlier components that remain in the codebase.

### 6.5 Implemented integration tests (current repository)

Integration tests execute real CLIs to validate routing and artifact contracts:

- **`tests/test_vcf_mode_pipeline.py`**
  - End-to-end regression on the matrix-mode pipeline:
    - `dataset.simulate` → `algorithms.cli.phase diploid-whats` (VCF-mode) → output VCF/TSV checks
  - Verifies the phased VCF contains phased heterozygous records and that reconstructed haplotypes reach a conservative accuracy threshold.

- **`tests/test_pipeline_whatshap_cli.py`**
  - Exercises the unified phasing CLI with `diploid-whats`.
  - Verifies the driver uses the **vendored** WhatsHap module by checking paths recorded in `summary.json`.

- **`tests/test_pipeline_integration.py`**
  - Exercises an older end-to-end flow (simulate → convert → phase → evaluate accuracy JSON) for regression coverage.

### 6.6 System / end-to-end smoke validation (long-read pipeline)

Because the long-read pipeline relies on external tools and is more expensive than unit tests, it is validated as a **smoke test** using the orchestrators:

- **Single-run orchestrator**
  - `python -m benchmark.longread_pipeline_runner --prefix ...`
  - Validates the full practical path: reference → truth → read simulation → BAM → VCF → WhatsHap phasing → evaluation → `pipeline.json`.

- **Experiment driver**
  - `python -m benchmark.experiment_driver --outdir ... --seeds ...`
  - Validates:
    - repeated runs across seeds
    - correct aggregation from `*.pipeline.json` to `aggregate.csv`
    - plot generation under `plots/`

**Pass criteria for smoke validation**
- All expected artifacts for the selected regime exist (oracle/called/both).
- `*.pipeline.json` is produced and parseable.
- `*.eval.json` contains required keys (e.g., `switch_error_rate`, `num_phase_sets`).
- Aggregation produces a non-empty `aggregate.csv` with expected columns.

### 6.7 Validation for realism knobs (correctness criteria)

Realism knobs are validated using both metadata checks and metric-based sanity checks:

- **Duplications (reference stressor)**
  - `*.ref.meta.json` contains a `duplications` list with valid, non-overlapping coordinate pairs.
  - Expected effect: decreased calling recall and/or increased phase set fragmentation under higher duplication counts.

- **Coverage dropout (read start model)**
  - `*.reads.meta.json` records dropout parameters.
  - Expected effect: more phase sets (fragmentation) and reduced shared heterozygous recall in called regime.

- **Correlated error bursts**
  - Read simulator parameters are recorded in `*.reads.meta.json`.
  - Expected effect: reduced calling recall or increased unphased rate due to fewer reliable allele observations.

- **Truth indels**
  - Indel introduction is recorded in `*.truth.meta.json`.
  - To keep evaluation meaningful, SNP-only phasing/evaluation is used (`--phase-snps-only`, `--eval-snps-only`) when indels are enabled.
  - Expected effect: if SNP-only policy is applied, SNP phasing metrics should remain stable across indel severity sweeps (within noise), indicating indels are not corrupting SNP evaluation.

### 6.8 Coverage limitations and recommended additional tests

**Current limitations**
- The repository’s automated `pytest` suite primarily covers matrix-mode and driver correctness; the long-read toolchain workflow is validated via smoke runs (requires external binaries).

**Recommended additions (future work / if time permits)**
1. **Long-read runner smoke test (skipped if tools missing)**
   - Run `benchmark.longread_pipeline_runner` with very small parameters (oracle-only) and assert `pipeline.json` schema.

2. **Schema/contract unit tests**
   - Validate that `pipeline.json` contains required fields and that aggregation produces mandatory columns.

3. **Reference meta integrity tests**
   - Deterministically generate a reference with duplications and assert coordinates are consistent (0-based half-open, within bounds, non-overlapping given min-gap).

4. **SNP-only policy regression**
   - Enable indels and assert that enabling `eval_snps_only` prevents spuriously low block-flip accuracy caused by indel representation mismatches.

These additions would strengthen automated coverage for the long-read realism pipeline without requiring large experiment runtimes.