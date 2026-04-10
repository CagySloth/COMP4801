## 7. Software Configuration Management Plan (SCMP-style)

This section describes how source code, experimental configurations, and generated artifacts are managed to ensure traceability and reproducibility of results reported in this project.

### 7.1 Configuration items (CIs)

The primary configuration items are:

1. **Source code repository**
   - All Python modules under `dataset/`, `algorithms/`, and `benchmark/`
   - Vendored phasing core under `vendor/whatshap_core/`
   - Documentation under `docs/` and project report drafts under `report/` (if applicable)

2. **Experiment definitions**
   - Experiment suites and parameter grids defined in `benchmark/experiment_driver.py`
   - Single-run parameterization defined by `benchmark/longread_pipeline_runner.py`

3. **Run artifacts and reports**
   - Prefix-scoped run artifacts (FASTA/FASTQ/BAM/VCF)
   - Machine-readable run reports: `*.pipeline.json`
   - Aggregated datasets: `aggregate.csv`
   - Plots derived from aggregates: `plots/*.png`

### 7.2 Version control strategy

The project uses Git to track changes to all code and documentation.

Recommended branch roles (reflecting actual workflow in this project):
- `main`: stable baseline used to generate report results and final experiments
- `feature/*`: isolated development of new realism knobs, evaluation features, or drivers
- `cleanup/*`: refactoring, removal of obsolete scripts/tests, documentation reorganization

A change is considered “ready for main” once:
- unit/integration tests pass (`pytest`)
- the end-to-end smoke run is successful (when the external toolchain is available)
- the change does not break the report artifact contracts (`*.pipeline.json` and `aggregate.csv` schemas)

### 7.3 Change control and traceability

To ensure results are traceable to code versions:

- Each experiment run is associated with:
  - a commit hash (recorded manually in experiment logs, or optionally stored in `pipeline.json.params.git_commit`)
  - a deterministic seed
  - a complete parameter record in `*.pipeline.json`

- All plots and tables included in the report must be derivable from:
  - the `*.pipeline.json` files produced by the runner, aggregated into
  - `aggregate.csv` via `benchmark.aggregate_pipeline_reports`

This design avoids dependence on ad-hoc terminal logs and makes the report reproducible.

### 7.4 Artifact naming and storage policy

Artifacts are generated using prefix-based naming to guarantee consistency and simplify cleanup.

For a run with prefix `P`:
- `P.ref.fasta`, `P.truth.vcf.gz`, `P.reads.fastq`, `P.bam`, `P.called.vcf.gz`
- WhatsHap outputs: `P.ws*.phased.vcf`, `P.ws*.summary.json`, `P.ws*.eval.json`
- Run report: `P.pipeline.json`

Experiment directories are structured as:
- `<outdir>/<experiment_section>/...` containing per-seed runs and final `aggregate.csv` + `plots/`.

Artifacts are treated as:
- **source-of-truth**: `*.pipeline.json` (and `*.eval.json`, `*.summary.json`)
- **derived**: `aggregate.csv`, plots

### 7.5 Baseline and release management

Key “release” points in the project correspond to:
- stable versions used for report figures and tables
- tagged commits (optional) such as `v1-report-baseline`

A recommended approach:
- tag the commit used to generate final report experiments
- store experiment outputs in a dedicated `output/<tag>/` directory or archive them as a zip for submission

### 7.6 Backup and recovery

To avoid data loss:
- code and docs are backed up via the remote Git repository
- final experiment outputs used in the report should be archived (zip) and stored separately from transient debugging outputs