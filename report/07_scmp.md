## 7. Software Configuration Management Plan (SCMP-style)

The management of source code, experiment configurations, and generated results will be covered in this section.

### 7.1 Configuration items (CIs)

The primary configuration items are:

1. **Source code repository**
   - All Python modules and scripts under `dataset/`, `algorithms/`, and `benchmark/`
   - Vendored WhatsHap phasing core algorithm under `vendor/whatshap_core/`
   - Documentation under `docs/`

2. **Experiment configurations**
   - Experiment sets and parameter grids are specified in `benchmark/experiment_driver.py`
   - Single-run parameterization are specified by `benchmark/longread_pipeline_runner.py`

3. **Run artifacts and reports**
   - Prefix-scoped run artifacts (FASTA/FASTQ/BAM/VCF)
   - Machine-readable run reports: `*.pipeline.json`
   - Aggregated datasets: `aggregate.csv`
   - Plots derived from aggregates: `plots/*.png`

### 7.2 Version control strategy

Git is used to track changes to all code and documentation.

Roles of different branches:

- `main`: Stable and validated baseline for executing experiments and generating final results.
- `feature/*`: Isolated branches for development and validation of newly implemented features, including realism knobs, evaluation features, or phasing drivers.
- `cleanup/*`: Intermediate branches for refactoring of the codebase, removal of obsolete components, and documentation reorganization.

A change is considered to be ready for merging into the "main" branch once:

- Unit tests and integration tests are passed.
- End-to-end smoke run is successful and correct.
- Report artifact contracts are not affected.

### 7.3 Change traceability and controll

To ensure results are traceable to code versions:

- Experiment runs are associated with:
  - A commit hash.
  - A set of predetermined seeds.
  - A comprehensive record of parameter configuration in `*.pipeline.json`.

- All plots and tables from experiment runs that are included in the report must be derivable and reproducible from:
  - Parameter configuration `*.pipeline.json` files.

Dependency on terminal logs can be eliminated with this approach, maximizing report reproducibility.

### 7.4 Artifact naming and storage policy

To ensure consistency and facilitate cleanups, artifacts are generated with a well-defined prefix-based naming system.

For a run with prefix `P`:

- `P.ref.fasta`, `P.truth.vcf.gz`, `P.reads.fastq`, `P.bam`, `P.called.vcf.gz`.
- WhatsHap outputs: `P.ws*.phased.vcf`, `P.ws*.summary.json`, `P.ws*.eval.json`.
- Run report: `P.pipeline.json`.

Experiment directories are structured as:

- `output/<experiment_section>/...` containing per-seed runs and final `aggregate.csv` + `plots/`.

Artifacts are categorized as:

- **source-of-truth**: `*.pipeline.json`, `*.eval.json`, and `*.summary.json`.
- **derived**: `aggregate.csv`, visualizations.

### 7.5 Baseline and milestone management

Key milestones correspond to:

- Stable versions used for experiments, report figures and tables.
- Tagged commits.

### 7.6 Backup and recovery

To avoid data loss:

- Backups of scripts and documentations are created through remote Git repositories.
- Archives of experiment outputs used in the report are created and stored separately from the repositories.