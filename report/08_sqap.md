## 8. Software Quality Assurance Plan (SQAP-style)

Quality assurance activites used to ensure the correctness and maintability of the platform will be covered in this section.

### 8.1 Quality assurance goals

The following quality properties are targeted by the quality assurance plan:

- **Correctness:** Evaluated metrics, phased outputs, and reports represent actual pipeline behavior accurately.
- **Reproducibility:** Seed and parameter configurations of experiment runs are recorded in `*.pipeline.json`.
- **Traceability:** Plots and tables in the report are generated from archived machine-readable data.
- **Maintainability:** The codebase is modularized and incremental implementation of realism knobs and experiments are supported.
- **Robustness:** Failures in external tools are detected and presented clearly.

### 8.2 Standards and conventions

- A consistent CLI-driven structure with `python -m <module>` entrypoints are followed.
- A prefix-based naming system is followed in output generation, ensuring deterministic file searching.
- Module-level documentation, under `docs/`, is maintained and updated regularly.

### 8.3 Verification activities

The following verification activities are used:

1. **Unit and integration tests using `pytest`**
   - Validation of metric computation, adapter correctness, and CLI routing.
   - Regression coverage for core evaluation logic is provided.

2. **System smoke validation**
   - Validates the full long-read workflow through:
     - `benchmark.longread_pipeline_runner` (single-run).
     - `benchmark.experiment_driver` (multi-run).
   - Existence of expected files and JSON/CSV contracts are verified.

3. **Realism knob validation**
   - Validation of each realism knob is performed through:
     - Integrity check for metadata (`*.ref.meta.json`, `*.reads.meta.json`, `*.truth.meta.json`).
     - Sanity check on output metrics.

### 8.4 Management of defects and prevention of regression

- Defects are addressed by:
  - Replicating with a fixed seed and saving the configuration `*.pipeline.json` of the failed run.
  - A regression test (unit test or small integration test) will be added when feasible.
- Changes and implementations are considered complete when:
  - Relevant tests are passed.
  - Reporting contract is preserved.

### 8.5 Documentation QA

Documentation is kept up-to-date with the codebase by:

- Runnable CLIs are referenced.
- New CLI flags are documented when new features are implemented.
- Ensuring documentation describe both oracle and called evaluation frameworks.

### 8.6 Limitations

- As full system validation would require validation of external tools, such as `minimap2`, `samtools`, and `bcftools`, only unit and integration tests are covered in CI-style testing.
- Some realism knobs are approximations of real genome sequencing characteristics. Validation of these realism features will focuses on correctness of implementation and expected effects, rather than biological truth.