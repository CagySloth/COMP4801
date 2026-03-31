## 8. Software Quality Assurance Plan (SQAP-style)

This section describes the quality assurance activities used to ensure that the platform is correct, maintainable, and suitable for generating report-grade experimental results.

### 8.1 Quality goals

The system targets the following quality properties:

- **Correctness:** Metrics, phasing outputs, and run reports reflect the true pipeline behavior.
- **Reproducibility:** Runs are seed-controlled and configurations are recorded in `*.pipeline.json`.
- **Traceability:** Report plots/tables are derived from stored machine-readable artifacts.
- **Maintainability:** The codebase is modular (stage-based) and supports incremental addition of realism knobs and experiments.
- **Robustness:** Failures in external tools are detected early and surfaced clearly.

### 8.2 Standards and conventions

- Code follows a consistent CLI-driven structure using `python -m <module>` entrypoints.
- Outputs follow prefix-based naming conventions, enabling deterministic file discovery.
- Module-level documentation is maintained under `docs/` and updated when new knobs/drivers are introduced.

### 8.3 Verification activities

The following verification activities are used:

1. **Unit and integration tests (`pytest`)**
   - Validate metric computation, adapter correctness, and CLI routing.
   - Provide regression coverage for core evaluation logic.

2. **System smoke validation**
   - Validates the full long-read workflow through:
     - `benchmark.longread_pipeline_runner` (single-run)
     - `benchmark.experiment_driver` (multi-run)
   - Confirms that expected artifacts exist and that JSON/CSV contracts are satisfied.

3. **Realism knob validation**
   - Each realism knob is validated using:
     - metadata integrity checks (`*.ref.meta.json`, `*.reads.meta.json`, `*.truth.meta.json`)
     - sanity expectations on metrics (e.g., duplication tends to reduce calling recall and/or increase fragmentation)

### 8.4 Defect handling and regression prevention

- Defects are addressed by:
  - reproducing with a fixed seed and saving the failing run’s `*.pipeline.json`
  - adding a regression test when feasible (unit test or small integration test)
- A change is not considered complete until it:
  - passes the relevant tests, and
  - preserves the reporting contract (pipeline JSON schema and aggregation columns)

### 8.5 Documentation QA

Documentation is kept consistent with the codebase by:
- referencing runnable CLIs (`python -m ...`) rather than informal scripts
- documenting new flags when new realism knobs are added
- ensuring docs describe both oracle and called evaluation regimes

### 8.6 Limitations

- Full system validation requires external binaries (`minimap2/samtools/bcftools`), so CI-style testing may only cover unit/integration tests.
- Some realism properties are approximations; validation focuses on correctness of implementation and expected qualitative effects rather than perfect biological fidelity.