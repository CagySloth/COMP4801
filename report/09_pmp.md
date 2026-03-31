## 9. Project Management Plan (PMP-style)

This section summarizes the project’s deliverables, milestones, and risk management approach, focusing on the work completed toward the final report and experimental study.

### 9.1 Deliverables

The final deliverables of this project are:

1. **Benchmarking platform implementation**
   - End-to-end pipeline:
     - reference → truth/oracle → read simulation → alignment → calling → phasing → evaluation → report
   - Vendored WhatsHap integration for controlled benchmarking

2. **Realism knob implementations**
   - Reference-level stressors (e.g., duplications / complexity presets)
   - Read-level stressors (dropout, length models, correlated error bursts)
   - Truth-level stressors (indels) with SNP-only policy for meaningful evaluation

3. **Experiment infrastructure**
   - Experiment driver for running systematic studies across seeds and knobs
   - Aggregation and plotting utilities producing report-ready figures and CSVs

4. **Validation artifacts**
   - Unit/integration tests and smoke-run outputs
   - `*.pipeline.json` reports and aggregated `aggregate.csv` datasets used by the report

5. **Final report**
   - Documentation of requirements, design, implementation, testing, experiments, and conclusions

### 9.2 Work plan and milestones (final-phase focus)

A practical final-phase work plan:

- **M1: Platform stabilization**
  - ensure smoke validation passes on target environment
  - confirm `pipeline.json` and `aggregate.csv` contain required columns for report tables/plots

- **M2: Realism knob completion and validation**
  - validate each realism knob with targeted runs and sanity checks
  - ensure SNP-only policy is consistently applied when indels are enabled

- **M3: Experiment execution**
  - run the planned experiment suites (oracle vs called regimes)
  - aggregate results and generate plots

- **M4: Report completion**
  - integrate experiment plots/tables into the report
  - finalize discussion of optimization opportunities and limitations

### 9.3 Risks and mitigations

| Risk | Impact | Mitigation |
|---|---|---|
| External toolchain not available / inconsistent versions | prevents end-to-end experiments | document dependencies; fail early in runner; test on target machine |
| Evaluation artifacts due to indel representation mismatch | misleading accuracy/switch metrics | SNP-only phasing/evaluation flags; validate with indel sweeps |
| Time constraints in final phase | reduced experiment coverage | prioritize experiments that directly support report claims; use smaller parameter sweeps |
| Overfitting to synthetic data | reduced realism | use realism knobs to approximate failure modes; discuss threats to validity |
| Codebase complexity | harder to maintain | cleanup branch; isolate legacy scripts/tests; keep `pipeline.json` schema stable |

### 9.4 Tracking and reporting

Progress is tracked through:
- Git commits and branch history for code changes
- prefix-scoped outputs and `pipeline.json` reports for experiment evidence
- aggregated CSVs and plots for report integration