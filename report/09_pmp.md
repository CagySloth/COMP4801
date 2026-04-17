## 9. Project Management Plan (PMP-style)

This project's deliverables, milestones, and risk management approach, will be covered in this section.

### 9.1 Deliverables

The expected final deliverables of this project includes:

1. **Implementation of the benchmarking platform**
   - End-to-end pipeline:
     - Reference → Truth/oracle → Read simulation → Alignment → Variant-calling → Phasing → Evaluation → Visualization
   - Vendored WhatsHap integration for controlled benchmarking and experimenting.

2. **Realism knob implementations**
   - Reference-level stressors, such as duplications and complexity presets.
   - Read-level stressors, such as dropout, length models, and correlated error bursts.
   - Truth-level stressors, such as indel insertion, with SNP-only policy for meaningful evaluation.

3. **Experiment infrastructure**
   - Experiment driver to facilitate systematic investigation on phasing performance across seeds, setups, and realism configurations.
   - Aggregation and visualization tools to produce report-ready figures and tables.

4. **Validation artifacts**
   - Outputs from unit tests, integration tests, and end-to-end smoke-runs.
   - `*.pipeline.json` reports and aggregated `aggregate.csv` datasets used by this report.

5. **Final report**
   - A detailed documentation of requirements, design, implementation, testing, experiments, findings, and conclusions for this project.

### 9.2 Work plan and milestones

A practical work plan for the second-phase (semester 2):

- **M1: Platform stabilization**
  - Ensure that the platform operates correctly on passes on target environment through smoke runs.
  - Ensure that output from the platform contain required information for report plots and tables.

- **M2: Realism knob implementation and validation**
  - Validation of each realism knob through targeted tests and sanity checks.
  - Confirm that the SNP-only policy is applied when indels are enabled.

- **M3: Experiment execution**
  - Execute the designed experiment suites.
  - Aggregate and visualize results, and explore optimization opportunity.

- **M4: Report completion**
  - Incorporation of experiment plots and tables into the final report.
  - Finalization of findings on optimization opportunities and limitations.

### 9.3 Risks and mitigations

| Risk | Impact | Mitigation |
|---|---|---|
| Inconsistent availability and versioning of external tools | Inability and inconsistency in end-to-end experiments | Document external tools dependencies; Validation on targeted environment |
| Mismatch in indel representations leading to inaccurate evaluation | Misleading metrics | SNP-only phasing and evaluation flags; Validation with indel sweep experiment runs |
| Time constraints in final phase | Reduced experiment coverage; Immature optimization findings | Design experimentation plan; Prioritize critical experiments; Reduce parameter sweep size |
| Overfitting to synthetic data | Reduced realism | Approximate failure modes through realism knobs; Ensure randomness in experiments |
| Codebase complexity | Harder to maintain | Regular cleanup and refactoring; Isolate legacy scripts; Maintain consistency in `pipeline.json` schema |

### 9.4 Tracking and reporting

The research process can be tracked through:

- Git commits and branch history documenting code modifications and feature implementations.
- Experiment outputs and `pipeline.json` reports for research traceability.
- Aggregated CSVs and plots for report integration.