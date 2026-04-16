## 4. System Design (SDD-style)

This section describes the design of the WhatsHap benchmarking platform, focusing on the architecture, component responsibilities, interfaces, and design decisions that support reproducible end-to-end long-read phasing experiments.

### 4.1 Design objectives

The system is designed to satisfy five main objectives.

1. **Benchmark WhatsHap in practical workflows**  
   Evaluate phasing in an end-to-end long-read pipeline (`FASTA → FASTQ → BAM → VCF → phased VCF`), not only on idealized matrices.

2. **Separate error sources**  
   Support both oracle and called variant regimes so that phasing-limited behaviour can be distinguished from end-to-end caller-limited behaviour.

3. **Be reproducible and traceable**  
   Ensure that every run is seed-controlled and produces a machine-readable report (`*.pipeline.json`) recording parameters, file paths, and metrics.

4. **Be modular and scriptable**  
   Implement each pipeline stage as a dedicated `python -m ...` entrypoint so stages can be run independently or orchestrated through runners and experiment drivers.

5. **Support realism knobs and systematic sweeps**  
   Expose configurable stressors such as duplicated regions, coverage dropout, burst errors, and indels so they can be varied independently in controlled experiments. :contentReference[oaicite:2]{index=2}

### 4.2 High-level pipeline structure

The system is organized as a sequential pipeline with clearly defined stage boundaries. Each stage has a single responsibility and produces artifacts that become inputs to the next stage.

#### 4.2.1 Conceptual pipeline stages

1. **Reference synthesis**  
   Generate a synthetic reference FASTA with optional complexity presets and duplicated segments.

2. **Truth synthesis**  
   Generate phased truth variants, haplotypes, and an oracle VCF.

3. **Read simulation**  
   Simulate ONT-like reads with configurable length, coverage, and error models.

4. **Alignment and preprocessing**  
   Align reads to the reference and produce a sorted, indexed BAM.

5. **Variant calling**  
   Produce a called VCF from the aligned BAM.

6. **Phasing**  
   Phase oracle and/or called variants using the vendored WhatsHap core.

7. **Evaluation and reporting**  
   Compare phased output to truth and record metrics and provenance.

8. **Experiment orchestration**  
   Run controlled sweeps across seeds and realism knobs, aggregate per-run reports, and generate CSV summaries and plots. 

#### 4.2.2 Codebase mapping

The conceptual pipeline is implemented through three main repository areas:

- `dataset/longread/`  
  Reference generation, truth generation, and ONT-like read simulation.

- `algorithms/`  
  Unified phasing CLI and WhatsHap-based phasing drivers.

- `benchmark/`  
  End-to-end orchestration, evaluation, aggregation, and experiment control. 

This separation keeps data generation, phasing, and benchmarking concerns distinct while still allowing the full workflow to be executed through a single pipeline runner.

### 4.3 Component decomposition and responsibilities

The main components and their roles are as follows.

#### 4.3.1 `dataset.longread.reference`
Generates the synthetic reference FASTA and reference metadata, including optional complexity presets and duplicated-region annotations.

#### 4.3.2 `dataset.longread.truth`
Generates truth variants, haplotype FASTAs, truth metadata, and the oracle VCF used to isolate phasing behaviour.

#### 4.3.3 `dataset.longread.readsim`
Simulates ONT-like reads and supports configurable realism knobs including length distribution, coverage dropout, and correlated error bursts.

#### 4.3.4 `benchmark.longread_pipeline_runner`
Acts as the canonical single-run orchestrator. It validates external dependencies, executes the full workflow from reference generation through evaluation, and writes the final `*.pipeline.json` report.

#### 4.3.5 `algorithms.cli.phase`
Provides a unified CLI entrypoint for phasing backends so that phasing can be invoked in a consistent way from the runner.

#### 4.3.6 `algorithms.diploid.whatshap_bam_driver`
Implements WhatsHap-based phasing from BAM + VCF using the vendored core. It performs readset construction, read selection, solving, PS assignment, and phased VCF writing. Reads covering fewer than two variants are excluded to satisfy WhatsHap read-selection requirements. :contentReference[oaicite:5]{index=5}

#### 4.3.7 `benchmark.vcf_phase_eval`
Compares phased VCF output against truth and computes correctness, completeness, and fragmentation metrics.

#### 4.3.8 `benchmark.aggregate_pipeline_reports`
Converts many per-run `*.pipeline.json` files into a single `aggregate.csv`, which acts as the standard interface between experiment execution and downstream plotting.

#### 4.3.9 `benchmark.experiment_driver`
Runs predefined experiment suites across multiple seeds and parameter settings, then aggregates results and generates report-ready plots. :contentReference[oaicite:6]{index=6}

### 4.4 Data contracts and file-level interfaces

The design relies on stable file contracts between stages. The main exchanged artifact types are:

- **FASTA** for reference and haplotypes  
  (`*.ref.fasta`, `*.hap1.fasta`, `*.hap2.fasta`)

- **FASTQ** for simulated reads  
  (`*.reads.fastq`)

- **BAM** for sorted and indexed alignments  
  (`*.bam`, `*.bam.bai`)

- **VCF** for truth, oracle, called, and phased variant sets  
  (`*.truth.vcf.gz`, `*.oracle.vcf.gz`, `*.called.vcf.gz`, `*.ws*.phased.vcf`)

- **JSON** for metadata and machine-readable reports  
  (`*.ref.meta.json`, `*.truth.meta.json`, `*.summary.json`, `*.eval.json`, `*.pipeline.json`) :contentReference[oaicite:7]{index=7}

Prefix-based naming ensures that all artifacts from a run can be discovered from a single run identifier, which simplifies orchestration, aggregation, and cleanup. :contentReference[oaicite:8]{index=8}

### 4.5 Realism knobs as configurable stressors

The platform treats realism knobs as controlled stressors. Each knob is designed to approximate a practical difficulty mode and to be measurable through changes in calling metrics, phasing metrics, or both.

#### 4.5.1 Stressor catalogue

The main stressors are:

- **Reference-level stressors**
  - reference complexity preset
  - duplicated segments

- **Truth-level stressors**
  - optional indels
  - truth-region constraints

- **Read-level stressors**
  - read length distribution
  - coverage dropout model
  - correlated error bursts

- **Runner / evaluation policy stressors**
  - oracle vs called VCF selection
  - SNP-only phasing / SNP-only evaluation when indels are present 

These knobs are intentionally parameterized, stage-local, and fully recorded in `*.pipeline.json`, so that experiment results remain reproducible and traceable.

#### 4.5.2 Measurement expectations

Increasing stressor severity is expected to affect the pipeline in different ways:

- **Caller-focused effects**
  - lower `call_recall`
  - possible precision changes depending on thresholding and error severity

- **Phaser-focused effects**
  - more phase sets under fragmentation-heavy stressors
  - possible switch-error increases when fewer reliable allele observations survive

- **End-to-end effects**
  - the strongest visible degradation is usually in effective phased recall, because it depends on both variant overlap and phasing completeness :contentReference[oaicite:10]{index=10}

These expectations justify treating realism knobs as independently testable stressors before combining them in interaction and optimization studies.

### 4.6 Key design decisions

Several design decisions are central to the system.

1. **Prefix-based artifact naming**  
   A single prefix determines the artifact family of one run, simplifying discovery, cleanup, and aggregation.

2. **Oracle vs called regimes**  
   The platform explicitly separates phaser-limited from end-to-end behaviour, enabling attribution of performance loss to calling vs phasing.

3. **Vendored WhatsHap core**  
   The WhatsHap core is vendored to avoid dependency drift and to keep benchmarking behaviour controlled and reproducible.

4. **SNP-only phasing/evaluation when indels are present**  
   This prevents indel representation differences from dominating SNP phasing evaluation.

5. **Machine-readable single-run report as the source of truth**  
   Aggregation and plotting are driven from `*.pipeline.json`, not from logs, which improves traceability and reproducibility. :contentReference[oaicite:11]{index=11}

### 4.7 Extensibility

The design supports incremental extension in two main directions.

- **New realism knobs** can be added in the relevant data-generation stage, surfaced through the pipeline runner, and then exposed to the experiment driver and aggregation layer.

- **Alternative phasers** can be integrated by adding a new driver under `algorithms/`, registering a new subcommand in `algorithms.cli.phase`, and extending the pipeline runner to collect the corresponding outputs and metrics. :contentReference[oaicite:12]{index=12}

This makes the platform suitable not only for the present WhatsHap study, but also for future controlled comparisons across stress conditions, phasing backends, and evaluation policies.