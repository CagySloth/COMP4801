## 4. System Design (SDD-style)

The design of the WhatsHap benchmarking platform, including architecture, component decomposition, interfaces, design decisions, and extensibility, are covered in this section.

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
   Expose configurable stressors such as duplicated regions, coverage dropout, burst errors, and indels so they can be varied independently in controlled experiments.

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
   Compare phased VCF output against truth VCF and measure performance metrics.

8. **Experiment orchestration**
   Conduct controlled experiments that systematically sweeps across realism settings and seed, then build reports for each benchmarking run and generate summary tables and plots.

#### 4.2.2 Codebase mapping

The pipeline is implemented through three primary areas:

- `dataset/longread/`
  Generation of reference genome and truth files, and simulation of ONT-like long-read sequencing.

- `algorithms/`
  Command line interface for phasing and phasing drivers for WhatsHap.

- `benchmark/`
  End-to-end organization, evaluation, summarization, and experiment control.

This separation of the pipeline maintain independency and modularity of data generation, phasing, and benchmarking, while pipeline runners are still able to execute the full pipeline.

### 4.3 Component decomposition

#### 4.3.1 `dataset.longread.reference`
Generates synthetic reference FASTA metadata where optional complexity presets and duplicated-region indications can be included for simulating realism.

#### 4.3.2 `dataset.longread.truth`
Generates truth variants, haplotype FASTAs and truth metadata for performance evaluation, and oracle VCFs that allows isolation of phasing performance from variant-caller.

#### 4.3.3 `dataset.longread.readsim`
Simulates ONT-like long-read sequencing with fully configurable realism knobs, including length distribution, coverage dropout, and correlated error bursts.

#### 4.3.4 `benchmark.longread_pipeline_runner`
Single-run end-to-end pipeline runner. Validates reqruired external dependencies, executes the full pipeline from data generation to phasing performance evaluation and create a summarised `*.pipeline.json` report.

#### 4.3.5 `algorithms.cli.phase`
Provides a unified CLI entrypoint for phasing backends so that phasing can be invoked in a consistent way from the runner.

#### 4.3.6 `algorithms.diploid.whatshap_bam_driver`
Implements WhatsHap-based phasing from BAM + VCF using the vendored core. It performs readset construction, read selection, solving, PS assignment, and phased VCF writing. Reads covering fewer than two variants are excluded to satisfy WhatsHap read-selection requirements.

#### 4.3.7 `benchmark.vcf_phase_eval`
Compares phased VCF output against truth and computes correctness, completeness, and fragmentation metrics.

#### 4.3.8 `benchmark.aggregate_pipeline_reports`
Converts many per-run `*.pipeline.json` files into a single `aggregate.csv`, which acts as the standard interface between experiment execution and downstream plotting.

#### 4.3.9 `benchmark.experiment_driver`
Runs predefined experiment suites across multiple seeds and parameter settings, then aggregates results and generates report-ready plots.

### 4.4 Data contracts and file-level interfaces

The design relies on stable file contracts between stages. The main exchanged artifact types are:

- **FASTA** for reference and haplotypes.
  (`*.ref.fasta`, `*.hap1.fasta`, `*.hap2.fasta`)

- **FASTQ** for simulated reads.
  (`*.reads.fastq`)

- **BAM** for sorted and indexed alignments.
  (`*.bam`, `*.bam.bai`)

- **VCF** for truth, oracle, called, and phased variant sets.
  (`*.truth.vcf.gz`, `*.oracle.vcf.gz`, `*.called.vcf.gz`, `*.ws*.phased.vcf`)

- **JSON** for metadata and machine-readable reports.
  (`*.ref.meta.json`, `*.truth.meta.json`, `*.summary.json`, `*.eval.json`, `*.pipeline.json`)

Prefix-based naming ensures that all artifacts from a run can be discovered from a single run identifier, which simplifies orchestration, aggregation, and cleanup.

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
  - Lower `call_recall`.
  - Possible precision changes depending on thresholding and error severity.

- **Phaser-focused effects**
  - More phase sets under fragmentation-heavy stressors.
  - Possible switch-error increases when fewer reliable allele observations survive.

- **End-to-end effects**
  - The strongest performance degradation is in effective phased recall, which depends on both variant overlap and phasing completeness.

Realism knobs will be constructed as independently variable and testable stressors to validate these assumptions, before combining them for interaction and optimization studies.

### 4.6 Key design decisions

1. **Prefix-based artifact naming**
   To simplify discovery, cleanup and aggreation, a single prefix determines the artifact family of each run.

2. **Oracle vs called framework**
   To enable attribution and breakdown of phasing performance loss from variant calling and phasing, the platform separates phaser-limited behavior from end-to-end behavior.

3. **Vendored WhatsHap core**
   To maintain full control and understanding over benchmarking behavior and avoid dependency drift, the core algorithms of WhatsHap are vendored instead of imported.

4. **SNP-only phasing and evaluation when indels are present**
   To better evaluate the performance on SNP-phasing, the pipeline prevents indel representation differences from dominating phasing evaluation.

5. **Machine-readable single-run report as the source of truth**
   To ensure traceability and reproducibility, aggregation and plotting are driven from `*.pipeline.json` summaries rather than logs.

### 4.7 Extensibility

Two main incremental extension directions are supported by the platform design.

- **New realism knobs** can be implemented in the any data-generation stage, accessible through the CLI and pipeline runner.
- **Alternative phasers** can be integrated and evaluated by building drivers or adaptors under the `algorithms/` module and registering new subcommands in `algorithms.cli.phase`.

The extensibility of this platform makes the platform capable of future evaluation and analsysis beyond WhatsHap, with additional realism knobs and evaluation policies.