## 3. Functional Requirements (FRS-style)

The functional requirements of the implemented benchamrking platform are covered in this section.

### 3.1 Perspective and scope

The system serves as a benchmarking platform for analysing WhatsHap's phasing capabilities in an ONT-like long-read phasing workflow. The pipeline follows this practical sequencing and phasing workflow:

**Reference FASTA → reads FASTQ → alignment BAM → variants VCF → phased VCF → evaluation reports**

The platform supports two evaluation frameworks:

- **Oracle (ground-truth) VCF:** isolates phasing performance by using the truth variant set.
- **Called VCF:** measures end-to-end performance including variant-calling limitations.

### 3.2 Users and operating environment

**Primary users:** developers or reserchers conducting controlled experiments for WhatsHap and generating plots and tables for analysis.

**Operating environment constraints:**

- Python modules are executed via `python -m ...` entrypoints.
- The end-to-end pipeline requires external tools on `PATH`:
  - `minimap2`
  - `samtools`
  - `bcftools`
- The pipeline runner checks tool availability at runtime and fails early if dependencies are missing.

### 3.3 Naming conventions and artifacts (prefix-based)

Most pipeline outputs are derived from a single `--prefix` argument used by `benchmark.longread_pipeline_runner`. Given `--prefix output/exp1/runA`, the system produces a consistent family of artifacts, including:

- Reference artifacts: `*.ref.fasta`, `*.ref.meta.json`.
- Truth artifacts: `*.truth.vcf`, `*.truth.vcf.gz`, `*.truth.meta.json`.
- Oracle VCF artifacts: `*.oracle.vcf`, `*.oracle.vcf.gz`.
- Simulated reads: `*.reads.fastq`, `*.reads.truth.tsv`.
- Alignment artifacts: `*.bam`, `*.bam.bai`.
- Called variants: `*.called.vcf.gz`, `*.called.vcf.gz.tbi`.
- WhatsHap outputs: `*.ws*.phased.vcf`, `*.ws*.eval.json`, `*.ws*.summary.json`.
- Run-level report: `*.pipeline.json`.

Aggregated experiment directories additionally produce:

- `aggregate.csv`
- plots under `plots/`

This prefix-based scheme is a functional requirement because it allows single-run reproducibility, experiment aggregation, and straightforward discovery of related artifacts.

### 3.4 Functional requirements

#### FR-1 Reference generation
The system shall generate a synthetic reference genome FASTA with configurable realism presets and optional duplicated segments.

**Implementation mapping:** implemented by `dataset.longread.reference`, invoked through `benchmark.long` `read_pipeline_runner`, producing `*.ref.fasta` and `*.ref.meta.json`.

#### FR-2 Ground-truth variant and haplotype generation
The system shall generate diploid ground-truth variants and two haplotype sequences, and export both a truth VCF and haplotype FASTAs.

**Implementation mapping:** implemented by `dataset.longread.truth`, invoked through `benchmark.longread_` `pipeline_runner`, producing `*.truth.vcf`, `*.truth.vcf.gz`, `*.hap1.fasta`, `*.hap2.fasta`, `*.truth.meta.json`, and `*.oracle.vcf.gz`.

#### FR-3 ONT-like read simulation
The system shall simulate long reads from the diploid haplotypes, producing FASTQ with quality strings and read-truth metadata.

**Implementation mapping:** implemented by `dataset.longread.readsim`, invoked through `benchmark.longrea` `d_pipeline_runner`, producing `*.reads.fastq` and `*.reads.truth.tsv`, with support for configurable read-length, coverage, and burst-error models.

#### FR-4 Read alignment
The system shall align simulated reads to the reference and output a sorted, indexed BAM suitable for downstream variant calling and WhatsHap phasing.

**Implementation mapping:** performed by `benchmark.longread_pipeline_runner` using `minimap2` and `samtools`, producing `*.bam` and `*.bam.bai`.

#### FR-5 Variant calling
The system shall call variants from the aligned BAM against the reference and produce a called VCF.

**Implementation mapping:** performed by `benchmark.longread_pipeline_runner` using `bcftools mpileup | ` `bcftools call`, producing `*.called.vcf.gz` and its index.

#### FR-6 WhatsHap phasing (vendored integration)
The system shall phase variants using WhatsHap with BAM+VCF inputs and produce phased VCF output, using the vendored WhatsHap core for controlled benchmarking.

**Implementation mapping:** exposed through `algorithms.cli.phase diploid-whats-bam`, implemented by `algorithms.diploid.whatshap_bam_driver`, and invoked through `benchmark.longread_pipeline_runner`, producing phased VCFs and run summaries.

#### FR-7 Phasing evaluation against truth
The system shall evaluate phased output against the truth VCF and output evaluation metrics focused on phasing correctness and fragmentation.

**Implementation mapping:** implemented by `benchmark.vcf_phase_eval`, invoked through `benchmark.longread_pipeline_runner`, producing `*.eval.json` with metrics such as phase accuracy, switch error rate, phase-set count, and shared-site counts.

**SNP-only policy for indels:** when indels are enabled, the system shall support `--phase-snps-only` and `--eval-snps-only` so that SNP phasing evaluation remains valid despite possible indel-representation mismatches.

#### FR-8 End-to-end reporting and aggregation
The system shall produce a single machine-readable report per run and provide aggregation utilities for experiment directories.

**Implementation mapping:** `benchmark.longread_pipeline_runner` produces `*.pipeline.json`; `benchmark.` `aggregate_pipeline_reports` converts many run reports into `aggregate.csv` for downstream analysis and plotting.

The pipeline report shall record:

- Parameter values used for the run.
- File paths for produced artifacts.
- Calling precision and recall when called-mode evaluation is enabled.
- Phasing evaluation metrics for oracle and/or called runs.

#### FR-9 Experiment driver for systematic studies
The system shall provide a driver to execute predefined experiment suites across multiple seeds, aggregate results, and generate plots.

**Implementation mapping:** implemented by `benchmark.experiment_driver`, which runs seeded experiment suites and produces per-experiment `aggregate.csv` files and report-ready plots.

### 3.5 Non-functional requirements

- **Reproducibility:** runs must be seed-controlled and record parameters in `pipeline.json`.
- **Traceability:** all plots and tables must be derived from stored machine-readable outputs (`pipeline.json`, `aggregate.csv`).
- **Maintainability:** core evaluation and aggregation logic must be testable and separated from legacy or exploratory scripts.
- **Portability:** code must run on a standard Linux environment with documented external tool dependencies.