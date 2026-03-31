## 3. Functional Requirements (FRS-style)

This section specifies the functional requirements for the implemented benchmarking platform.

### 3.1 Product perspective and scope

The system is a reproducible benchmarking platform for evaluating WhatsHap in an ONT-like long-read phasing workflow. The pipeline targets the practical path:

**Reference FASTA → reads FASTQ → alignment BAM → variants VCF → phased VCF → evaluation reports**

The platform supports two evaluation regimes:

- **Oracle (ground-truth) VCF**: isolates phasing performance given the correct variant set.
- **Called VCF**: measures end-to-end performance including variant calling limitations.

### 3.2 Users and operating environment

**Primary user:** a developer/research user running controlled experiments and producing plots/tables for analysis and reporting.

**Operating environment constraints:**
- Python modules are executed via `python -m ...` entrypoints.
- The end-to-end pipeline requires external tools available on `PATH`:
  - `minimap2`, `samtools`, `bcftools`
- The pipeline runner checks tool availability at runtime (and fails early if missing).

### 3.3 Naming conventions and artifacts (prefix-based)

Most pipeline outputs are derived from a single `--prefix` argument used consistently by `benchmark.longread_pipeline_runner`. Given `--prefix output/exp1/runA`, the system produces artifacts following this convention:

- Reference: `output/exp1/runA.ref.fasta`, `output/exp1/runA.ref.meta.json`
- Truth: `output/exp1/runA.truth.vcf`, `output/exp1/runA.truth.vcf.gz` (+ index), `output/exp1/runA.truth.meta.json`
- Oracle VCF: `output/exp1/runA.oracle.vcf`, `output/exp1/runA.oracle.vcf.gz` (+ index)
- Reads: `output/exp1/runA.reads.fastq`, `output/exp1/runA.reads.truth.tsv`, optional read meta JSON
- Alignment: `output/exp1/runA.bam` (+ `.bai`)
- Called variants: `output/exp1/runA.called.vcf.gz` (+ `.tbi`)
- WhatsHap outputs (called/oracle runs): `*.ws*.phased.vcf`, `*.ws*.eval.json`, `*.ws*.summary.json`
- Pipeline report: `output/exp1/runA.pipeline.json`

Aggregated experiments additionally produce:
- `aggregate.csv` under each experiment directory
- plots under `plots/` within the experiment directory

### 3.4 Functional requirements

#### FR-1 Reference generation
The system shall generate a synthetic reference genome FASTA with configurable realism presets and optional duplicated segments.

**Implementation mapping:**
- CLI: `python -m dataset.longread.reference`
- Called by: `python -m benchmark.longread_pipeline_runner`
- Parameters include:
  - `--ref-preset {plain,toy,realistic}`
  - `--ref-length`
  - duplication knobs: `--dup-segments`, `--dup-len`, `--dup-min-gap`
- Outputs: `*.ref.fasta`, `*.ref.meta.json`

#### FR-2 Ground-truth variant and haplotype generation
The system shall generate diploid ground-truth variants and two haplotype sequences, and export both a truth VCF and haplotype FASTAs.

**Implementation mapping:**
- CLI: `python -m dataset.longread.truth`
- Called by: `python -m benchmark.longread_pipeline_runner`
- Parameters include:
  - SNP knobs: `--num-snps`, `--het-rate`
  - optional indel knobs: `--num-indels`, `--indel-min-len`, `--indel-max-len`, `--indel-het-rate`
- Outputs:
  - `*.truth.vcf` (+ `.gz` and index)
  - `*.hap1.fasta`, `*.hap2.fasta`
  - `*.truth.meta.json`
  - `*.oracle.vcf` (+ `.gz` and index)

#### FR-3 ONT-like read simulation
The system shall simulate long reads from the diploid haplotypes, producing FASTQ with quality strings and read-truth metadata.

**Implementation mapping:**
- CLI: `python -m dataset.longread.readsim`
- Called by: `python -m benchmark.longread_pipeline_runner`
- Parameters include:
  - `--num-reads`, `--min-len`, `--max-len`
  - `--platform {ont,perfect}`, `--ont-profile {classic,q20}`
  - length model knobs: `--len-model {uniform,lognormal}`, `--ln-mean`, `--ln-sigma`
  - start/coverage model knobs: `--start-model {uniform,dropout}`, `--dropout-fraction`, `--dropout-block-len`
  - burst error knobs (passed via pipeline runner): `--burst-prob`, `--burst-count`, `--burst-len`, `--burst-mult`
- Outputs:
  - `*.reads.fastq`
  - `*.reads.truth.tsv` (read origin / truth mapping)
  - optional `*.reads.meta.json` if enabled by the generator

#### FR-4 Read alignment
The system shall align simulated reads to the reference and output a sorted, indexed BAM suitable for downstream variant calling and WhatsHap phasing.

**Implementation mapping:**
- Invoked by: `python -m benchmark.longread_pipeline_runner`
- Tools: `minimap2`, `samtools sort`, `samtools index`
- Parameter: `--map-preset` (default `map-ont`)
- Outputs: `*.bam`, `*.bam.bai`

#### FR-5 Variant calling
The system shall call variants from the aligned BAM against the reference and produce a called VCF.

**Implementation mapping:**
- Invoked by: `python -m benchmark.longread_pipeline_runner`
- Tools: `bcftools mpileup | bcftools call`
- Parameters include:
  - `--call-min-mapq` (bcftools `-q`)
  - `--call-min-baseq` (bcftools `-Q`)
- Outputs: `*.called.vcf.gz`, `*.called.vcf.gz.tbi`

#### FR-6 WhatsHap phasing (vendored integration)
The system shall phase variants using WhatsHap with BAM+VCF inputs and produce phased VCF output, using the vendored WhatsHap core for controlled benchmarking.

**Implementation mapping:**
- CLI entrypoint: `python -m algorithms.cli.phase diploid-whats-bam`
- Driver: `algorithms/diploid/whatshap_bam_driver.py`
- Called by: `python -m benchmark.longread_pipeline_runner`
- Parameters include:
  - `--bam`, `--vcf`, `--output-vcf`, `--output-prefix`
  - `--max-coverage`, `--min-mapq`, `--min-baseq`
- Outputs:
  - `*.phased.vcf`
  - `*.summary.json` (instrumentation + counts)

#### FR-7 Phasing evaluation against truth
The system shall evaluate phased output against the truth VCF and output evaluation metrics focused on phasing correctness and fragmentation.

**Implementation mapping:**
- CLI: `python -m benchmark.vcf_phase_eval`
- Called by: `python -m benchmark.longread_pipeline_runner`
- Outputs: `*.eval.json` including:
  - `phase_accuracy_blockflip`
  - `switch_error_rate`
  - `num_phase_sets`
  - counts of shared/het/phased records

**SNP-only policy for indels:**
When indels are enabled (`--num-indels > 0`), the system supports filtering to biallelic SNPs for phasing and/or evaluation:
- Pipeline runner flags: `--phase-snps-only`, `--eval-snps-only`
- These ensure evaluations remain valid when indel representations differ across truth/called/predicted VCFs.

#### FR-8 End-to-end reporting and aggregation
The system shall produce a single machine-readable report per run and provide aggregation utilities for experiment directories.

**Implementation mapping:**
- Run-level report: `*.pipeline.json` produced by `benchmark.longread_pipeline_runner`
- Aggregation:
  - `python -m benchmark.aggregate_pipeline_reports_full --root <dir> --out <csv>`
  - Used by `benchmark.experiment_driver`

The pipeline report shall record:
- parameter values used for the run
- file paths for produced artifacts
- calling precision/recall (when `--vcf-source` includes called)
- phasing evaluation metrics for oracle and/or called runs

#### FR-9 Experiment driver for systematic studies
The system shall provide a driver to execute predefined experiment suites across multiple seeds, aggregate results, and generate plots.

**Implementation mapping:**
- CLI: `python -m benchmark.experiment_driver --outdir ... --seeds ...`
- Selective execution: `--only <section1,section2,...>`
- Outputs per experiment section:
  - per-seed run artifacts
  - `aggregate.csv` and `plots/` for report-ready figures

### 3.5 Non-functional requirements (supporting the FYP goals)

- **Reproducibility:** runs must be seed-controlled and record parameters in `pipeline.json`.
- **Traceability:** all plots/tables must be derived from stored machine-readable outputs (`pipeline.json`, `aggregate.csv`).
- **Maintainability:** core evaluation and aggregation logic must be unit tested; legacy scripts/tests should be isolated from the main pipeline.
- **Portability:** code must run on a standard Linux environment with documented external tool dependencies.