## 5. Implementation

This section describes how the system is implemented, focusing on the repository structure, key modules, the end-to-end runner, and the WhatsHap integration choices that enable controlled benchmarking.

### 5.1 Repository structure (implementation map)

The repository is organized around three concerns: data generation, phasing algorithms, and benchmarking/orchestration.

- `dataset/longread/`  
  Implements long-read pipeline data generation:
  - `reference.py`: reference FASTA synthesis + complexity presets + duplication metadata
  - `truth.py`: truth variants + haplotypes (truth VCF + hap1/hap2 FASTA; optional indels)
  - `readsim.py`: ONT-like read simulation (FASTQ + read-truth TSV + meta JSON; length/dropout/bursts)

- `algorithms/`  
  Implements phasing algorithm CLIs and drivers:
  - `cli/phase.py`: unified entrypoint `python -m algorithms.cli.phase`
  - `diploid/whatshap_bam_driver.py`: WhatsHap phasing from BAM+VCF using vendored core
  - `vendor/whatshap_vendor.py` (or equivalent): ensures imports resolve to vendored WhatsHap

- `benchmark/`  
  Implements end-to-end orchestration, evaluation, experiments, and aggregation:
  - `longread_pipeline_runner.py`: single-run orchestrator producing `*.pipeline.json`
  - `vcf_phase_eval.py`: phasing evaluation vs truth VCF
  - `aggregate_pipeline_reports.py`: produces `aggregate.csv` from many `pipeline.json`
  - `experiment_driver.py`: runs experiment suites across seeds, aggregates and plots

- `tests/`  
  Unit tests for evaluation/aggregation and specific realism implementations.

- `vendor/whatshap_core/`  
  Vendored WhatsHap core module (`whatshap`) with compiled extensions, used for controlled benchmarking.

This structure keeps the “single source of truth” for pipeline behavior in `benchmark.longread_pipeline_runner` and keeps experiments as wrappers around it.

---

### 5.2 End-to-end orchestration: `benchmark.longread_pipeline_runner`

The canonical workflow is implemented as a single-run pipeline runner:

- Entry point: `python -m benchmark.longread_pipeline_runner`
- Primary argument: `--prefix <path-prefix>`

The runner enforces two design contracts:
1. **Prefix-based artifact naming**: all produced files share the same prefix, enabling easy discovery.
2. **Machine-readable run report**: each run produces a `*.pipeline.json` capturing parameters, file outputs, and metrics.

#### 5.2.1 Pipeline execution steps

For a run with prefix `P`, the runner performs the following steps:

1. **Reference**
   - Executes: `python -m dataset.longread.reference -o P --length ... --seed ... --preset ... --dup-segments ...`
   - Produces: `P.ref.fasta`, `P.ref.meta.json`

2. **Truth + haplotypes**
   - Executes: `python -m dataset.longread.truth --ref P.ref.fasta -o P --seed ... --num-snps ... --het-rate ... [--num-indels ...]`
   - Produces: `P.truth.vcf`, `P.hap1.fasta`, `P.hap2.fasta`, `P.truth.meta.json`
   - Compress/index truth VCF using `bcftools view/index` → `P.truth.vcf.gz`

3. **Oracle VCF**
   - Generates an oracle VCF from truth and compress/indexes it:
     - `P.oracle.vcf.gz` (+ index)
   - Oracle keeps the truth variant set and is used to isolate phasing behavior.

4. **Read simulation**
   - Executes: `python -m dataset.longread.readsim --hap1 P.hap1.fasta --hap2 P.hap2.fasta -o P ...`
   - Produces: `P.reads.fastq`, `P.reads.truth.tsv`, `P.reads.meta.json`

5. **Alignment**
   - Executes:
     - `minimap2 -a -x <map_preset> P.ref.fasta P.reads.fastq | samtools sort -o P.bam`
     - `samtools index P.bam`
   - Produces: `P.bam`, `P.bam.bai`

6. **Variant calling (called regime)**
   - Executes:
     - `bcftools mpileup -Ou -f P.ref.fasta -q <MAPQ> -Q <BASEQ> P.bam | bcftools call -mv -Oz -o P.called.vcf.gz`
     - `bcftools index -t P.called.vcf.gz`
   - Produces: `P.called.vcf.gz`, `P.called.vcf.gz.tbi`

7. **Phasing (oracle and/or called)**
   - Executes (for each chosen VCF input):
     - `python -m algorithms.cli.phase diploid-whats-bam --bam P.bam --vcf <VCF> --output-vcf <OUT> --output-prefix <PREFIX> ...`
   - Produces:
     - phased VCF: `*.phased.vcf`
     - WhatsHap run summary: `*.summary.json`

8. **Evaluation**
   - Executes:
     - `python -m benchmark.vcf_phase_eval --truth <truth> --pred <phased.vcf> --out <eval.json>`
   - Produces: `*.eval.json`

9. **Run report**
   - Writes: `P.pipeline.json` containing:
     - `params` (all effective parameters, including realism knobs and flags)
     - `steps` (file paths)
     - `callset` (calling precision/recall when called regime is enabled)
     - `phasing_runs` (oracle and/or called metrics + derived measures)

#### 5.2.2 Toolchain validation and robustness

The runner performs early validation:
- Ensures `minimap2`, `samtools`, and `bcftools` exist on `PATH` before executing.
- Uses piped-command helpers to propagate failures from upstream tools.
- Stores intermediate and final artifacts even when only partial regimes are enabled (oracle only, called only, or both).

---

### 5.3 WhatsHap integration: `diploid-whats-bam`

Phasing is implemented as a dedicated algorithm target:

- CLI: `python -m algorithms.cli.phase diploid-whats-bam`
- Driver: `algorithms/diploid/whatshap_bam_driver.py`

This driver is designed to benchmark WhatsHap while maintaining full control of:
- input selection (which variants are phased),
- read filtering (mapQ/baseQ thresholds),
- read selection (max coverage),
- output formatting (phased VCF + summary JSON).

#### 5.3.1 Vendored WhatsHap core

To avoid version drift and ensure controlled benchmarking, the project uses a vendored WhatsHap core under `vendor/whatshap_core`. The code ensures that `import whatshap` resolves to the vendored module and records provenance in `*.summary.json`:

- `whatshap_module`
- `whatshap_core_module`
- `whatshap_readselect_module`

This provides traceability between reported results and the exact phasing implementation used.

#### 5.3.2 Read filtering and selection

The driver applies filters to approximate practical usage:
- discard reads below `--min-mapq`
- discard allele observations below `--min-baseq`
- discard reads that cover fewer than two variants (WhatsHap read selection expects ≥2 covered variants)
- apply read selection with `whatshap.readselect.readselection` at `--max-coverage`

#### 5.3.3 Phased VCF output

The driver writes a phased VCF by:
- streaming the input VCF records,
- replacing `GT` with phased genotypes for phased heterozygous sites, and
- writing `PS` phase set tags to label blocks.

Homozygous genotypes do not require phasing and are generally preserved as unphased (`0/0`, `1/1`) or left unchanged depending on configuration.

---

### 5.4 Evaluation implementation: `benchmark.vcf_phase_eval`

Evaluation compares predicted phasing to truth using biallelic SNP records and reports:

- shared record counts (total, SNP, heterozygous)
- number of phased heterozygous records
- number of phase sets (fragmentation)
- phase accuracy under best block flip
- switch error rate (switches / adjacency denominator)
- per-phase-set breakdown

---

### 5.5 Indel-aware policy: SNP-only phasing/evaluation flags

When indels are introduced, representation differences between truth, called, and phased VCFs can break “exact match” comparisons. To keep evaluation meaningful, the platform supports:

- `--phase-snps-only`: filter the phasing input VCF to biallelic SNPs
- `--eval-snps-only`: filter the truth VCF to biallelic SNPs for evaluation

These flags are recorded in `*.pipeline.json` and surfaced into `aggregate.csv` so downstream plots always reflect whether SNP-only policy was applied.

---

### 5.6 Experiment automation and aggregation

#### 5.6.1 Experiment driver

- Entry point: `python -m benchmark.experiment_driver`
- Runs multiple pipeline instances across seeds and parameter sets.
- Produces per-experiment:
  - run artifacts (`*.pipeline.json` per run),
  - `aggregate.csv`,
  - plots under `plots/`.

#### 5.6.2 Aggregation contract

`benchmark.aggregate_pipeline_reports` reads each `*.pipeline.json` and outputs a standardized `aggregate.csv` containing:
- configuration columns (including realism knobs and SNP-only flags)
- callset metrics (precision/recall where applicable)
- phasing metrics for oracle/called regimes
- derived metrics (effective phased recall, shared het recall, etc.)

This CSV is the primary input for plotting and for generating tables in the report.

---

### 5.7 Realism knobs: where they are implemented and how they are surfaced

The platform’s realism features are not “one monolithic simulator.” Instead, realism is decomposed into stage-appropriate knobs, implemented where the phenomenon naturally arises:

1. **Reference-level stressors (`dataset.longread.reference`)**
   - Presets: `--preset {plain,toy,realistic}`
   - Duplications: `--dup-segments`, `--dup-len`, `--dup-min-gap`
   - Export: `P.ref.meta.json` containing `regions` and `duplications`

2. **Truth-level stressors (`dataset.longread.truth`)**
   - Indels: `--num-indels`, `--indel-min-len`, `--indel-max-len`, `--indel-het-rate`
   - Placement control: `--avoid-regions` (requires `--ref-meta`)
   - Export: `P.truth.meta.json` and a phased `P.truth.vcf(.gz)` when `--phased-truth` is used

3. **Read-level stressors (`dataset.longread.readsim`)**
   - Length distribution: `--len-model {uniform,lognormal}`, `--ln-mean`, `--ln-sigma`
   - Coverage dropout: `--start-model dropout`, `--dropout-fraction`, `--dropout-block-len`
   - Correlated error bursts: `--burst-prob`, `--burst-count`, `--burst-len`, `--burst-mult`
   - Export: `P.reads.meta.json` capturing the effective parameters used

4. **Runner-level policy stressors (`benchmark.longread_pipeline_runner`)**
   - VCF selection: `--vcf-source {oracle,called,both}`
   - Indel-safe evaluation: `--phase-snps-only`, `--eval-snps-only`
   - Recorded in: `P.pipeline.json.params` and aggregated into `aggregate.csv`

This structure makes it possible to validate each realism knob independently (unit validation), and later construct controlled sweeps where only one stressor changes at a time.