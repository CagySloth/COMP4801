# COMP4801: WhatsHap Benchmarking Suite

A final-year project codebase for **controlled benchmarking of WhatsHap-style phasing** under both idealized and practical long-read conditions.

The repository has two tracks:

1. **Legacy matrix track**
   - Simulate read/variant matrices (`.reads.npz`)
   - Run diploid/polyploid phasing algorithms
   - Score against truth haplotypes in TSV form

2. **Primary long-read end-to-end track**
   - Generate a synthetic reference and diploid truth
   - Simulate ONT-like reads (FASTQ)
   - Align reads (BAM)
   - Call variants (VCF)
   - Phase with the vendored WhatsHap core
   - Evaluate and aggregate results

The long-read track is the main research contribution of the project. It is designed to support:

- **oracle vs called attribution**
- **realism stressors** such as duplication, coverage dropout, burst errors, and indels
- **seed-controlled reproducibility**
- **traceable experiment outputs** through prefix-based artifact naming and machine-readable reports

---

## What is in this repository?

### Data generation
- `dataset.simulate` — legacy matrix-track simulator
- `dataset.longread.reference` — synthetic reference FASTA + metadata
- `dataset.longread.truth` — truth VCF + haplotype FASTAs + oracle VCF
- `dataset.longread.readsim` — ONT-like FASTQ simulation + read metadata

### Phasing
- `algorithms.cli.phase` — unified CLI for all phasing backends
- `algorithms.diploid.whatshap_driver` — matrix / VCF-mode vendored WhatsHap adapter
- `algorithms.diploid.whatshap_bam_driver` — BAM + VCF long-read phasing path
- `vendor/whatshap_core` — vendored WhatsHap core used by the project

### Benchmarking and experiments
- `benchmark.longread_pipeline_runner` — canonical single-run long-read pipeline
- `benchmark.vcf_phase_eval` — phased-VCF evaluation against truth
- `benchmark.aggregate_pipeline_reports` — aggregate many `*.pipeline.json` files into CSV
- `benchmark.experiment_driver` — run full experiment suites and generate plots

### Documentation and validation
- `docs/` — repository documentation
- `tests/` — unit and integration tests for the core project logic
- `report/` — final project report source and generated PDF

---

## Installation

### Prerequisites

- Python **3.11 or 3.12**
- `minimap2`
- `samtools`
- `bcftools`

Example installation on Ubuntu / Debian:

```bash
sudo apt-get update
sudo apt-get install -y minimap2 samtools bcftools
```

Check that the external tools are available:

```bash
minimap2 --version
samtools --version
bcftools --version
```

### Python environment

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -e .
pip install -e vendor/whatshap_core
```

Verify that the vendored WhatsHap package is the one being imported:

```bash
python - <<'PY'
import whatshap as wh
from whatshap import core, readselect
print('whatshap:', wh.__file__)
print('core:', core.__file__)
print('readselect:', readselect.__file__)
PY
```

---

## Quick start

### 1. Legacy matrix track

Generate a small matrix-track dataset:

```bash
python -m dataset.simulate \
  -p 2 -n 200 -r 200 -l 60 \
  -e 0.01 -m 0.0 \
  --seed 0 \
  -o output/demo
```

Phase with the vendored WhatsHap adapter:

```bash
python -m algorithms.cli.phase diploid-whats \
  -i output/demo.reads.npz \
  --vcf output/demo.vcf \
  --output-prefix output/demo_phased \
  --solver whatshap
```

Score TSV haplotypes:

```bash
python -m benchmark.benchmark_accuracy \
  --truth output/demo.haplotypes.tsv \
  --pred  output/demo_phased.haplotypes.tsv \
  --output output/demo_phased.accuracy.json
```

### 2. Long-read end-to-end track

Run the full long-read pipeline:

```bash
python -m benchmark.longread_pipeline_runner \
  --prefix output/lr_demo \
  --seed 0 \
  --ref-length 20000 \
  --num-snps 200 \
  --het-rate 0.8 \
  --num-reads 200 \
  --min-len 2000 \
  --max-len 6000 \
  --platform ont \
  --ont-profile q20 \
  --vcf-source both
```

This produces a family of prefix-based artifacts such as:

- `output/lr_demo.ref.fasta`
- `output/lr_demo.truth.vcf.gz`
- `output/lr_demo.oracle.vcf.gz`
- `output/lr_demo.reads.fastq`
- `output/lr_demo.bam`
- `output/lr_demo.called.vcf.gz`
- `output/lr_demo.ws*.phased.vcf`
- `output/lr_demo.ws*.eval.json`
- `output/lr_demo.pipeline.json`

### 3. Indel-enabled mode

If truth indels are enabled, restrict phasing and evaluation to SNPs:

```bash
python -m benchmark.longread_pipeline_runner \
  --prefix output/lr_indels \
  --seed 0 \
  --ref-length 80000 \
  --num-snps 800 \
  --num-indels 80 \
  --indel-min-len 1 \
  --indel-max-len 5 \
  --indel-het-rate 0.5 \
  --num-reads 200 \
  --min-len 2000 \
  --max-len 6000 \
  --platform ont \
  --ont-profile q20 \
  --vcf-source both \
  --phase-snps-only \
  --eval-snps-only
```

---

## Why the pipeline is structured this way

The long-read pipeline is intentionally:

- **modular** — each stage has a single responsibility
- **reproducible** — every run is seed-controlled
- **traceable** — all artifacts share a prefix and are summarized in `*.pipeline.json`
- **experiment-ready** — the experiment driver repeatedly calls the canonical runner and aggregates outputs into `aggregate.csv`

This design supports both software-engineering goals (maintainability, testing, debugging) and research goals (controlled attribution, realism sweeps, parameter tuning).

---

## Documentation map

Start with `docs/index.md`.

Recommended reading order:

1. `docs/longread_pipeline.md`
2. `docs/realism_knobs.md`
3. `docs/benchmark.md`
4. `docs/experiments.md`
5. `docs/validation.md`
6. `docs/io.md`

The remaining docs cover the dataset generators, algorithms, CLIs, and example commands.

---

## Repository layout

```text
COMP4801/
├── algorithms/          # Diploid & polyploid phasing algorithms + CLI
├── benchmark/           # Runners, evaluators, aggregation, experiments
├── dataset/             # Matrix and long-read data generation
├── docs/                # Repository documentation
├── examples/            # Helper scripts used during report / plot generation
├── report/              # Final report source and generated PDF
├── tests/               # Tests for core logic and adapters
└── vendor/              # Vendored WhatsHap core
```

---

## Notes

- The long-read research track is the current primary focus of the repository.
- The matrix track is still useful for fast debugging and algorithm checks.
- Structural variants are not currently a full first-class phasing target in the evaluation workflow; indels mainly act as realism stressors and are handled through SNP-only phasing/evaluation when needed.
