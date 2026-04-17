\clearpage
\pagenumbering{arabic}
\setcounter{page}{1}

## 1. Introduction

Haplotype phasing assigns heterozygous variants to their chromosome of origin, enabling downstream analyses such as compound heterozygosity detection, allele-specific expression, and improved variant interpretation. Long-read sequencing provides reads that span multiple variants and can therefore improve both the contiguity and correctness of phased haplotypes compared to short-read data. In practice, however, phasing quality is not determined by the phasing algorithm alone: it is affected by read alignment quality, variant calling accuracy, and sequencing characteristics such as coverage bias, repetitive/duplicated regions, and structured error patterns.

WhatsHap is one of the most widely adopted phasing tools for long reads. While its algorithmic foundations are well studied, its behavior in practical end-to-end pipelines can vary significantly depending on upstream errors and data realism. To understand WhatsHap's capabilities and limitations in realistic and practical situations, it is necessary to have a benchmarking platform, where individual factors and parameters can be controlled and varied independently, benchmarking runs are efficient, reliable and reproducible.

### 1.1 Project goal

This final year project aims to:

1. Build a data synthesizing pipeline that generate genome data that simulates ONT long-read characteristics and errors.
2. Build a reproducible end-to-end benchmarking pipeline that simulates ONT-like long-read phasing workflows.
3. Evaluate WhatsHap’s phasing capability under increasingly realistic and difficult conditions and quantify how different factors degrade performance.
4. Explore optimization directions or parameter configurations for practical WhatsHap usage.

### 1.2 Contributions

- An end-to-end long-read simulation and benchmarking pipeline capable of generating reference, truth variants/haplotypes, reads, alignments, variants, phased outputs, and standardized JSON/CSV reports.
- A dual evaluation approach using **oracle** (ground-truth) and **called** variant sets to separate phasing limitations from variant-calling limitations.
- Implemented realism “knobs” (e.g., duplicated regions, dropout, correlated error bursts, indel-rich truth) to isolate impacts on alignment/calling/phasing.
- Automated experiment running, aggregation, and plotting to support systematic studies and report-ready results.

### 1.3 Report organization

Section 2 covers background on long-read phasing, WhatsHap, and evaluation metrics. Sections 3–9 cover functional requirements, system design, implementation, validation, configuration management, quality assurance, and project management. Section 10 covers experimental results and optimization investigations, followed by discussion and conclusions in section 11 and 12.