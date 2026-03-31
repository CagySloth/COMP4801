## 1. Introduction

Haplotype phasing assigns heterozygous variants to their chromosome of origin, enabling downstream analyses such as compound heterozygosity detection, allele-specific expression, and improved variant interpretation. Long-read sequencing provides reads that span multiple variants and can therefore improve both the contiguity and correctness of phased haplotypes compared to short-read data. In practice, however, phasing quality is not determined by the phasing algorithm alone: it is affected by read alignment quality, variant calling accuracy, and sequencing characteristics such as coverage bias, repetitive/duplicated regions, and structured error patterns.

WhatsHap is one of the most widely adopted phasing tools for long reads. While its algorithmic foundations are well studied, its behavior in practical end-to-end pipelines can vary significantly depending on upstream errors and data realism. A controlled benchmarking environment—where individual factors can be varied independently and runs can be reproduced reliably—is necessary to understand WhatsHap’s capabilities and limitations in realistic settings, and to identify opportunities for optimization.

### 1.1 Project goal

This project aims to:

1. Build a reproducible end-to-end benchmarking platform that mimics key characteristics of ONT-like long-read phasing workflows.
2. Evaluate WhatsHap’s phasing capability under increasingly realistic conditions and quantify how different factors degrade performance.
3. Identify optimization directions or configuration recommendations for practical WhatsHap usage.

### 1.2 Contributions

- An end-to-end long-read simulation and benchmarking pipeline producing reference, truth variants/haplotypes, reads, alignments, variants, phased outputs, and standardized JSON/CSV reports.
- A dual evaluation approach using **oracle** (ground-truth) and **called** variant sets to separate phasing limitations from variant-calling limitations.
- Implemented realism “knobs” (e.g., duplicated regions, dropout, correlated error bursts, indel-rich truth) to isolate impacts on alignment/calling/phasing.
- Automated experiment running, aggregation, and plotting to support systematic studies and report-ready results.

### 1.3 Report organization

Section 2 introduces background on long-read phasing, WhatsHap, and evaluation metrics. Sections 3–6 describe requirements, design, implementation, and validation. Section 10 presents experimental results and optimization investigations, followed by discussion and conclusions.