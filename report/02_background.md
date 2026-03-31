## 2. Background

### 2.1 Long-read sequencing and ONT-like characteristics

Long-read sequencing reads can span kilobases to megabases, often covering many heterozygous variants within a single read. This enables long-range phasing and reduces ambiguity compared to short reads. ONT data, however, exhibits properties that complicate downstream analysis, including higher raw error rates than short reads, a mixture of substitution and indel errors, potential error clustering (bursts), and non-uniform coverage due to sequencing and mapping biases. Repetitive and duplicated regions further introduce ambiguity in read mapping, which can degrade both variant calling and phasing.

### 2.2 WhatsHap in practical pipelines (VCF + BAM → phased VCF)

In a typical long-read workflow, reads are aligned to a reference to produce BAM/CRAM, variants are called to produce a VCF containing genotypes, and a phasing tool assigns heterozygous variants to haplotypes. WhatsHap takes aligned reads (BAM) and variants (VCF) and extracts read-backed allele observations at variant positions. It then partitions heterozygous variants into phase blocks and computes a phasing solution within each block under an error-correction objective. The output is a phased VCF where heterozygous genotypes are phased (`0|1` / `1|0`) and annotated with phase set identifiers (`PS`) indicating blocks.

### 2.3 Why end-to-end benchmarking is necessary

Phasing performance depends on multiple stages:

- **Alignment:** mapping quality, ambiguous mappings in repeats/duplications, and alignment artifacts around indels.
- **Variant calling:** missing true variants (recall loss), spurious variants (precision loss), and genotype errors.
- **Phasing:** block fragmentation and phase inconsistencies (switch errors) among phased variants.

To understand where failures originate, this project evaluates WhatsHap under two regimes:

- **Oracle VCF:** uses ground-truth variants to isolate phasing behavior given a correct variant set.
- **Called VCF:** uses variants produced by a caller to measure realistic end-to-end performance, including caller limitations.

This separation makes it possible to distinguish “phaser-limited” regimes (phasing errors dominate) from “caller-limited” regimes (missing/incorrect variants dominate).

### 2.4 Metrics

This project reports metrics designed for long-read phasing:

- **Number of phase sets:** a measure of block fragmentation (fewer/larger blocks generally indicate better long-range phasing).
- **Switch error rate:** the fraction of adjacent phased heterozygous variants where predicted phase flips relative to truth.
- **Phasing rate on shared heterozygous sites:** the fraction of shared heterozygous variants that are phased.
- **Effective phased recall:** an end-to-end metric combining variant set overlap and correctness of phasing on truth heterozygous variants (used to compare oracle vs called regimes consistently).