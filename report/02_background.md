## 2. Background

### 2.1 Long-read sequencing and ONT-like characteristics

Long-read sequencing reads range from kilobases to megabases, by which is able to cover a large quantity of heterozygous variants within a single read. This could be more effective than short reads by enabling long-range phasing and reducing ambiguity. ONT data, however, exhibits properties including higher raw error rates than short reads, a mixture of substitution and indel errors, potential error clustering and burst, as well as non-uniform coverage due to sequencing and mapping biases, that could complicate the downstream analysis. Repetitive and duplicated regions further induce ambiguity in read mapping, which could degrade both variant calling and phasing.

### 2.2 WhatsHap in practical pipelines (VCF + BAM → phased VCF)

In a typical long-read workflow, reads are aligned to a reference to produce BAM/CRAM, variants are called to produce a VCF containing genotypes, and a phasing tool is utilised to assign heterozygous variants to haplotypes. WhatsHap extract aligned reads (BAM) and variants (VCF) as well as read-backed allele observations at variable positions. Its next procedure is to partition heterozygous variants into phase blocks and compute a phasing solution within each block under an error-correction objective. As a result, a phased VCF where heterozygous genotypes are phased (`0|1` / `1|0`) and annotated with phase set identifiers (`PS`) indicating blocks is produced. 

### 2.3 Necessity of end-to-end benchmarking

Phasing performance is undermined by unfavourable properties in multiple stages: 

- **Alignment:** Adverse mapping quality, ambiguous mappings in repeats and duplications, and alignment artifacts around indels. 
- **Variant calling:** Missing true variants (recall loss), spurious variants (precision loss), and genotype errors.
- **Phasing:** Block fragmentation and phase inconsistencies (switch errors) among phased variants.

To local the source of failures, this project evaluates WhatsHap under two framework: 

- **Oracle VCF:** Ground-truth variants used for phasing to isolate phasing behavior given a correct variant set.
- **Called VCF:** Variants produced by a variant caller used for phasing to measure realistic end-to-end performance, including caller limitations.

Attribution of inefficiency from phaser and caller is possible with this separation.

### 2.4 Metrics

This project reports metrics designed for long-read phasing:

- **Number of phase sets:** A measure of block fragmentation, fewer or larger blocks generally indicate better long-range phasing.
- **Switch error rate:** The fraction of adjacent phased heterozygous variants where predicted phase flips relative to truth.
- **Phasing rate on shared heterozygous sites:** The fraction of shared heterozygous variants that are phased.
- **Effective phased recall:** An end-to-end metric combining variant set overlap and correctness of phasing on truth heterozygous variants. Used to compare oracle vs called regimes consistently.