## I/O formats (quick reference)

This repo uses a mix of text and binary formats.

---

# Matrix track

## Truth haplotypes TSV
`<prefix>.haplotypes.tsv`

- Each line: `hap_id<TAB>{0/1 string}`
- Used for TSV-based accuracy scoring

## Reads NPZ
`<prefix>.reads.npz`

Dense allele matrix representation used by algorithms:
- `alleles`: shape (R, N), values in {-1,0,1}
- `positions`: optional site indices

## Reads sparse TSV
`<prefix>.reads.sparse.tsv`
Sparse fragment list representation.

---

# Long-read track

## Reference FASTA
`<prefix>.ref.fasta`

Single contig.

## Reference meta JSON
`<prefix>.ref.meta.json`

Contains:
- `regions`: list of complex intervals (start0,end0,type)
- `duplications`: list of duplication events

## Truth VCF
`<prefix>.truth.vcf` (+ `.gz/.tbi`)

Phased GT recommended (`0|1`, `1|0`).

If indels are enabled, the truth VCF includes indels; use SNP-only filters for phasing/eval.

## Oracle VCF
`<prefix>.oracle.vcf.gz`

Derived from truth:
- same sites
- GT bars replaced with slashes (`|` → `/`)
- PS removed

Used to isolate phasing from calling.

## Reads FASTQ
`<prefix>.reads.fastq`

Simulated long reads.

## Reads truth TSV
`<prefix>.reads.truth.tsv`

Per-read labels:
- hap (1/2)
- start0, ref length, observed length
- edit counts (subs/ins/dels)
- contig

## Reads meta JSON
`<prefix>.reads.meta.json`

Stores simulation parameters and summary stats.

## BAM
`<prefix>.bam` (+ `.bai`)

Alignment from minimap2/samtools.

## Called VCF
`<prefix>.called.vcf.gz` (+ `.tbi`)

bcftools mpileup/call output.

## Phased VCF
`<prefix>.ws*.phased.vcf`

Output of `diploid-whats-bam`.  
Includes `GT` and `PS` for phased het sites.

## Pipeline report JSON
`<prefix>.pipeline.json`

Single file capturing:
- inputs/outputs
- callset metrics
- phasing metrics (called/oracle)
- derived metrics (effective phased recall, etc.)

---

# SNP-only filtering

When working with indels, filter to biallelic SNPs using:

```bash
bcftools view -v snps -m2 -M2 -Oz -o out.snps.vcf.gz in.vcf.gz
bcftools index -t out.snps.vcf.gz
```

The pipeline runner’s `--phase-snps-only` automates this.
