## 📂 `benchmark` — Benchmarking and evaluation

This folder contains:

1) **Matrix-track benchmarking**
   - `benchmark.benchmark_runner`: parameter sweeps (simulate → phase → score)
   - `benchmark.benchmark_accuracy`: accuracy from TSV haplotypes

2) **VCF-track evaluation**
   - `benchmark.vcf_phase_eval`: phased VCF vs phased truth VCF (block-flip aware)

3) **Long-read end-to-end automation**
   - `benchmark.longread_pipeline_runner`: runs Steps 1–7 and writes a pipeline report

---

# 1) `python -m benchmark.benchmark_runner`

Runs many simulations and algorithms, writing a `benchmark_summary.json` + `.csv`.

Example:

```bash
python -m benchmark.benchmark_runner \
  --algorithms diploid-em diploid-mst diploid-whats \
  --ploidy 2 \
  --num-variants 500 \
  --num-reads 500 \
  --read-length 60 \
  --error-rate 0.01 \
  --missing-rate 0.0 \
  --num-runs 2 \
  --outdir output/bench_smoke
```

Notes:
- For `diploid-whats`, the runner enables VCF-mode when a simulated `{prefix}.vcf` exists.
- This runner targets the **matrix track** (NPZ/TSV), not the BAM long-read track.

---

# 2) `python -m benchmark.benchmark_accuracy`

Scores predicted TSV haplotypes against truth TSV haplotypes.

```bash
python -m benchmark.benchmark_accuracy \
  --truth output/demo.haplotypes.tsv \
  --pred  output/demo_phased.haplotypes.tsv \
  --output output/demo_phased.accuracy.json
```

---

# 3) `python -m benchmark.vcf_phase_eval`

Evaluates predicted phased VCF vs phased truth VCF (block-flip aware).

Inputs:
- `--truth`: truth VCF/VCF.GZ with phased GT (e.g. `0|1`, `1|0`)
- `--pred`: predicted VCF with phased GT for het sites and PS tags

Example:

```bash
python -m benchmark.vcf_phase_eval \
  --truth output/lr_demo.truth.vcf.gz \
  --pred  output/lr_demo.ws.phased.vcf \
  --out   output/lr_demo.ws.eval.json
```

---

# 4) `python -m benchmark.longread_pipeline_runner`

Runs the long-read end-to-end pipeline (Steps 1–7):

1. reference FASTA
2. truth phased VCF + hap FASTAs
3. reads FASTQ
4. align to BAM (minimap2 + samtools)
5. call variants (bcftools mpileup + call)
6. phase called VCF using BAM (`diploid-whats-bam`)
7. evaluate phased VCF vs truth VCF (`vcf_phase_eval`)

Example:

```bash
python -m benchmark.longread_pipeline_runner \
  --prefix output/lr_demo \
  --seed 0 \
  --ref-length 20000 \
  --num-snps 200 \
  --het-rate 0.8 \
  --num-reads 200 \
  --min-len 2000 \
  --max-len 6000
```

Outputs:
- `{prefix}.pipeline.json` (the pipeline report)
- plus all intermediate files under the same prefix
