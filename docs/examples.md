# Example commands

This page collects short, copyable commands for common tasks.

---

## Matrix-track demo

```bash
PREFIX=output/demo

python -m dataset.simulate \
  -p 2 -n 200 -r 200 -l 60 -e 0.01 -m 0.0 \
  --seed 0 \
  -o "$PREFIX"

python -m algorithms.cli.phase diploid-whats \
  -i "$PREFIX.reads.npz" \
  --vcf "$PREFIX.vcf" \
  --output-prefix "${PREFIX}_phased" \
  --solver whatshap

python -m benchmark.benchmark_accuracy \
  --truth "$PREFIX.haplotypes.tsv" \
  --pred  "${PREFIX}_phased.haplotypes.tsv" \
  --output "${PREFIX}_phased.accuracy.json"
```

---

## Long-read baseline run

```bash
python -m benchmark.longread_pipeline_runner \
  --prefix output/lr_demo \
  --seed 0 \
  --ref-length 80000 \
  --num-snps 800 \
  --het-rate 0.8 \
  --num-reads 200 \
  --min-len 2000 \
  --max-len 6000 \
  --platform ont \
  --ont-profile q20 \
  --vcf-source both
```

---

## Long-read hard scenario starter

```bash
python -m benchmark.longread_pipeline_runner \
  --prefix output/acc_hard \
  --seed 0 \
  --ref-length 120000 \
  --num-snps 1200 \
  --ref-preset realistic \
  --dup-segments 5 \
  --dup-len 3000 \
  --dup-min-gap 500 \
  --start-model dropout \
  --dropout-fraction 0.1 \
  --dropout-block-len 1000 \
  --burst-prob 0.6 \
  --burst-count 2 \
  --burst-len 300 \
  --burst-mult 8 \
  --num-indels 120 \
  --indel-min-len 1 \
  --indel-max-len 5 \
  --indel-het-rate 0.5 \
  --num-reads 300 \
  --min-len 2000 \
  --max-len 6000 \
  --platform ont \
  --ont-profile q20 \
  --vcf-source both \
  --phase-snps-only \
  --eval-snps-only
```

---

## Aggregate many run reports

```bash
python -m benchmark.aggregate_pipeline_reports \
  --root output/experiments/01_depth \
  --out  output/experiments/01_depth/aggregate.csv
```

---

## Run selected experiment families

```bash
python -m benchmark.experiment_driver \
  --outdir output/experiments_final \
  --seeds 0 1 2 \
  --only depth,dropout,interaction,optimize
```

---

## Reality check section

```bash
python -m benchmark.experiment_driver \
  --outdir output/experiments_final \
  --seeds 0 1 2 \
  --only reality \
  --real-fastq /path/to/real.fastq.gz \
  --real-bam /path/to/real.bam
```
