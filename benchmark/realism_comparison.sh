# 100 reads
for s in 0 1 2; do
  python -m benchmark.longread_pipeline_runner \
    --prefix output/exp_refpreset_compare/realistic_r100_s${s} \
    --seed ${s} \
    --ref-preset realistic \
    --vcf-source both \
    --platform ont \
    --ont-profile q20 \
    --num-reads 100 \
    --max-coverage 15
done

# 200 reads
for s in 0 1 2; do
  python -m benchmark.longread_pipeline_runner \
    --prefix output/exp_refpreset_compare/realistic_r200_s${s} \
    --seed ${s} \
    --ref-preset realistic \
    --vcf-source both \
    --platform ont \
    --ont-profile q20 \
    --num-reads 200 \
    --max-coverage 15
done