# Validation and acceptance checks

Validation in this project focuses on two goals:

1. **software correctness**
2. **behavioral sanity of the benchmark**

The project does not attempt to prove full biological realism for every knob. Instead, it checks that:

- the pipeline runs successfully end-to-end
- outputs are correctly structured
- realism knobs affect metrics in reasonable ways
- oracle vs called behavior remains interpretable
- runs are reproducible under fixed seeds

---

## Test layers in the repository

### Unit / integration tests
Located in `tests/`.

Examples include:
- `test_whatshap_adapter.py`
- `test_vcf_phase_eval_unit.py`
- `test_aggregate_pipeline_reports_flags.py`
- `test_reference_duplications.py`

These cover critical logic such as:
- adapter behavior
- evaluator behavior
- aggregation behavior
- duplication metadata generation

### Smoke validation
The long-read pipeline is also validated through small end-to-end runs that exercise:
- reference generation
- truth/oracle creation
- read simulation
- alignment
- calling
- phasing
- evaluation
- report generation

---

## Recommended acceptance checks

### 1. Baseline smoke run

```bash
python -m benchmark.longread_pipeline_runner \
  --prefix output/acc_baseline \
  --seed 0 \
  --ref-length 80000 \
  --num-snps 800 \
  --num-reads 200 \
  --min-len 2000 \
  --max-len 6000 \
  --platform ont \
  --ont-profile q20 \
  --vcf-source both
```

Check that:
- all expected artifacts exist
- `*.pipeline.json` is produced
- `*.ws*.eval.json` contains expected keys
- aggregation can produce a non-empty `aggregate.csv`

### 2. One-knob realism checks
Run small sanity cases for:
- duplication
- dropout
- bursts
- indels

Confirm that:
- metadata files reflect the chosen settings
- the pipeline still completes
- metric changes are at least directionally plausible

### 3. Indel mode sanity
When `num_indels > 0`, verify that SNP-only phasing/evaluation is enabled and that filtered truth/prediction VCFs are created as expected.

---

## What to inspect in `*.pipeline.json`

### Callset sanity
- `call_precision`
- `call_recall`

### Phasing sanity
- `oracle_effective_phased_recall`
- `called_effective_phased_recall`
- `oracle_num_phase_sets`
- `called_num_phase_sets`
- `oracle_switch_error`
- `called_switch_error`

### Derived narrative metrics
- `shared_het_recall`
- `phasing_rate_on_shared_het`
- `phase_accuracy`

---

## Reproducibility check

Run the same command twice with the same `--seed` and compare key metrics.

For deterministic runs, the numerical metrics should match unless you intentionally include timestamps or machine-specific metadata in the output.

---

## Practical validation philosophy

The benchmark is accepted when:
- the code is correct
- the outputs are traceable
- the experiments are reproducible
- realism knobs behave as controlled stressors
- oracle vs called attribution remains meaningful
