## 📂 `dataset` — Synthetic Data Generator

The `dataset` package generates synthetic haplotypes and sequencing reads for benchmarking and regression tests.

The main entry point is:

- `dataset/simulate.py` (module: `dataset.simulate`)

---

## `python -m dataset.simulate`

Generate a synthetic dataset (diploid or polyploid) and write it to disk.

### Parameters

- `-p, --ploidy` (int, required): number of haplotypes (2 = diploid)
- `-n, --num-variants` (int, required): number of variant sites
- `-r, --num-reads` (int, required): number of reads/fragments
- `-l, --read-length` (int, required): length of each read (in variant sites, not bp)
- `-e, --error-rate` (float): flip probability for an observed allele (default `0.01`)
- `-m, --missing-rate` (float): probability that an allele is missing (`-1`) (default `0.0`)
- `--maf-alpha`, `--maf-beta` (float): Beta distribution parameters for minor allele frequency (default `1.0, 1.0`)
- `--allow-monomorphic`: allow sites that are always 0 or always 1
- `--seed` (int): random seed for reproducibility
- `-o, --output-prefix` (str, required): output prefix (directory will be created if needed)

### Output files

For all ploidies:

- `{prefix}.haplotypes.tsv` — ground truth haplotypes
- `{prefix}.reads.sparse.tsv` — sparse read fragments (`read_id \t idx:allele \t idx:allele ...`)
- `{prefix}.reads.npz` — dense reads matrix for algorithm drivers

For **diploid only** (`--ploidy 2`), additionally:

- `{prefix}.vcf` — a minimal, single-sample VCF with **unphased GT** derived from truth haplotypes
  - Homozygous sites: `0/0` or `1/1`
  - Heterozygous sites: `0/1`

### Example

```bash
python -m dataset.simulate       -p 2 -n 100 -r 50 -l 30       -e 0.01 -m 0.0       --seed 0       -o output/demo
```

---

## Notes on formats

- The `.reads.npz` file stores:
  - `reads`: an `R x N` dense matrix with values `{0, 1, -1}` (`-1` = missing)
  - `positions`: an `R x N` matrix of variant positions
    - For simulated data, this is typically a tiled `0..N-1` index grid.
    - Some drivers (e.g., WhatsHap VCF-mode) may override this to preserve original VCF coordinates after subsetting.
