## 📂 `algorithms/cli` — Command-line interfaces

This repo exposes a small set of CLIs to make simulation, conversion, phasing, and benchmarking reproducible.

All commands can be run as modules (recommended):

- `python -m dataset.simulate`
- `python -m algorithms.cli.convert`
- `python -m algorithms.cli.phase`
- `python -m benchmark.benchmark_runner`
- `python -m benchmark.benchmark_accuracy`

---

## `python -m algorithms.cli.phase`

Run one phasing algorithm on a `.reads.npz` file.

### Subcommands

- `diploid-em`
- `diploid-mst`
- `diploid-whats`
- `polyploid-em`
- `polyploid-spectral`

### Common arguments

- `-i, --input` (required): input `.reads.npz`
- `--output-prefix` (required): prefix for outputs
- `--max-coverage`: (supported by some algorithms; relevant for WhatsHap read selection)

### Outputs

Each run typically writes:

- `{output-prefix}.haplotypes.tsv`
- `{output-prefix}.assignments.tsv`
- `{output-prefix}.summary.json` (when supported)

---

## `diploid-whats` (vendored WhatsHap)

### Matrix-mode (no VCF)

```bash
python -m algorithms.cli.phase diploid-whats       -i output/demo.reads.npz       --output-prefix output/demo_whats
```

Notes:
- This mode does not write a VCF.
- Internally it still performs read selection, then phases using `HapChatCore`.

### VCF-mode

```bash
python -m algorithms.cli.phase diploid-whats       -i output/demo.reads.npz       --vcf output/demo.vcf       --output-prefix output/demo_phased       --solver whatshap
```

VCF-mode options:

- `--vcf PATH`: enable VCF-mode (phase only heterozygous variants from GT)
- `--sample NAME`: which sample column in the VCF to use (default: first sample)
- `--output-vcf PATH`: output VCF path (default: `<output-prefix>.phased.vcf`)
- `--solver {whatshap,hapchat}`:
  - `whatshap` → PedigreeDPTable (default WhatsHap)
  - `hapchat` → HapChatCore
- `--recomb-rate FLOAT`: recombination rate (cM/Mb) used to compute recombination costs (PedigreeDPTable path)

---

## `python -m algorithms.cli.convert`

Convert between sparse TSV and dense NPZ read formats.

### Usage

```bash
python -m algorithms.cli.convert <input.{tsv|npz}> --output <output.{npz|tsv}>
```

Examples:

```bash
# sparse TSV -> NPZ
python -m algorithms.cli.convert output/demo.reads.sparse.tsv --output output/demo.reads.npz

# NPZ -> sparse TSV
python -m algorithms.cli.convert output/demo.reads.npz --output output/demo.reads.sparse.tsv
```
