## 📂 CLIs — Running the repo reproducibly

This repo supports running everything as Python modules (recommended):

Matrix track:
- `python -m dataset.simulate`
- `python -m algorithms.cli.convert`
- `python -m algorithms.cli.phase`
- `python -m benchmark.benchmark_runner`
- `python -m benchmark.benchmark_accuracy`

Long-read track:
- `python -m dataset.longread.reference`
- `python -m dataset.longread.truth`
- `python -m dataset.longread.readsim`
- `python -m benchmark.longread_pipeline_runner`
- `python -m benchmark.vcf_phase_eval`

If you installed with `pip install -e .`, you also get console scripts:

- `simulate`, `phase`, `convert`, `benchmark`, `accuracy`
- `longread-reference`, `longread-truth`, `longread-readsim`

---

## `python -m algorithms.cli.phase`

Run one phasing algorithm.

### Subcommands

- Diploid:
  - `diploid-em`
  - `diploid-mst`
  - `diploid-whats`
  - `diploid-whats-bam`
- Polyploid:
  - `polyploid-em`
  - `polyploid-spectral`

### Common outputs

Most algorithms write:

- `{output-prefix}.haplotypes.tsv`
- `{output-prefix}.assignments.tsv`
- `{output-prefix}.summary.json` (when supported)

---

## `diploid-whats` (vendored WhatsHap)

### Matrix-mode (no VCF)

```bash
python -m algorithms.cli.phase diploid-whats \
  -i output/demo.reads.npz \
  --output-prefix output/demo_whats
```

### VCF-mode

```bash
python -m algorithms.cli.phase diploid-whats \
  -i output/demo.reads.npz \
  --vcf output/demo.vcf \
  --output-prefix output/demo_phased \
  --solver whatshap
```

VCF-mode options:

- `--vcf PATH`: enable VCF-mode (phase only heterozygous variants from GT)
- `--sample NAME`: which sample column in the VCF to use (default: first sample)
- `--output-vcf PATH`: output VCF path (default: `<output-prefix>.phased.vcf`)
- `--solver {whatshap,hapchat}`:
  - `whatshap` → PedigreeDPTable (default WhatsHap DP)
  - `hapchat` → HapChatCore (MEC)
- `--recomb-rate FLOAT`: recombination rate (cM/Mb) used for recombination costs (PedigreeDPTable path)

---

## `diploid-whats-bam` (BAM + VCF)

```bash
python -m algorithms.cli.phase diploid-whats-bam \
  --bam output/lr_demo.bam \
  --vcf output/lr_demo.called.vcf.gz \
  --output-prefix output/lr_demo.ws \
  --output-vcf output/lr_demo.ws.phased.vcf
```

Options:

- `--max-coverage INT` (read selection)
- `--min-mapq INT` (filter low MAPQ reads)
- `--min-baseq INT` (filter low base quality observations)

---

## `python -m algorithms.cli.convert`

Convert between sparse TSV and dense NPZ read formats.

```bash
python -m algorithms.cli.convert <input.{tsv|npz}> --output <output.{npz|tsv}>
```

Examples:

```bash
python -m algorithms.cli.convert output/demo.reads.sparse.tsv --output output/demo.reads.npz
python -m algorithms.cli.convert output/demo.reads.npz --output output/demo.reads.sparse.tsv
```
