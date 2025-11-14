## 📂 `algorithms/cli` – Command-Line Interfaces

This module provides scriptable CLI entry points for phasing and data format conversion, enabling reproducible pipelines and benchmark integration.

---

### 📄 `phase.py`

Runs a specified phasing algorithm on an input `.npz` read file and outputs predicted haplotypes and assignments.

#### Usage

```bash
python -m algorithms.cli.phase [ALGORITHM] -i reads.npz -o output_prefix [options]
```

#### Positional Arguments

| Argument    | Description                                                                                       |
| ----------- | ------------------------------------------------------------------------------------------------- |
| `ALGORITHM` | Algorithm name. Must be one of: `diploid_em`, `diploid_mst`, `polyploid_em`, `polyploid_spectral` |

#### Options

| Option           | Type | Default              | Description                           |
| ---------------- | ---- | -------------------- | ------------------------------------- |
| `-i`, `--input`  | str  | —                    | Path to input `.npz` read file        |
| `-o`, `--output` | str  | —                    | Prefix for output files               |
| `-k`, `--ploidy` | int  | *Only for polyploid* | Number of haplotypes (polyploid only) |

#### Outputs

* `output_prefix.haplotypes.tsv` — predicted haplotypes
* `output_prefix.assignments.tsv` — read-to-haplotype assignments

#### Example

```bash
python -m algorithms.cli.phase diploid_em -i data/example.reads.npz -o results/em_run
```

---

### 📄 `convert.py`

Convert between `.tsv` and `.npz` read data formats.

#### Usage

```bash
python -m algorithms.cli.convert -i input_file [--to-tsv output.tsv | --to-npz output.npz]
```

#### Options

| Option          | Type | Description                                     |
| --------------- | ---- | ----------------------------------------------- |
| `-i`, `--input` | str  | Path to the input file (`.tsv` or `.npz`)       |
| `--to-tsv`      | str  | Output path for `.tsv` format (dense or sparse) |
| `--to-npz`      | str  | Output path for `.npz` format                   |

#### Notes

* If input is `.tsv`, this will save as `.npz`
* If input is `.npz`, it will save as `.reads.sparse.tsv`

#### Example

```bash
python -m algorithms.cli.convert -i data/example.reads.sparse.tsv --to-npz data/example.reads.npz
```