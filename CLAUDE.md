# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a haplotype phasing benchmark platform for COMP4801. It evaluates multiple phasing algorithms (diploid and polyploid) using synthetic read data. The project uses Python 3.11+ with NumPy-based dense matrix operations and includes WhatsHap integration via Cython bindings.

## Development Setup

```bash
# Create virtual environment
python3.11 -m venv .venv
source .venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Or install as editable package
pip install -e .

# Install dev dependencies
pip install -e ".[dev]"
```

### Building WhatsHap Core (Required for diploid-whats algorithm)

The WhatsHap integration requires building Cython extensions:

```bash
cd vendor/whatshap_core
pip install -e .
cd ../..
```

This compiles the C++ backend and makes `whatshap.core` importable. The build uses:
- Cython to wrap C++ code
- C++11 standard
- Multiple source files from `vendor/whatshap_core/src/`

## Testing

```bash
# Run all tests
pytest tests/

# Run specific test file
pytest tests/test_diploid_algorithms.py

# Run specific test
pytest tests/test_diploid_algorithms.py::test_diploid_em_basic

# Verbose output with print statements
pytest tests/ -v -s

# Run only WhatsHap-related tests
pytest tests/test_diploid_whatshap.py tests/test_whatshap_adapter.py tests/test_pipeline_whatshap_cli.py
```

Tests are located in `tests/` and cover:
- Data simulation and I/O (`test_data_simulation_io.py`)
- Algorithm correctness (`test_diploid_algorithms.py`, `test_polyploid_algorithms.py`)
- WhatsHap integration (`test_diploid_whatshap.py`, `test_whatshap_adapter.py`, `test_pipeline_whatshap_cli.py`)
- End-to-end pipelines (`test_pipeline_integration.py`)

## Code Architecture

### Data Flow Pipeline

The canonical workflow follows this sequence:

1. **Simulation** → generates ground truth haplotypes and synthetic reads
2. **Conversion** → TSV (sparse) ↔ NPZ (dense) format conversion
3. **Phasing** → algorithm runs on read data to infer haplotypes
4. **Evaluation** → compare predicted vs ground truth haplotypes

### Core Data Structure: ReadsData

All algorithms operate on `ReadsData` ([algorithms/io/reads_data.py](algorithms/io/reads_data.py)):

```python
@dataclass
class ReadsData:
    reads: np.ndarray       # shape (num_reads, read_length), values: {0, 1, -1}
    positions: np.ndarray   # shape (num_reads, read_length), variant indices
    num_variants: int       # total variant positions
```

**Key conventions:**
- `-1` represents missing data (gaps in read coverage)
- `0/1` are allele values
- `reads` is a dense matrix; sparse data is filled with `-1`

### Algorithm Modules

All phasing algorithms follow a standard pattern:

**Location:** `algorithms/diploid/` or `algorithms/polyploid/`

**Structure:**
1. Load data via `ReadsData.from_npz()`
2. Initialize haplotypes (k-means++, random, or graph-based)
3. Iterative refinement (EM, MST, spectral clustering)
4. Write outputs via `algorithms/io/writer.py`

**Expected outputs:**
- `<prefix>.haplotypes.tsv` - inferred haplotype matrix
- `<prefix>.assignments.tsv` - read-to-haplotype assignments
- `<prefix>.summary.json` - metadata (iterations, MEC score, runtime)

### CLI Entry Points

The unified CLI is in [algorithms/cli/phase.py](algorithms/cli/phase.py) with subcommands:

- `diploid-em` → [algorithms/diploid/em.py](algorithms/diploid/em.py)
- `diploid-mst` → [algorithms/diploid/mst.py](algorithms/diploid/mst.py)
- `diploid-whats` → [algorithms/diploid/whatshap_driver.py](algorithms/diploid/whatshap_driver.py)
- `polyploid-em` → [algorithms/polyploid/em.py](algorithms/polyploid/em.py)
- `polyploid-spectral` → [algorithms/polyploid/spectral.py](algorithms/polyploid/spectral.py)

Each algorithm module defines a `main(args)` function called by the dispatcher.

### WhatsHap Integration Architecture

**Purpose:** Leverage WhatsHap's production-grade DP solver for MEC optimization.

**Components:**
1. **Adapter** ([algorithms/diploid/whatshap_adapter.py](algorithms/diploid/whatshap_adapter.py))
   - Converts `ReadsData` → WhatsHap's `ReadSet` format
   - Handles dense-to-sparse conversion (skip `-1` values)
   - **Critical:** Must call `read.sort()` and `readset.sort()` before phasing

2. **Driver** ([algorithms/diploid/whatshap_driver.py](algorithms/diploid/whatshap_driver.py))
   - Orchestrates the phasing workflow
   - Calls WhatsHap's `HapChatCore` solver
   - Extracts superreads and converts back to standard format

3. **Vendor Code** ([vendor/whatshap_core/](vendor/whatshap_core/))
   - Full WhatsHap source with custom `setup.py`
   - Built via Cython + C++11 backend
   - Exposes `whatshap.core` module with `ReadSet`, `Read`, `HapChatCore`

**Debugging WhatsHap:**
- Ensure reads are sorted by position before phasing
- Check that quality scores are reasonable (typically 40)
- Verify `mapq` values are set (default 60)
- Use `readset.sort()` after building to avoid cryptic errors

### I/O Modules

**Parser** ([algorithms/io/parser.py](algorithms/io/parser.py)):
- `parse_sparse_tsv()` - reads `.reads.sparse.tsv` format
- `load_reads()` - dispatches to NPZ or TSV loader

**Writer** ([algorithms/io/writer.py](algorithms/io/writer.py)):
- `write_haplotypes_tsv()` - outputs phased haplotypes
- `write_assignments_tsv()` - outputs read assignments
- `write_summary_json()` - outputs metadata

**Formats:**
- **Sparse TSV:** `read_id \t pos:allele \t pos:allele ...`
- **NPZ:** Compressed NumPy archive with `reads` and `positions` arrays
- **Haplotypes TSV:** One haplotype per row, space-separated 0/1 values

### Evaluation Metrics

[algorithms/eval/metrics.py](algorithms/eval/metrics.py) implements:

- **MEC (Minimum Error Correction):** Count of read positions that disagree with assigned haplotype
- **Switch Error Rate:** Detects phase switches in predictions
- **Hamming Distance:** Position-wise differences between haplotypes

The benchmark runner ([benchmark/benchmark_accuracy.py](benchmark/benchmark_accuracy.py)) compares predictions against ground truth with automatic label permutation to handle phase ambiguity (e.g., swapping haplotype 0 and 1).

## Common Commands

### Running Algorithms

```bash
# Simulate data
python dataset/simulate.py -p 2 -n 1000 -r 5000 -o sim/test

# Convert to NPZ
python -m algorithms.cli.convert -i sim/test.reads.sparse.tsv --to-npz sim/test.npz

# Phase with diploid EM
python -m algorithms.cli.phase diploid-em -i sim/test.npz --output-prefix sim/test.out

# Phase with WhatsHap
python -m algorithms.cli.phase diploid-whats -i sim/test.npz --output-prefix sim/test.wh

# Evaluate accuracy
python benchmark/benchmark_accuracy.py \
  --truth sim/test.haplotypes.tsv \
  --pred sim/test.out.haplotypes.tsv \
  --output sim/test.accuracy.json
```

### Benchmarking

```bash
# Full parameter sweep
python benchmark/benchmark_runner.py \
  --algorithms diploid_em diploid_mst polyploid_em \
  --ploidy 4 \
  --num-variants 1000 \
  --num-reads 3000 \
  --vary error_rate \
  --vary-values 0.005 0.01 0.02 \
  --num-runs 3 \
  --outdir results/sweep

# Results saved to results/sweep/benchmark_summary.json
```

### Shell Wrappers

```bash
# Quick simulation
./scripts/simulate.sh

# Full pipeline (simulate + phase + evaluate)
./scripts/run_pipeline.sh

# Benchmark all algorithms
./scripts/benchmark_all.sh
```

## Package Configuration

The project is configured via `pyproject.toml` with CLI entry points:

- `simulate` → `dataset.simulate:main`
- `phase` → `algorithms.cli.phase:main`
- `convert` → `algorithms.cli.convert:main`
- `benchmark` → `benchmark.benchmark_runner:main`
- `accuracy` → `benchmark.benchmark_accuracy:main`

After `pip install .`, these become available as shell commands.

## Important Implementation Notes

### Algorithm Convergence

EM-based algorithms use two stopping criteria:
- `--max-iters`: Maximum iterations (default 20)
- `--tol-iters`: Stop if assignments unchanged for N iterations (default 5)

### Missing Data Handling

All algorithms must handle `-1` (missing) values:
- Skip when computing distances/likelihoods
- Only count overlapping positions between reads and haplotypes
- Use global majority vote to fill missing positions in final haplotypes

### Random Seeds

All stochastic algorithms accept `--seed` parameter for reproducibility. This is critical for benchmarking and testing.

### WhatsHap-Specific Notes

- WhatsHap expects **sorted reads** by first variant position
- Use `read.sort()` after adding variants to each read
- Use `readset.sort()` after adding all reads
- Quality scores (Q) are Phred-scaled: Q=40 means 0.01% error probability
- MAPQ represents mapping quality; use 60 for high-confidence alignments

## Git Workflow

Recent work focuses on WhatsHap integration. When modifying WhatsHap code:

1. Changes to adapter logic → [algorithms/diploid/whatshap_adapter.py](algorithms/diploid/whatshap_adapter.py)
2. Changes to phasing workflow → [algorithms/diploid/whatshap_driver.py](algorithms/diploid/whatshap_driver.py)
3. Vendor updates → rebuild with `cd vendor/whatshap_core && pip install -e .`
4. Always run `pytest tests/test_diploid_whatshap.py` after changes
