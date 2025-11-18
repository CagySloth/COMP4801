# tests/test_data_simulation_io.py

import os
import numpy as np
import io
import json
import tempfile

from dataset import simulate   # assuming simulate.py is in package path
from algorithms.io import parser, writer, ReadsData

def test_generate_diploid_data_basic():
    """Generate a small diploid dataset and verify haplotypes differ at least one position and reads align."""
    num_variants = 10
    num_reads = 5
    read_length = 5
    # Generate data
    haplotypes, reads = simulate.generate_diploid_data(
        num_variants=num_variants,
        num_reads=num_reads,
        read_length=read_length,
        error_rate=0.0,
        missing_rate=0.0,
        allow_monomorphic=False
    )
    # Check output shapes
    assert haplotypes.shape == (2, num_variants), "Diploid haplotypes shape should be (2, N)."
    # Ensure haplotypes are not identical (monomorphic), since allow_monomorphic=False
    assert not np.array_equal(haplotypes[0], haplotypes[1]), "Haplotypes should differ at >=1 variant when not allowing monomorphic."
    # Each read is a dict with 'indices' and 'values'
    assert isinstance(reads, list) and reads, "Reads should be a non-empty list of fragments."
    for read in reads:
        idx, vals = read["indices"], read["values"]
        # Read length consistency
        assert len(idx) == len(vals) == read_length, "Each read fragment should have length equal to read_length."
        # Indices should be sorted and within range
        assert idx == sorted(idx), "Read indices should be sorted."
        assert max(idx) < num_variants, "Read indices must be < num_variants."

def test_generate_polyploid_data_missing():
    """Generate a small polyploid dataset and ensure missing data is present and handled."""
    num_variants = 20
    ploidy = 3
    num_reads = 10
    read_length = 10
    # Use a higher missing rate
    haps, reads = simulate.generate_polyploid_data(
        num_variants=num_variants,
        ploidy=ploidy,
        num_reads=num_reads,
        read_length=read_length,
        error_rate=0.0,
        missing_rate=0.5,
        alpha=0.4,
        beta=0.4,
        allow_monomorphic=False
    )
    # Check shapes
    assert haps.shape == (ploidy, num_variants), "Polyploid haplotypes shape should be (P, N)."
    # Ensure at least one variant differs among haplotypes
    all_same = np.all(haps == haps[0], axis=0)
    assert not np.all(all_same), "Not all variants should be identical across haplotypes (monomorphic should be prevented)."
    # Check that some missing values (-1) are present in reads
    any_missing = any(-1 in frag["values"] for frag in reads)
    assert any_missing, "With high missing_rate, some reads should contain -1 values."

def test_parse_and_convert_roundtrip(tmp_path):
    """Test that sparse TSV to NPZ and back to TSV retains data."""
    # Create a dummy sparse TSV content
    tsv_content = "0\t0:1\t5:0\t7:1\n1\t5:0\t7:1\n2\t7:1\n"
    tsv_file = tmp_path / "test.sparse.tsv"
    tsv_file.write_text(tsv_content)
    # Parse sparse TSV
    fragments = parser.parse_sparse_tsv(str(tsv_file))
    # Fragments should be a list of dicts
    assert isinstance(fragments, list) and len(fragments) == 3, "Should parse 3 reads from TSV."
    # Convert fragments to ReadsData and save to NPZ
    reads_data = ReadsData.from_fragments(fragments)
    npz_path = tmp_path / "test.npz"
    reads_data.to_npz(npz_path)  # assuming ReadsData has to_npz method
    assert npz_path.exists(), "NPZ file was not created."
    # Load back the NPZ and convert to sparse TSV
    loaded = ReadsData.from_npz(npz_path)
    out_tsv_path = tmp_path / "roundtrip.tsv"
    loaded.to_sparse_tsv(out_tsv_path)
    # Now compare the original TSV content with round-tripped content (ignoring line order)
    roundtrip_content = out_tsv_path.read_text().strip().splitlines()
    original_lines = tsv_content.strip().splitlines()
    assert set(roundtrip_content) == set(original_lines), "Round-trip conversion TSV->NPZ->TSV did not preserve content."

def test_reads_loading_modes(tmp_path):
    """Test that load_reads auto-detects sparse vs dense TSV correctly."""
    # Prepare a dense TSV (each line is a full binary string of equal length)
    dense_content = "01\n11\n00\n"
    dense_file = tmp_path / "reads.dense.tsv"
    dense_file.write_text(dense_content)
    data_dense = parser.load_reads(dense_file)
    # Should result in a ReadsData with a dense matrix of shape (3 reads, 2 variants)
    assert data_dense.reads.shape == (3, 2), "Dense TSV with 3 lines of 2 chars should yield a 3x2 matrix."
    # Prepare a sparse TSV (with colons)
    sparse_content = "0\t1:1\t3:0\n1\t0:1\n"
    sparse_file = tmp_path / "reads.sparse.tsv"
    sparse_file.write_text(sparse_content)
    data_sparse = parser.load_reads(sparse_file)
    # The sparse loader should produce a ReadsData as well
    assert data_sparse.reads.shape[0] == 2, "Sparse TSV had 2 reads."
    assert data_sparse.reads.shape[1] == 4, "Sparse TSV with highest index 3 should have 4 variant columns."
    # Check that the content matches expectations (fill missing with -1)
    # For read 0: index 1=1, index 3=0 => positions 0 and 2 are missing (-1)
    first_read = data_sparse.reads[0]
    assert first_read.tolist() == [-1, 1, -1, 0], "Sparse parse didn't produce expected read vector."
