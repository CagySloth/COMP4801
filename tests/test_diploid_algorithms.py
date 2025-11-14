# tests/test_diploid_algorithms.py

import numpy as np
import os, json
import tempfile
from algorithms.diploid import em_cluster, mst_phase

def _write_npz(tempdir, reads_matrix, filename="input.npz"):
    """Helper to write a reads matrix to NPZ in a temp directory."""
    path = os.path.join(tempdir, filename)
    np.savez_compressed(path, reads=reads_matrix)
    return path

def _load_hap_file(path):
    """Helper to load haplotypes.tsv into list of strings."""
    with open(path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    return lines

def test_diploid_em_perfect_phase(tmp_path):
    """Diploid EM: Test on a trivial case where two reads are exact haplotypes."""
    # Construct a simple scenario with 2 variants:
    # True haplotypes: [0,1] and [1,0]
    # Reads: one corresponds exactly to hap0, the other to hap1.
    reads_matrix = np.array([[0, 1],
                              [1, 0]], dtype=int)
    # Save NPZ
    npz_file = _write_npz(tmp_path, reads_matrix, "em_input.npz")
    out_prefix = os.path.join(tmp_path, "em_output")
    # Run EM clustering algorithm
    # Simulate command-line args
    args = type("Args", (), {})()  # empty object to set attributes
    args.input = npz_file
    args.output_prefix = out_prefix
    args.max_iters = 10
    args.tol_iters = 2
    args.seed = 42
    em_cluster.main(args)
    # Check output files
    hap_file = out_prefix + ".haplotypes.tsv"
    assign_file = out_prefix + ".assignments.tsv"
    summary_file = out_prefix + ".summary.json"
    assert os.path.exists(hap_file) and os.path.exists(assign_file), "EM outputs not found."
    # Load haplotypes
    haps = _load_hap_file(hap_file)
    # We expect 2 haplotype lines of length 2
    assert len(haps) == 2 and all(len(line) == 2 for line in haps), "Haplotype output format incorrect."
    # The two haplotypes should be complementary (e.g., "01" and "10" or vice versa)
    assert haps[0] != haps[1], "Diploid haplotypes should differ."
    combined = {haps[0], haps[1]}
    assert combined == {"01", "10"} or combined == {"10", "01"}, "Unexpected haplotype sequences for trivial case."
    # Check assignments
    with open(assign_file) as f:
        assigns = [line.strip() for line in f if line.strip()]
    # Should have 2 assignment lines (for 2 reads), and they should be different (one read to hap0, one to hap1)
    assert len(assigns) == 2
    assert assigns[0] != assigns[1], "In trivial case, reads should split between the two haplotypes."
    # MEC should be zero in summary (no errors, as each read exactly matches one haplotype)
    with open(summary_file) as f:
        summary = json.load(f)
    assert summary.get("MEC_total", None) == 0, "MEC should be 0 for perfect phasing."
    assert summary.get("assignment_accuracy", None) in (1.0, None), "Accuracy should be perfect (1.0) if truth known or not present."

def test_diploid_mst_basic(tmp_path):
    """Diploid MST: Test phasing on a simple 3-variant case."""
    # Design haplotypes: [0,0,1], [1,1,0] which differ at all three sites.
    # Construct reads that each cover all 3 variants (so overlap = 3):
    # Use 4 reads: 2 matching hap0 exactly, 2 matching hap1 exactly.
    hap0 = np.array([0,0,1], dtype=int)
    hap1 = np.array([1,1,0], dtype=int)
    reads_matrix = np.array([hap0, hap0, hap1, hap1], dtype=int)
    npz_file = _write_npz(tmp_path, reads_matrix, "mst_input.npz")
    out_prefix = os.path.join(tmp_path, "mst_output")
    # Prepare args for MST main
    args = type("Args", (), {})()
    args.input = npz_file
    args.output_prefix = out_prefix
    args.min_overlap = 1  # since reads span all, any overlap >0 is fine
    args.min_het_minor = 1
    mst_phase.main(args)
    # Load results
    hap_file = out_prefix + ".haplotypes.tsv"
    assign_file = out_prefix + ".assignments.tsv"
    summary_file = out_prefix + ".summary.json"
    assert os.path.exists(hap_file), "MST haplotypes output missing."
    haps = _load_hap_file(hap_file)
    # Expect 2 haplotypes of length 3
    assert len(haps) == 2 and all(len(h) == 3 for h in haps), "MST output haplotypes format incorrect."
    # The output haplotypes (as strings) should be some permutation of the true haplotypes "001" and "110"
    out_set = {haps[0], haps[1]}
    true_set = {"001", "110"}
    assert out_set == true_set, f"MST failed to recover the true haplotypes. Got {out_set}"
    # Assignments: expect first two reads assigned to one haplotype, next two to the other
    with open(assign_file) as f:
        assigns = [int(x.strip()) for x in f if x.strip()]
    assert assigns.count(0) == 2 and assigns.count(1) == 2, "MST should assign two reads to each haplotype cluster."
    # MEC check: since all reads exactly match a haplotype, MEC should be 0
    with open(summary_file) as f:
        summary = json.load(f)
    assert summary.get("MEC_total", None) == 0, "MST MEC should be 0 for perfectly clustered reads."
