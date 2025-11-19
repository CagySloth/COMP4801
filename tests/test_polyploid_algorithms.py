# tests/test_polyploid_algorithms.py

import numpy as np
import os, json
from algorithms.polyploid import em as poly_em, spectral as poly_spec

def test_polyploid_em_trivial(tmp_path):
    """Polyploid EM: Test on a trivial 3-ploidy case where each read is exactly one haplotype."""
    # True haplotypes for 3-ploidy, 2 variants each (distinct patterns):
    true_haps = np.array([[0,0], [0,1], [1,1]], dtype=int)
    # Create one read per haplotype covering all variants (so each read equals a true haplotype)
    reads_matrix = true_haps.copy()
    npz_file = os.path.join(tmp_path, "poly3_input.npz")
    np.savez_compressed(npz_file, reads=reads_matrix)
    out_prefix = os.path.join(tmp_path, "polyem_output")
    # Run polyploid EM
    args = type("Args", (), {})()
    args.input = npz_file
    args.ploidy = 3
    args.output_prefix = out_prefix
    args.max_iters = 20
    args.tol_iters = 5
    args.seed = 1
    poly_em.main(args)
    # Check outputs
    hap_file = out_prefix + ".haplotypes.tsv"
    assign_file = out_prefix + ".assignments.tsv"
    summary_file = out_prefix + ".summary.json"
    assert os.path.exists(hap_file), "Polyploid EM haplotypes output missing."
    haps_out = []
    with open(hap_file) as f:
        for line in open(hap_file):
            if line.strip():
                seq = line.strip().split("\t")[1]
                haps_out.append(seq)
    # Expect 3 haplotypes of length 2
    assert len(haps_out) == 3 and all(len(h) == 2 for h in haps_out), "Polyploid EM output format incorrect."
    # They should match the true haplotypes (order can differ)
    out_set = set(haps_out)
    true_set = {"00", "01", "11"}
    assert out_set == true_set, f"Polyploid EM did not find correct haplotypes. Got {out_set}"
    # Assignments: with one read per haplotype, each should be in its own cluster
    with open(assign_file) as f:
        assigns = [line.strip() for line in f if line.strip()]
    # We have 3 reads, expecting assignments like 0,1,2 in some order
    assert set(assigns) == {"0","1","2"}, "Each read should be assigned to a unique haplotype cluster."
    # MEC should be 0 (each read perfectly matches a haplotype)
    with open(summary_file) as f:
        summary = json.load(f)
    assert summary.get("MEC_total", None) == 0, "Polyploid EM MEC should be 0 for perfectly partitioned reads."

def test_polyploid_spectral_basic(tmp_path):
    """Polyploid Spectral: Test clustering on a simple 2-ploidy scenario (should behave like diploid spectral)."""
    # Use a simple diploid-like scenario for spectral (K=2) to validate clustering:
    # True haplotypes: [0,1,1,0], [1,0,0,1] (exact opposites for clarity)
    hap0 = np.array([0,1,1,0], dtype=int)
    hap1 = np.array([1,0,0,1], dtype=int)
    # Create 4 reads, two fully matching hap0, two matching hap1
    reads_matrix = np.array([hap0, hap0, hap1, hap1], dtype=int)
    npz_file = os.path.join(tmp_path, "spec2_input.npz")
    np.savez_compressed(npz_file, reads=reads_matrix)
    out_prefix = os.path.join(tmp_path, "spec_output")
    # Run spectral clustering with k=2
    args = type("Args", (), {})()
    args.input = npz_file
    args.ploidy = 2
    args.output_prefix = out_prefix
    args.min_overlap = 1
    args.seed = 42
    poly_spec.main(args)
    # Load outputs
    hap_file = out_prefix + ".haplotypes.tsv"
    assign_file = out_prefix + ".assignments.tsv"
    assert os.path.exists(hap_file), "Spectral haplotypes output missing."
    haps_out = [line.strip() for line in open(hap_file) if line.strip()]
    # Parse haplotypes from TSV
    haps_out = []
    for line in open(hap_file):
        if line.strip():
            parts = line.strip().split("\t")
            assert len(parts) == 2, f"Unexpected haplotype format: {line}"
            haps_out.append(parts[1])
    # Expect 2 haplotypes of length 4
    assert len(haps_out) == 2 and all(len(h) == 4 for h in haps_out), "Spectral output format incorrect."
    # They should correspond to the two true haplotypes (or their complements, since consensus might pick the opposite allele for each cluster)
    out_set = set(haps_out)
    # Either we get the exact true haplotypes, or possibly the bitwise complements ("1001" vs "0110") depending on majority definition.
    true_set1 = {"0110", "1001"}  # corresponds to hap0="0110", hap1="1001"
    true_set2 = {"1001", "0110"}  # just order swap (set covers it anyway)
    complement_set = {"..."}  # not needed actually, above covers both possibilities in set form
    assert out_set == true_set1 or out_set == true_set2, f"Spectral clustering haplotypes incorrect: {out_set}"
    # Assignments: expect two reads in cluster0, two in cluster1
    assigns = [line.strip() for line in open(assign_file) if line.strip()]
    assert assigns.count("0") == 2 and assigns.count("1") == 2, "Spectral should cluster reads into two even groups for this data."
