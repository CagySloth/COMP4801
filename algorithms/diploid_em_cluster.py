#!/usr/bin/env python3
"""
diploid_em_cluster.py

Simple diploid phasing via hard EM clustering:
- Initialize two haplotypes.
- Repeat: assign reads to nearest haplotype; update haplotypes by per-site majority.
- Stop when assignments stabilize or after max iterations.

Input:
  NPZ from read_reads_tsv.py (dense or sparse) via --to-npz/--to-dense-npz.

Outputs:
  <prefix>.haplotypes.tsv
  <prefix>.assignments.tsv
  <prefix>.summary.json
"""

from __future__ import annotations

import argparse
import os
from typing import Tuple

import numpy as np

from phase_io import (
    load_reads,
    write_haplotypes_tsv,
    write_assignments_tsv,
    write_summary_json,
    global_majority,
    compute_mec,
    hap_truth_accuracy,
)


def init_haplotypes_kpp(alleles: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """
    k-means++-like initialization for K=2 using read rows as seeds, then consensus per seed.
    """
    R, N = alleles.shape
    valid = alleles >= 0

    # Pick first seed uniformly
    i0 = int(rng.integers(0, R))
    seed0 = alleles[i0].copy()

    # Pick second seed: farthest by overlap-normalized Hamming (ignoring missing)
    best_i = i0
    best_d = -1.0
    for i in range(R):
        v = valid[i] & valid[i0]
        ov = int(v.sum())
        if ov == 0:
            continue
        d = np.sum((alleles[i][v] != alleles[i0][v])) / ov
        if d > best_d:
            best_d = d
            best_i = i

    seed1 = alleles[best_i].copy() if best_d >= 0 else (1 - seed0.clip(min=0))  # fallback

    # Convert seeds into 0/1 haplotypes; fill missing with global majority
    gmaj = global_majority(alleles)
    def to_hap(s):
        h = s.copy()
        m = h < 0
        if np.any(m):
            h[m] = gmaj[m]
        return h.astype(np.uint8)

    H = np.stack([to_hap(seed0), to_hap(seed1)], axis=0)
    return H


def assign_reads(alleles: np.ndarray, H: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Assign reads to the nearest haplotype by mismatch count (ignoring missing).
    Returns (assignments, dist_to_assigned).
    """
    R, N = alleles.shape
    K = H.shape[0]
    valid = alleles >= 0
    dists = np.zeros((R, K), dtype=np.int64)
    for k in range(K):
        mismatch = valid & (alleles != H[k][None, :])
        dists[:, k] = mismatch.sum(axis=1)
    assignments = np.argmin(dists, axis=1).astype(np.int32)
    return assignments, dists[np.arange(R), assignments]


def update_haplotypes(alleles: np.ndarray, assignments: np.ndarray, H_prev: np.ndarray) -> np.ndarray:
    """
    Update each haplotype per site by majority of assigned reads (ignoring missing).
    Empty cluster: keep previous haplotype.
    """
    R, N = alleles.shape
    K = H_prev.shape[0]
    H = H_prev.copy()
    for k in range(K):
        idx = np.where(assignments == k)[0]
        if idx.size == 0:
            continue
        A = alleles[idx]  # (r_k, N)
        valid = A >= 0
        ones = (A == 1)
        ones_count = ones.sum(axis=0)
        obs_count = valid.sum(axis=0)
        zeros_count = obs_count - ones_count
        h = (ones_count > zeros_count).astype(np.uint8)  # ties -> 0
        # If a site has no observations in this cluster, keep previous
        no_obs = obs_count == 0
        h[no_obs] = H_prev[k, no_obs]
        H[k] = h
    return H


def main():
    ap = argparse.ArgumentParser(description="Diploid phasing via hard EM clustering.")
    ap.add_argument("-i", "--input", required=True, help="NPZ from read_reads_tsv.py")
    ap.add_argument("-o", "--output-prefix", required=True, help="Output prefix")
    ap.add_argument("--max-iters", type=int, default=30, help="Max iterations (default: 30)")
    ap.add_argument("--tol-iters", type=int, default=2, help="Stop if no assignment change for T iterations (default: 2)")
    ap.add_argument("-s", "--seed", type=int, default=None, help="Random seed")
    args = ap.parse_args()

    data = load_reads(args.input)
    alleles = data.alleles
    rng = np.random.default_rng(args.seed)
    if data.R == 0:
        # No reads: output two haplotypes (all 0s) and empty assignments
        H = np.zeros((2, data.N), dtype=np.uint8)
        assign = np.zeros(0, dtype=np.int32)
        mec, per_read = 0, np.zeros(0, dtype=np.int64)
        acc_info = None
        # Proceed to output generation
    else:
        # Initialize haplotypes
        H = init_haplotypes_kpp(alleles, rng)

        prev_assign = None
        stable_count = 0

        for it in range(1, args.max_iters + 1):
            assign, _ = assign_reads(alleles, H)
            if prev_assign is not None and np.array_equal(assign, prev_assign):
                stable_count += 1
            else:
                stable_count = 0
            prev_assign = assign.copy()

            H_new = update_haplotypes(alleles, assign, H)
            H = H_new

            if stable_count >= args.tol_iters:
                break

        # Final evaluation
        mec, per_read = compute_mec(alleles, H, prev_assign)
        acc_info = hap_truth_accuracy(data.hap_truth, prev_assign)

    # Outputs
    outdir = os.path.dirname(os.path.abspath(args.output_prefix))
    os.makedirs(outdir, exist_ok=True)
    hap_path = f"{args.output_prefix}.haplotypes.tsv"
    asg_path = f"{args.output_prefix}.assignments.tsv"
    sum_path = f"{args.output_prefix}.summary.json"
    write_haplotypes_tsv(H, hap_path)
    write_assignments_tsv(data.read_ids, prev_assign, asg_path)

    summary = {
        "R": data.R,
        "N": data.N,
        "algorithm": "diploid_hard_em",
        "iterations": it,
        "MEC_total": int(mec),
        "MEC_mean_per_read": float(np.mean(per_read)) if per_read.size else 0.0,
        "cluster_sizes": {int(k): int(np.sum(prev_assign == k)) for k in range(2)},
        "assignment_accuracy": acc_info["accuracy"] if acc_info else None,
        "label_mapping_pred_to_true": acc_info["mapping_pred_to_true"] if acc_info else None,
    }
    write_summary_json(summary, sum_path)
    print(f"Wrote:\n  {hap_path}\n  {asg_path}\n  {sum_path}")


if __name__ == "__main__":
    main()