#!/usr/bin/env python3
"""
polyploid_em_cluster.py

Polyploid phasing via hard EM with K haplotypes:
- Initialize K haplotypes using k-means++-like seeding from reads.
- Repeat: assign reads to nearest hap; update each hap by majority over its assigned reads.
- Stop when assignments stabilize or max iterations reached.

Input:
  NPZ from read_reads_tsv.py (dense or sparse)

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
    compute_mec,
    hap_truth_accuracy,
    global_majority,
)


def dist_reads_to_hap(alleles: np.ndarray, hap: np.ndarray) -> np.ndarray:
    """
    Hamming distance ignoring missing, per read to given hap.
    """
    valid = alleles >= 0
    return np.sum(valid & (alleles != hap[None, :]), axis=1)


def kmeanspp_init_haps(alleles: np.ndarray, K: int, rng: np.random.Generator) -> np.ndarray:
    """
    k-means++-like seeding from reads, then fill missing with global majority.
    """
    R, N = alleles.shape
    gmaj = global_majority(alleles)

    # Choose first center uniformly
    centers = []
    first = int(rng.integers(0, R))
    centers.append(alleles[first].copy())

    # Distance function (overlap-normalized)
    def norm_dist(a, b) -> float:
        v = (a >= 0) & (b >= 0)
        ov = int(v.sum())
        if ov == 0:
            return 1.0
        return float(np.sum(a[v] != b[v]) / ov)

    # Choose remaining centers
    d2 = np.array([norm_dist(alleles[i], centers[0]) ** 2 for i in range(R)], dtype=np.float64)
    for _ in range(1, K):
        probs = d2 / (d2.sum() if d2.sum() > 0 else 1.0)
        i = int(rng.choice(R, p=probs))
        centers.append(alleles[i].copy())
        # Update d2
        for r in range(R):
            d = norm_dist(alleles[r], centers[-1]) ** 2
            if d < d2[r]:
                d2[r] = d

    # Convert centers with missing to haplotypes by filling with global majority
    H = np.zeros((K, N), dtype=np.uint8)
    for k in range(K):
        c = centers[k]
        m = c < 0
        if np.any(m):
            c = c.copy()
            c[m] = gmaj[m]
        H[k] = c.astype(np.uint8)
    return H


def update_haplotypes(alleles: np.ndarray, assignments: np.ndarray, H_prev: np.ndarray) -> np.ndarray:
    R, N = alleles.shape
    K = H_prev.shape[0]
    H = H_prev.copy()
    for k in range(K):
        idx = np.where(assignments == k)[0]
        if idx.size == 0:
            # Keep previous hap; optionally reseed from a far read
            continue
        A = alleles[idx]
        valid = A >= 0
        ones = (A == 1)
        ones_count = ones.sum(axis=0)
        obs_count = valid.sum(axis=0)
        zeros_count = obs_count - ones_count
        h = (ones_count > zeros_count).astype(np.uint8)  # ties->0
        no_obs = obs_count == 0
        h[no_obs] = H_prev[k, no_obs]
        H[k] = h
    return H


def main():
    ap = argparse.ArgumentParser(description="Polyploid phasing via hard EM clustering (K haplotypes).")
    ap.add_argument("-i", "--input", required=True, help="NPZ from read_reads_tsv.py")
    ap.add_argument("-k", "--ploidy", type=int, required=True, help="Number of haplotypes (K)")
    ap.add_argument("-o", "--output-prefix", required=True, help="Output prefix")
    ap.add_argument("--max-iters", type=int, default=40, help="Max iterations (default: 40)")
    ap.add_argument("--tol-iters", type=int, default=3, help="Stop if assignments unchanged for T iterations (default: 3)")
    ap.add_argument("-s", "--seed", type=int, default=None, help="Random seed")
    args = ap.parse_args()

    data = load_reads(args.input)
    A = data.alleles
    K = args.ploidy
    rng = np.random.default_rng(args.seed)

    # Initialize haplotypes
    H = kmeanspp_init_haps(A, K, rng)

    prev_assign = None
    stable = 0
    it = 0

    while it < args.max_iters:
        it += 1
        # Assign
        dists = np.stack([dist_reads_to_hap(A, H[k]) for k in range(K)], axis=1)  # (R, K)
        assign = np.argmin(dists, axis=1).astype(np.int32)

        if prev_assign is not None and np.array_equal(assign, prev_assign):
            stable += 1
        else:
            stable = 0
        prev_assign = assign.copy()

        # Update
        H = update_haplotypes(A, assign, H)

        if stable >= args.tol_iters:
            break

    mec, per_read = compute_mec(A, H, prev_assign)
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
        "algorithm": "polyploid_hard_em",
        "K": K,
        "iterations": it,
        "MEC_total": int(mec),
        "MEC_mean_per_read": float(np.mean(per_read)) if per_read.size else 0.0,
        "cluster_sizes": {int(k): int(np.sum(prev_assign == k)) for k in range(K)},
        "assignment_accuracy": acc_info["accuracy"] if acc_info else None,
        "label_mapping_pred_to_true": acc_info["mapping_pred_to_true"] if acc_info else None,
    }
    write_summary_json(summary, sum_path)
    print(f"Wrote:\n  {hap_path}\n  {asg_path}\n  {sum_path}")


if __name__ == "__main__":
    main()