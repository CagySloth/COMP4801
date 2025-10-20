#!/usr/bin/env python3
"""
polyploid_spectral_cluster.py

Polyploid phasing via spectral clustering of reads:
- Build read agreement graph W_ij = (#agree / #overlap) for read pairs with overlap >= min_overlap.
- Compute normalized Laplacian L = I - D^{-1/2} W D^{-1/2} and its k smallest eigenvectors.
- Cluster rows of eigenvector embedding with k-means; derive consensus haplotypes per cluster.

Input:
  NPZ from read_reads_tsv.py (dense or sparse)

Outputs:
  <prefix>.haplotypes.tsv
  <prefix>.assignments.tsv
  <prefix>.summary.json

Caveats: O(R^2 * N) to build W; intended for small/moderate R.
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
)


def build_agreement_matrix(A: np.ndarray, min_overlap: int) -> np.ndarray:
    """
    W_ij = (#agree / #overlap) for overlap >= min_overlap; else 0.
    Symmetric, with zeros on diagonal.
    """
    R, N = A.shape
    W = np.zeros((R, R), dtype=np.float64)
    for i in range(R):
        ai = A[i]
        for j in range(i + 1, R):
            aj = A[j]
            v = (ai >= 0) & (aj >= 0)
            ov = int(np.sum(v))
            if ov < min_overlap:
                continue
            agree = int(np.sum(ai[v] == aj[v]))
            W[i, j] = agree / ov
            W[j, i] = W[i, j]
    return W


def spectral_kmeans(W: np.ndarray, K: int, rng: np.random.Generator, iters: int = 50) -> Tuple[np.ndarray, np.ndarray]:
    """
    Return assignments (R,) to K clusters using normalized spectral clustering + k-means.
    Also returns eigenvector embedding Y (R, K).
    """
    R = W.shape[0]
    d = np.sum(W, axis=1)
    # Handle zero-degree nodes
    d_safe = np.where(d > 0, d, 1.0)
    Dm12 = 1.0 / np.sqrt(d_safe)
    Dm12_mat = np.diag(Dm12)
    L = np.eye(R) - Dm12_mat @ W @ Dm12_mat

    # Compute k smallest eigenvectors of L
    vals, vecs = np.linalg.eigh(L)  # full eigen-decomp (OK for small R)
    idx = np.argsort(vals)[:K]
    Y = vecs[:, idx]
    # Row-normalize Y
    row_norm = np.linalg.norm(Y, axis=1, keepdims=True)
    row_norm[row_norm == 0] = 1.0
    Y = Y / row_norm

    # k-means on Y
    assign = kmeans(Y, K, rng, iters=iters)
    return assign, Y


def kmeans(X: np.ndarray, K: int, rng: np.random.Generator, iters: int = 50) -> np.ndarray:
    """
    Simple k-means with k-means++ init.
    """
    R, D = X.shape
    # Init
    centers = np.empty((K, D), dtype=np.float64)
    i0 = int(rng.integers(0, R))
    centers[0] = X[i0]
    # distances to nearest center
    dist2 = np.sum((X - centers[0]) ** 2, axis=1)
    for k in range(1, K):
        probs = dist2 / (dist2.sum() if dist2.sum() > 0 else 1.0)
        i = int(rng.choice(R, p=probs))
        centers[k] = X[i]
        dist2 = np.minimum(dist2, np.sum((X - centers[k]) ** 2, axis=1))

    # Lloyd iterations
    assign = np.zeros(R, dtype=np.int32)
    for _ in range(iters):
        # Assign
        d2 = np.sum((X[:, None, :] - centers[None, :, :]) ** 2, axis=2)
        new_assign = np.argmin(d2, axis=1)
        if np.array_equal(new_assign, assign):
            break
        assign = new_assign
        # Update
        for k in range(K):
            idx = np.where(assign == k)[0]
            if idx.size > 0:
                centers[k] = X[idx].mean(axis=0)
            else:
                # Reseed empty center
                centers[k] = X[int(rng.integers(0, R))]
    return assign


def consensus_haplotypes(A: np.ndarray, assign: np.ndarray, K: int) -> np.ndarray:
    """
    Per-cluster per-site majority; ties->0; if no obs, fill with 0.
    """
    R, N = A.shape
    H = np.zeros((K, N), dtype=np.uint8)
    for k in range(K):
        idx = np.where(assign == k)[0]
        if idx.size == 0:
            continue
        Ak = A[idx]
        valid = Ak >= 0
        ones = (Ak == 1)
        ones_count = ones.sum(axis=0)
        obs_count = valid.sum(axis=0)
        zeros_count = obs_count - ones_count
        H[k] = (ones_count > zeros_count).astype(np.uint8)  # ties->0
    return H


def main():
    ap = argparse.ArgumentParser(description="Polyploid spectral clustering phasing.")
    ap.add_argument("-i", "--input", required=True, help="NPZ from read_reads_tsv.py")
    ap.add_argument("-k", "--ploidy", type=int, required=True, help="Number of haplotypes (K)")
    ap.add_argument("-o", "--output-prefix", required=True, help="Output prefix")
    ap.add_argument("--min-overlap", type=int, default=3, help="Min overlapping positions to connect reads (default: 3)")
    ap.add_argument("-s", "--seed", type=int, default=None, help="Random seed")
    args = ap.parse_args()

    data = load_reads(args.input)
    A = data.alleles
    K = args.ploidy
    rng = np.random.default_rng(args.seed)

    # Build agreement matrix
    W = build_agreement_matrix(A, min_overlap=args.min_overlap)
    # Cluster
    assign, _ = spectral_kmeans(W, K=K, rng=rng, iters=60)
    # Consensus haplotypes
    H = consensus_haplotypes(A, assign, K)
    # MEC
    mec, per_read = compute_mec(A, H, assign)
    acc_info = hap_truth_accuracy(data.hap_truth, assign)

    # Outputs
    outdir = os.path.dirname(os.path.abspath(args.output_prefix))
    os.makedirs(outdir, exist_ok=True)
    hap_path = f"{args.output_prefix}.haplotypes.tsv"
    asg_path = f"{args.output_prefix}.assignments.tsv"
    sum_path = f"{args.output_prefix}.summary.json"
    write_haplotypes_tsv(H, hap_path)
    write_assignments_tsv(data.read_ids, assign, asg_path)

    summary = {
        "R": data.R,
        "N": data.N,
        "algorithm": "polyploid_spectral",
        "K": K,
        "min_overlap": args.min_overlap,
        "MEC_total": int(mec),
        "MEC_mean_per_read": float(np.mean(per_read)) if per_read.size else 0.0,
        "cluster_sizes": {int(k): int(np.sum(assign == k)) for k in range(K)},
        "assignment_accuracy": acc_info["accuracy"] if acc_info else None,
        "label_mapping_pred_to_true": acc_info["mapping_pred_to_true"] if acc_info else None,
    }
    write_summary_json(summary, sum_path)
    print(f"Wrote:\n  {hap_path}\n  {asg_path}\n  {sum_path}")


if __name__ == "__main__":
    main()