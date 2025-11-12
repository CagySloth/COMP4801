#!/usr/bin/env python3
"""
phase_io.py

Shared I/O and utilities for simple phasing baselines.

Supports loading NPZ files produced by read_reads_tsv.py:
- dense NPZ keys: mode="dense", read_ids, hap_truth, alleles (int8, R x N)
- sparse NPZ keys: mode="sparse", read_ids, hap_truth, read_ptr, positions, alleles, num_variants
and converts to a dense alleles matrix with -1 for missing.

Also includes:
- write_haplotypes_tsv(haplotypes, path)
- write_assignments_tsv(read_ids, assignments, path)
- compute_mec(alleles, haplotypes, assignments)
- hap_truth_accuracy(hap_truth, assignments) with best label mapping
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np


@dataclass
class ReadsData:
    alleles: np.ndarray      # (R, N) int8, -1 for missing
    read_ids: List[str]      # length R
    hap_truth: Optional[np.ndarray]  # (R,) int32 or None
    R: int
    N: int


def _load_npz_dense(npz) -> ReadsData:
    alleles = np.array(npz["alleles"], dtype=np.int8)
    read_ids = npz["read_ids"].astype(str).tolist() if "read_ids" in npz else [f"r{i}" for i in range(alleles.shape[0])]
    hap_truth = np.array(npz["hap_truth"], dtype=np.int32) if "hap_truth" in npz else None
    R, N = alleles.shape
    return ReadsData(alleles=alleles, read_ids=read_ids, hap_truth=hap_truth, R=R, N=N)


def _load_npz_sparse(npz) -> ReadsData:
    read_ids = npz["read_ids"].astype(str).tolist() if "read_ids" in npz else None
    hap_truth = np.array(npz["hap_truth"], dtype=np.int32) if "hap_truth" in npz else None
    read_ptr = np.array(npz["read_ptr"], dtype=np.int32)
    positions = np.array(npz["positions"], dtype=np.int32)
    alle = np.array(npz["alleles"], dtype=np.int8)
    num_variants = int(np.array(npz["num_variants"]).item())
    R = read_ptr.size - 1
    N = num_variants
    alleles = np.full((R, N), -1, dtype=np.int8)
    for i in range(R):
        s, e = read_ptr[i], read_ptr[i + 1]
        if e > s:
            pos = positions[s:e]
            val = alle[s:e]
            alleles[i, pos] = val
    if read_ids is None:
        read_ids = [f"r{i}" for i in range(R)]
    return ReadsData(alleles=alleles, read_ids=read_ids, hap_truth=hap_truth, R=R, N=N)


def load_reads(input_path: str) -> ReadsData:
    """
    Load reads from NPZ produced by read_reads_tsv.py (dense or sparse) and return dense alleles.
    """
    if not input_path.lower().endswith(".npz"):
        raise ValueError("Please provide an NPZ produced by read_reads_tsv.py (--to-npz).")
    npz = np.load(input_path, allow_pickle=False)
    mode = None
    if "mode" in npz:
        try:
            mode = str(npz["mode"].item())
        except Exception:
            # Some NumPy versions may return a 0-d array of dtype='<U'
            mode = str(npz["mode"])
    if mode == "dense":
        return _load_npz_dense(npz)
    elif mode == "sparse":
        return _load_npz_sparse(npz)
    else:
        # Fallback: try dense keys
        if "alleles" in npz:
            return _load_npz_dense(npz)
        raise ValueError("Unrecognized NPZ format. Ensure you used read_reads_tsv.py --to-npz.")


def write_haplotypes_tsv(haplotypes: np.ndarray, path: str) -> None:
    """
    haplotypes: (K, N) uint8 with 0/1
    """
    K, N = haplotypes.shape
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"# haplotypes\t{K}\n")
        f.write(f"# variants\t{N}\n")
        f.write("# format: hap_id<TAB>alleles_as_0_1_string\n")
        for k in range(K):
            bitstr = "".join("1" if x else "0" for x in haplotypes[k])
            f.write(f"h{k}\t{bitstr}\n")


def write_assignments_tsv(read_ids: List[str], assignments: np.ndarray, path: str) -> None:
    """
    assignments: (R,) integers in [0, K-1]
    """
    with open(path, "w", encoding="utf-8") as f:
        f.write("# read_id\tpredicted_hap\n")
        for rid, a in zip(read_ids, assignments.tolist()):
            f.write(f"{rid}\t{int(a)}\n")


def compute_mec(alleles: np.ndarray, haplotypes: np.ndarray, assignments: np.ndarray) -> Tuple[int, np.ndarray]:
    """
    Minimum Error Correction (given assignments):
    - Sum over reads of mismatches to assigned haplotype (ignore missing).
    Returns (total_mec, per_read_mec).
    """
    R, N = alleles.shape
    K, N2 = haplotypes.shape
    if N != N2:
        raise ValueError("alleles and haplotypes must have the same number of variants.")
    per_read = np.zeros(R, dtype=np.int64)
    valid = alleles >= 0
    # Loop over clusters to avoid building a huge (R,K,N) tensor
    for k in range(K):
        idx = np.where(assignments == k)[0]
        if idx.size == 0:
            continue
        a = alleles[idx]       # (r_k, N)
        v = valid[idx]
        h = haplotypes[k][None, :]  # (1, N)
        per_read[idx] = np.sum(v & (a != h), axis=1)
    return int(per_read.sum()), per_read


def global_majority(alleles: np.ndarray) -> np.ndarray:
    """
    Return per-site majority allele across all reads (ignoring missing), shape (N,), uint8.
    Ties break to 0. Sites with no coverage default to 0.
    """
    N = alleles.shape[1]
    maj = np.zeros(N, dtype=np.uint8)
    valid = alleles >= 0
    ones = (alleles == 1)
    ones_count = ones.sum(axis=0)
    obs_count = valid.sum(axis=0)
    zeros_count = obs_count - ones_count
    maj = (ones_count > zeros_count).astype(np.uint8)
    # ties => 0 (already)
    return maj


def hap_truth_accuracy(hap_truth: Optional[np.ndarray], assignments: np.ndarray) -> Optional[dict]:
    """
    Compute assignment accuracy if hap_truth is available (>=0).
    Finds best label mapping. If K<=8, uses exhaustive permutation; else greedy.
    Returns dict with keys: accuracy, mapping (pred->true), confusion (list of lists).
    """
    if hap_truth is None:
        return None
    mask = hap_truth >= 0
    if not np.any(mask):
        return None
    ht = hap_truth[mask]
    pr = assignments[mask]
    K_true = int(ht.max()) + 1
    K_pred = int(assignments.max()) + 1
    # Build confusion matrix C[true, pred]
    C = np.zeros((K_true, K_pred), dtype=np.int64)
    for t, p in zip(ht.tolist(), pr.tolist()):
        if 0 <= t < K_true and 0 <= p < K_pred:
            C[t, p] += 1

    # If sizes differ, pad to square
    K = max(K_true, K_pred)
    Cpad = np.zeros((K, K), dtype=np.int64)
    Cpad[:K_true, :K_pred] = C

    mapping = list(range(K))  # pred -> true
    best = -1

    if K <= 8:
        import itertools
        for perm in itertools.permutations(range(K)):
            total = sum(Cpad[perm[p], p] for p in range(K))
            if total > best:
                best = total
                mapping = list(perm)
    else:
        # Greedy mapping
        used_true = set()
        mapping = [-1] * K
        for p in range(K):
            t_best = int(np.argmax(Cpad[:, p]))
            if t_best in used_true:
                # pick next best unused
                sorted_idx = np.argsort(-Cpad[:, p]).tolist()
                for t in sorted_idx:
                    if t not in used_true:
                        t_best = int(t)
                        break
            mapping[p] = t_best
            used_true.add(t_best)
        best = sum(Cpad[mapping[p], p] for p in range(K))

    total_used = int(mask.sum())
    acc = best / total_used if total_used > 0 else 0.0
    return {
        "accuracy": float(acc),
        "mapping_pred_to_true": mapping,
        "confusion": C.tolist(),
    }


def write_summary_json(obj: dict, path: str) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)