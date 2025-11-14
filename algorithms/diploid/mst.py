#!/usr/bin/env python3
"""
diploid_mst_phase.py

Diploid phasing via SNP–SNP co-occurrence and maximum spanning forest:
- Genotype sites by majority across reads (homozygous vs heterozygous).
- For heterozygous sites, compute pairwise co-occurrence weights:
    w = (#agree) - (#disagree) over overlapping reads.
- Build a maximum spanning forest on |w| (Kruskal), storing edge parity (equal/different).
- Pick a root per component; set its allele by local majority; propagate along edges.
- Assign each read to the nearer haplotype; report MEC.

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
from typing import List, Tuple

import numpy as np

from algorithms.io.parser import parse_sparse_tsv, load_reads
from algorithms.io.writer import write_haplotypes_tsv, write_summary_json, write_assignments_tsv
from algorithms.io.reads_data import ReadsData
from algorithms.eval.metrics import compute_mec, hap_truth_accuracy


class DSU:
    def __init__(self, n: int):
        self.p = list(range(n))
        self.r = [0] * n

    def find(self, x: int) -> int:
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a: int, b: int) -> bool:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        if self.r[ra] < self.r[rb]:
            ra, rb = rb, ra
        self.p[rb] = ra
        if self.r[ra] == self.r[rb]:
            self.r[ra] += 1
        return True


def genotype_sites(alleles: np.ndarray, min_minor: int = 1) -> Tuple[np.ndarray, np.ndarray]:
    """
    Determine per-site genotype from reads.
    Returns:
      hom_mask: N bool (True if homozygous)
      major:    N uint8 majority allele (0/1), ties->0
    het_mask = ~hom_mask
    """
    valid = alleles >= 0
    ones = (alleles == 1)
    ones_count = ones.sum(axis=0)
    obs_count = valid.sum(axis=0)
    zeros_count = obs_count - ones_count
    major = (ones_count > zeros_count).astype(np.uint8)
    minor_count = np.minimum(ones_count, zeros_count)
    het_mask = minor_count >= min_minor
    hom_mask = ~het_mask
    return hom_mask, major


def build_edges(alleles: np.ndarray, het_idx: np.ndarray, min_overlap: int) -> List[Tuple[int, int, int, int]]:
    """
    Build edges among heterozygous sites.
    Returns list of (abs_w, u_idx, v_idx, parity) with parity=0 (equal) or 1 (different).
    u_idx, v_idx are indices into het_idx array (0..M-1).
    """
    A = alleles
    R = A.shape[0]
    edges = []
    M = het_idx.size
    for ui in range(M):
        j = het_idx[ui]
        col_j = A[:, j]
        for vi in range(ui + 1, M):
            k = het_idx[vi]
            col_k = A[:, k]
            vmask = (col_j >= 0) & (col_k >= 0)
            ov = int(np.sum(vmask))
            if ov < min_overlap:
                continue
            eq = int(np.sum((col_j[vmask] == col_k[vmask])))
            total = ov
            w = 2 * eq - total  # (#agree) - (#disagree)
            if w == 0:
                continue
            parity = 0 if w > 0 else 1
            edges.append((abs(w), ui, vi, parity))
    # sort desc by abs weight
    edges.sort(key=lambda x: (-x[0], x[1], x[2]))
    return edges


def mst_forest(M: int, edges: List[Tuple[int, int, int, int]]) -> List[List[Tuple[int, int, int]]]:
    """
    Build maximum spanning forest using Kruskal.
    Returns adjacency list: adj[u] = list of (v, parity), using het-indexed nodes 0..M-1.
    """
    dsu = DSU(M)
    adj = [[] for _ in range(M)]
    for wabs, u, v, parity in edges:
        if dsu.union(u, v):
            # Add undirected edge with parity
            adj[u].append((v, parity))
            adj[v].append((u, parity))
    return adj


def orient_component(adj: List[List[Tuple[int, int]]],
                     comp_nodes: List[int],
                     root: int,
                     root_value: int) -> np.ndarray:
    """
    BFS from root, propagate parity: parity=0 => same as parent; 1 => flip.
    Returns hap0 vector for nodes in comp (indexed by absolute SNP indices), but here we work on comp node indices.
    """
    H_local = {}
    from collections import deque
    q = deque([root])
    H_local[root] = int(root_value)
    while q:
        u = q.popleft()
        for v, parity in adj[u]:
            if v in H_local:
                continue
            if parity == 0:
                H_local[v] = H_local[u]
            else:
                H_local[v] = 1 - H_local[u]
            q.append(v)
    # Convert to dense array aligned with comp_nodes
    out = np.zeros(len(adj), dtype=np.uint8)  # we will only use entries for comp_nodes
    for u, val in H_local.items():
        out[u] = val
    return out


def main():
    ap = argparse.ArgumentParser(description="Diploid MST-based phasing.")
    ap.add_argument("-i", "--input", required=True, help="NPZ from read_reads_tsv.py")
    ap.add_argument("-o", "--output-prefix", required=True, help="Output prefix")
    ap.add_argument("--min-overlap", type=int, default=3, help="Min overlapping reads for a pair (default: 3)")
    ap.add_argument("--min-het-minor", type=int, default=1, help="Min minor count to call a site heterozygous (default: 1)")
    args = ap.parse_args()

    data = load_reads(args.input)
    A = data.alleles
    R, N = A.shape

    hom_mask, major = genotype_sites(A, min_minor=args.min_het_minor)
    het_mask = ~hom_mask
    het_idx = np.where(het_mask)[0]
    M = het_idx.size

    # Initialize haplotypes
    H = np.zeros((2, N), dtype=np.uint8)
    # Homozygous: both haplotypes = majority allele
    H[:, hom_mask] = major[hom_mask][None, :]

    if M == 0:
        # No heterozygous sites; trivial assignment
        assign = np.zeros(R, dtype=np.int32)
        mec, per_read = compute_mec(A, H, assign)
    else:
        # Build edges and forest
        edges = build_edges(A, het_idx, min_overlap=args.min_overlap)
        adj = mst_forest(M, edges)

        # Root per component and orient
        visited = np.zeros(M, dtype=bool)
        gmaj = global_majority(A)

        for u in range(M):
            if visited[u]:
                continue
            # BFS to collect component
            comp = []
            stack = [u]
            visited[u] = True
            while stack:
                x = stack.pop()
                comp.append(x)
                for v, _ in adj[x]:
                    if not visited[v]:
                        visited[v] = True
                        stack.append(v)
            # Pick root allele from global majority at that site
            root_node = u
            abs_root = het_idx[root_node]
            root_val = int(gmaj[abs_root])
            # Orient
            hap0_local = orient_component(adj, comp, root_node, root_val)
            # Fill H
            for node in comp:
                abs_idx = het_idx[node]
                val = hap0_local[node]
                H[0, abs_idx] = val
                H[1, abs_idx] = 1 - val

        # Assign reads to hap0/hap1
        valid = A >= 0
        d0 = np.sum(valid & (A != H[0][None, :]), axis=1)
        d1 = np.sum(valid & (A != H[1][None, :]), axis=1)
        assign = (d1 < d0).astype(np.int32)  # tie -> hap0

        mec, per_read = compute_mec(A, H, assign)

    acc_info = hap_truth_accuracy(data.hap_truth, assign)

    # Outputs
    outdir = os.path.dirname(os.path.abspath(args.output_prefix))
    os.makedirs(outdir, exist_ok=True)
    hap_path = f"{args.output_prefix}.haplotypes.tsv"
    asg_path = f"{args.output_prefix}.assignments.tsv"
    sum_path = f"{args.output_prefix}.summary.json"
    write_haplotypes_tsv(hap_path, H)
    write_assignments_tsv(data.read_ids, assign, asg_path)

    summary = {
        "R": data.R,
        "N": data.N,
        "algorithm": "diploid_mst",
        "MEC_total": int(mec),
        "MEC_mean_per_read": float(np.mean(per_read)) if per_read.size else 0.0,
        "heterozygous_sites": int(M),
        "min_overlap": args.min_overlap,
        "cluster_sizes": {int(k): int(np.sum(assign == k)) for k in range(2)},
        "assignment_accuracy": acc_info["accuracy"] if acc_info else None,
        "label_mapping_pred_to_true": acc_info["mapping_pred_to_true"] if acc_info else None,
    }
    write_summary_json(summary, sum_path)
    print(f"Wrote:\n  {hap_path}\n  {asg_path}\n  {sum_path}")


if __name__ == "__main__":
    main()