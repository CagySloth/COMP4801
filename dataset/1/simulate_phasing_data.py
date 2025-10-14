#!/usr/bin/env python3
"""
simulate_phasing_data.py

Generate synthetic data for benchmarking phasing algorithms.

Key parameters:
- ploidy (number of haplotypes)
- per-site error rate
- per-site missing rate
- read length (in number of variant sites)
- number of reads

Outputs:
- <prefix>.haplotypes.tsv       Ground truth haplotypes (P rows, N columns as a 0/1 string per haplotype)
- <prefix>.haplotypes.npz       Same haplotypes in a compact NumPy format (array shape: [P, N], dtype=uint8)
- <prefix>.reads.sparse.tsv     One read per row; columns: read_id, hap_truth, start, positions, alleles
- <prefix>.summary.json         Parameters and simple summary statistics

Notes:
- Read length is measured in number of variant sites (not base pairs).
- Missing is applied first, then (independently) substitution error is applied to non-missing bases.
- Alleles are biallelic 0/1; missing is denoted with '-'.
- Positions are 0-based variant indices.
"""

from __future__ import annotations

import argparse
import json
import os
from dataclasses import dataclass, asdict
from typing import Tuple

import numpy as np


@dataclass
class SimParams:
    ploidy: int
    num_variants: int
    num_reads: int
    read_length: int
    error_rate: float
    missing_rate: float
    maf_alpha: float
    maf_beta: float
    ensure_polymorphic: bool
    seed: int | None
    output_prefix: str


@dataclass
class SimSummary:
    ploidy: int
    num_variants: int
    num_reads: int
    read_length: int
    error_rate: float
    missing_rate: float
    maf_alpha: float
    maf_beta: float
    ensure_polymorphic: bool
    seed: int | None
    expected_coverage: float
    mean_coverage: float
    min_coverage: int
    max_coverage: int
    observed_error_rate_over_non_missing: float
    observed_missing_fraction_over_all_bases: float


def _validate_params(p: SimParams) -> None:
    if p.ploidy < 1:
        raise ValueError("ploidy must be >= 1")
    if p.num_variants < 1:
        raise ValueError("num_variants must be >= 1")
    if p.num_reads < 0:
        raise ValueError("num_reads must be >= 0")
    if p.read_length < 1:
        raise ValueError("read_length must be >= 1")
    if not (0.0 <= p.error_rate <= 1.0):
        raise ValueError("error_rate must be in [0, 1]")
    if not (0.0 <= p.missing_rate <= 1.0):
        raise ValueError("missing_rate must be in [0, 1]")
    if p.maf_alpha <= 0 or p.maf_beta <= 0:
        raise ValueError("maf_alpha and maf_beta must be > 0")


def simulate_haplotypes(
    ploidy: int,
    num_variants: int,
    maf_alpha: float,
    maf_beta: float,
    ensure_polymorphic: bool,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Simulate biallelic haplotypes as a matrix of shape (ploidy, num_variants),
    entries in {0,1}. Each site has an alternate-allele frequency p ~ Beta(alpha, beta)
    across haplotypes. If ensure_polymorphic=True, resamples sites that end up monomorphic.
    """
    H = np.empty((ploidy, num_variants), dtype=np.uint8)

    # Draw initial site-level allele frequencies
    p_site = rng.beta(maf_alpha, maf_beta, size=num_variants)

    # Sample haplotypes for all sites
    H[:, :] = (rng.random((ploidy, num_variants)) < p_site[None, :]).astype(np.uint8)

    if ensure_polymorphic:
        # Identify monomorphic columns and resample until all polymorphic
        # Polymorphic means: not all 0 and not all 1
        max_iters = 1000
        for _ in range(max_iters):
            col_sums = H.sum(axis=0)
            mono_mask = (col_sums == 0) | (col_sums == ploidy)
            if not np.any(mono_mask):
                break
            k = mono_mask.sum()
            # Resample the problematic columns
            p_resample = rng.beta(maf_alpha, maf_beta, size=k)
            H[:, mono_mask] = (rng.random((ploidy, k)) < p_resample[None, :]).astype(np.uint8)
        # A final check (should be polymorphic now)
        col_sums = H.sum(axis=0)
        mono_mask = (col_sums == 0) | (col_sums == ploidy)
        if np.any(mono_mask):
            raise RuntimeError("Failed to ensure polymorphic sites after resampling.")

    return H


def simulate_reads(
    H: np.ndarray,
    num_reads: int,
    read_length: int,
    error_rate: float,
    missing_rate: float,
    rng: np.random.Generator,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int, int]:
    """
    Simulate reads from haplotypes.

    Returns:
        hap_ids: (R,) int haplotype id for each read (truth)
        starts: (R,) int start variant index (0-based)
        lengths: (R,) int actual length per read (clamped to variant bounds)
        positions: concatenated positions for all reads (1D array)
        alleles_obs: concatenated observed alleles for all reads (-1=missing, 0/1 otherwise)
        error_count: number of mismatches (observed vs truth) excluding missing
        missing_count: number of missing observations
    """
    ploidy, num_variants = H.shape
    if num_reads == 0:
        return (
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int32),
            np.empty(0, dtype=np.int8),
            0,
            0,
        )

    # Choose a haplotype for each read (uniform)
    hap_ids = rng.integers(0, ploidy, size=num_reads, dtype=np.int32)

    # Determine start positions; clamp when read_length > num_variants
    max_start_inclusive = max(0, num_variants - 1)
    if read_length <= num_variants:
        starts = rng.integers(0, num_variants - read_length + 1, size=num_reads, dtype=np.int32)
        lengths = np.full(num_reads, read_length, dtype=np.int32)
    else:
        # All reads start at 0 and are truncated to num_variants
        starts = np.zeros(num_reads, dtype=np.int32)
        lengths = np.full(num_reads, num_variants, dtype=np.int32)

    # Precompute total emitted bases (sum of lengths)
    total_bases = int(lengths.sum())

    # Build concatenated positions array for all reads
    positions = np.empty(total_bases, dtype=np.int32)
    write_ptr = 0
    for s, L in zip(starts, lengths):
        positions[write_ptr : write_ptr + L] = np.arange(s, s + L, dtype=np.int32)
        write_ptr += L

    # Gather truth alleles per read segment into a flat array
    alleles_true = np.empty(total_bases, dtype=np.uint8)
    write_ptr = 0
    for h, s, L in zip(hap_ids, starts, lengths):
        if L == 0:
            continue
        alleles_true[write_ptr : write_ptr + L] = H[h, s : s + L]
        write_ptr += L

    # Apply missing and errors
    rng_missing = rng.random(total_bases)
    missing_mask = rng_missing < missing_rate

    rng_error = rng.random(total_bases)
    error_mask = (~missing_mask) & (rng_error < error_rate)

    alleles_obs = alleles_true.copy().astype(np.int8)
    # Flip alleles where error occurs
    flip_idx = np.where(error_mask)[0]
    alleles_obs[flip_idx] = 1 - alleles_obs[flip_idx]
    # Set missing
    alleles_obs[missing_mask] = -1  # use -1 to denote missing

    # Counters
    observed_mask = ~missing_mask
    error_count = int(np.sum((alleles_obs[observed_mask] != alleles_true[observed_mask])))
    missing_count = int(np.sum(missing_mask))

    return hap_ids, starts, lengths, positions, alleles_obs, error_count, missing_count


def write_haplotypes_tsv(H: np.ndarray, path: str) -> None:
    P, N = H.shape
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"# haplotypes\t{P}\n")
        f.write(f"# variants\t{N}\n")
        f.write("# format: hap_id<TAB>alleles_as_0_1_string\n")
        for h in range(P):
            bitstr = "".join("1" if x else "0" for x in H[h])
            f.write(f"h{h}\t{bitstr}\n")


def write_haplotypes_npz(H: np.ndarray, path: str) -> None:
    np.savez_compressed(path, haplotypes=H)


def write_reads_sparse_tsv(
    hap_ids: np.ndarray,
    starts: np.ndarray,
    lengths: np.ndarray,
    positions: np.ndarray,
    alleles_obs: np.ndarray,
    path: str,
) -> None:
    """
    Write reads to a sparse TSV with columns:
    read_id, hap_truth, start, positions, alleles
    - positions: comma-separated 0-based variant indices for the read segment
    - alleles: string over {0,1,-} aligned to positions
    """
    with open(path, "w", encoding="utf-8") as f:
        f.write("# read_id\thap_truth\tstart\tpositions\talleles\n")
        ptr = 0
        for rid, (h, s, L) in enumerate(zip(hap_ids, starts, lengths)):
            if L == 0:
                f.write(f"r{rid}\t{h}\t{s}\t\t\n")
                continue
            pos_chunk = positions[ptr : ptr + L]
            obs_chunk = alleles_obs[ptr : ptr + L]
            pos_str = ",".join(map(str, pos_chunk.tolist()))
            alleles_str = "".join("-" if a < 0 else ("1" if a == 1 else "0") for a in obs_chunk.tolist())
            f.write(f"r{rid}\t{h}\t{s}\t{pos_str}\t{alleles_str}\n")
            ptr += L


def build_coverage_from_reads(starts: np.ndarray, lengths: np.ndarray, num_variants: int) -> np.ndarray:
    """
    Compute per-variant coverage using a difference-array trick in O(R) time.
    """
    diff = np.zeros(num_variants + 1, dtype=np.int64)
    ends = starts + lengths  # exclusive
    np.add.at(diff, starts, 1)
    np.add.at(diff, ends, -1)
    cov = np.cumsum(diff[:-1])
    return cov


def main():
    parser = argparse.ArgumentParser(description="Generate synthetic phasing reads and haplotypes.")
    parser.add_argument("-p", "--ploidy", type=int, default=2, help="Number of haplotypes (default: 2)")
    parser.add_argument("-n", "--num-variants", type=int, default=1000, help="Number of variant sites (default: 1000)")
    parser.add_argument("-r", "--num-reads", type=int, default=20000, help="Number of reads to simulate (default: 20000)")
    parser.add_argument("-l", "--read-length", type=int, default=25, help="Read length in number of variants (default: 25)")
    parser.add_argument("-e", "--error-rate", type=float, default=0.01, help="Per-site substitution error rate (default: 0.01)")
    parser.add_argument("-m", "--missing-rate", type=float, default=0.05, help="Per-site missing probability (default: 0.05)")
    parser.add_argument("--maf-alpha", type=float, default=0.4, help="Beta(alpha,beta) for site alt-allele freq; alpha (default: 0.4)")
    parser.add_argument("--maf-beta", type=float, default=0.4, help="Beta(alpha,beta) for site alt-allele freq; beta (default: 0.4)")
    parser.add_argument("--allow-monomorphic", action="store_true", help="Allow monomorphic variant sites (default: off)")
    parser.add_argument("-s", "--seed", type=int, default=None, help="Random seed (default: None)")
    parser.add_argument("-o", "--output-prefix", type=str, default="sim", help="Output file prefix (default: sim)")

    args = parser.parse_args()

    params = SimParams(
        ploidy=args.ploidy,
        num_variants=args.num_variants,
        num_reads=args.num_reads,
        read_length=args.read_length,
        error_rate=args.error_rate,
        missing_rate=args.missing_rate,
        maf_alpha=args.maf_alpha,
        maf_beta=args.maf_beta,
        ensure_polymorphic=not args.allow_monomorphic,
        seed=args.seed,
        output_prefix=args.output_prefix,
    )
    _validate_params(params)

    rng = np.random.default_rng(params.seed)

    # Simulate haplotypes
    H = simulate_haplotypes(
        params.ploidy,
        params.num_variants,
        params.maf_alpha,
        params.maf_beta,
        params.ensure_polymorphic,
        rng,
    )

    # Simulate reads
    hap_ids, starts, lengths, positions, alleles_obs, error_count, missing_count = simulate_reads(
        H=H,
        num_reads=params.num_reads,
        read_length=params.read_length,
        error_rate=params.error_rate,
        missing_rate=params.missing_rate,
        rng=rng,
    )

    # Coverage
    coverage = build_coverage_from_reads(starts, lengths, params.num_variants)
    expected_cov = (params.num_reads * params.read_length) / max(1, params.num_variants)
    mean_cov = float(np.mean(coverage)) if coverage.size > 0 else 0.0
    min_cov = int(np.min(coverage)) if coverage.size > 0 else 0
    max_cov = int(np.max(coverage)) if coverage.size > 0 else 0

    # Observed rates
    total_bases = int(lengths.sum())
    non_missing_bases = total_bases - missing_count
    observed_error_rate = (error_count / non_missing_bases) if non_missing_bases > 0 else 0.0
    observed_missing_frac = (missing_count / total_bases) if total_bases > 0 else 0.0

    # Write outputs
    os.makedirs(os.path.dirname(os.path.abspath(params.output_prefix)), exist_ok=True)  # ensure directory exists if prefix includes path

    hap_tsv = f"{params.output_prefix}.haplotypes.tsv"
    hap_npz = f"{params.output_prefix}.haplotypes.npz"
    reads_tsv = f"{params.output_prefix}.reads.sparse.tsv"
    summary_json = f"{params.output_prefix}.summary.json"

    write_haplotypes_tsv(H, hap_tsv)
    write_haplotypes_npz(H, hap_npz)
    write_reads_sparse_tsv(hap_ids, starts, lengths, positions, alleles_obs, reads_tsv)

    summary = SimSummary(
        ploidy=params.ploidy,
        num_variants=params.num_variants,
        num_reads=params.num_reads,
        read_length=params.read_length,
        error_rate=params.error_rate,
        missing_rate=params.missing_rate,
        maf_alpha=params.maf_alpha,
        maf_beta=params.maf_beta,
        ensure_polymorphic=params.ensure_polymorphic,
        seed=params.seed,
        expected_coverage=float(expected_cov),
        mean_coverage=mean_cov,
        min_coverage=min_cov,
        max_coverage=max_cov,
        observed_error_rate_over_non_missing=float(observed_error_rate),
        observed_missing_fraction_over_all_bases=float(observed_missing_frac),
    )

    with open(summary_json, "w", encoding="utf-8") as f:
        json.dump(asdict(summary), f, indent=2)

    # Friendly stdout
    print("Wrote:")
    print(f"  haplotypes (tsv): {hap_tsv}")
    print(f"  haplotypes (npz): {hap_npz} (array key: 'haplotypes', shape={tuple(H.shape)})")
    print(f"  reads (sparse tsv): {reads_tsv}")
    print(f"  summary: {summary_json}")
    print("")
    print("Key summary:")
    print(json.dumps(asdict(summary), indent=2))


if __name__ == "__main__":
    main()