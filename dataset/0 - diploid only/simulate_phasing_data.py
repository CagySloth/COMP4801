#!/usr/bin/env python3
"""
simple_sim_phasing.py

Diploid-only synthetic data generator for phasing benchmarks.
Simplifications:
- Every read covers all variant sites (read length = num_variants).
- Direct control over per-site homozygosity via --homozygous-rate.

Outputs:
- <prefix>.haplotypes.tsv     Ground-truth haplotypes (2 rows, N columns as 0/1 string)
- <prefix>.reads.tsv          One read per row: read_id, hap_truth, alleles (0/1/-)
- <prefix>.summary.json       Basic summary stats

Allele conventions:
- Biallelic 0/1, missing denoted '-'.
- Missing applied first; substitution errors applied only to non-missing bases.
"""

import argparse
import json
import os
from typing import Tuple

import numpy as np


def simulate_haplotypes_diploid(
    num_variants: int,
    homozygous_rate: float,
    hom_alt_prob: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Create a 2 x num_variants haplotype matrix with entries in {0,1}.

    - With probability homozygous_rate, a site is homozygous: [0,0] or [1,1].
      [1,1] is chosen with probability hom_alt_prob; [0,0] otherwise.
    - Otherwise the site is heterozygous: [0,1] or [1,0] (50/50 orientation).
    """
    is_hom = rng.random(num_variants) < homozygous_rate

    # For homozygous sites, choose whether they are 1/1 or 0/0
    alt_when_hom = (rng.random(num_variants) < hom_alt_prob).astype(np.uint8)

    H0 = np.empty(num_variants, dtype=np.uint8)
    H1 = np.empty(num_variants, dtype=np.uint8)

    # Homozygous assignments
    H0[is_hom] = alt_when_hom[is_hom]
    H1[is_hom] = alt_when_hom[is_hom]

    # Heterozygous assignments: pick orientation randomly
    heter_mask = ~is_hom
    if np.any(heter_mask):
        h0_heter = rng.integers(0, 2, size=int(heter_mask.sum()), dtype=np.uint8)  # 0 or 1
        H0[heter_mask] = h0_heter
        H1[heter_mask] = 1 - h0_heter

    H = np.stack([H0, H1], axis=0)  # shape (2, N)
    return H


def simulate_reads_cover_all(
    H: np.ndarray,
    num_reads: int,
    error_rate: float,
    missing_rate: float,
    rng: np.random.Generator,
) -> Tuple[np.ndarray, np.ndarray, int, int]:
    """
    Simulate num_reads reads, each covering all variant sites.
    Returns:
        hap_ids: (R,) the true haplotype index (0 or 1) per read
        obs:     (R, N) int8 observed alleles with values in {-1,0,1} (-1=missing)
        error_count: number of mismatches (excluding missing)
        missing_count: number of missing observations
    """
    assert H.shape[0] == 2
    N = H.shape[1]
    R = num_reads

    if R == 0 or N == 0:
        return np.empty(0, dtype=np.int32), np.empty((0, N), dtype=np.int8), 0, 0

    # Choose haplotypes for each read
    hap_ids = rng.integers(0, 2, size=R, dtype=np.int32)

    # True alleles for each read (R x N)
    truth = H[hap_ids, :].astype(np.uint8)

    # Missing and errors
    miss_mask = rng.random((R, N)) < missing_rate
    err_mask = (~miss_mask) & (rng.random((R, N)) < error_rate)

    obs = truth.astype(np.int8)
    # Flip alleles where error occurs
    obs[err_mask] = 1 - obs[err_mask]
    # Apply missing
    obs[miss_mask] = -1

    error_count = int(err_mask.sum())
    missing_count = int(miss_mask.sum())
    return hap_ids, obs, error_count, missing_count


def write_haplotypes_tsv(H: np.ndarray, path: str) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"# haplotypes\t{H.shape[0]}\n")
        f.write(f"# variants\t{H.shape[1]}\n")
        f.write("# format: hap_id<TAB>alleles_as_0_1_string\n")
        for h in range(H.shape[0]):
            bitstr = "".join("1" if x else "0" for x in H[h])
            f.write(f"h{h}\t{bitstr}\n")


def write_reads_tsv(hap_ids: np.ndarray, obs: np.ndarray, path: str) -> None:
    """
    Writes: read_id, hap_truth, alleles_string
    alleles_string is length N over {0,1,-}.
    """
    def row_to_string(row: np.ndarray) -> str:
        return "".join("-" if a < 0 else ("1" if a == 1 else "0") for a in row.tolist())

    with open(path, "w", encoding="utf-8") as f:
        f.write("# read_id\thap_truth\talleles\n")
        for i in range(obs.shape[0]):
            f.write(f"r{i}\t{hap_ids[i]}\t{row_to_string(obs[i])}\n")


def main():
    ap = argparse.ArgumentParser(description="Simple diploid phasing simulator (reads cover all variants).")
    ap.add_argument("-n", "--num-variants", type=int, default=1000, help="Number of variant sites (default: 1000)")
    ap.add_argument("-r", "--num-reads", type=int, default=2000, help="Number of reads (default: 2000)")
    ap.add_argument("-e", "--error-rate", type=float, default=0.01, help="Per-site substitution error prob (default: 0.01)")
    ap.add_argument("-m", "--missing-rate", type=float, default=0.05, help="Per-site missing prob (default: 0.05)")
    ap.add_argument("--homozygous-rate", type=float, default=0.5, help="Per-site probability of homozygosity (default: 0.5)")
    ap.add_argument("--hom-alt-prob", type=float, default=0.5, help="If homozygous, prob of [1,1] (else [0,0]) (default: 0.5)")
    ap.add_argument("-o", "--output-prefix", type=str, default="simple", help="Output file prefix (default: simple)")
    ap.add_argument("-s", "--seed", type=int, default=None, help="Random seed (default: None)")
    args = ap.parse_args()

    # Basic validation
    if args.num_variants < 1:
        raise ValueError("num_variants must be >= 1")
    if args.num_reads < 0:
        raise ValueError("num_reads must be >= 0")
    for name, val in [("error_rate", args.error_rate), ("missing_rate", args.missing_rate),
                      ("homozygous_rate", args.homozygous_rate), ("hom_alt_prob", args.hom_alt_prob)]:
        if not (0.0 <= val <= 1.0):
            raise ValueError(f"{name} must be in [0,1]")

    rng = np.random.default_rng(args.seed)

    # Simulate haplotypes and reads
    H = simulate_haplotypes_diploid(
        num_variants=args.num_variants,
        homozygous_rate=args.homozygous_rate,
        hom_alt_prob=args.hom_alt_prob,
        rng=rng,
    )

    hap_ids, obs, error_count, missing_count = simulate_reads_cover_all(
        H=H,
        num_reads=args.num_reads,
        error_rate=args.error_rate,
        missing_rate=args.missing_rate,
        rng=rng,
    )

    # Summaries
    R, N = obs.shape
    total_bases = R * N
    non_missing = total_bases - missing_count
    observed_error_rate = (error_count / non_missing) if non_missing > 0 else 0.0
    observed_missing_frac = (missing_count / total_bases) if total_bases > 0 else 0.0
    heter_sites = int(np.sum(H[0] != H[1]))
    hom_sites = N - heter_sites

    # Ensure output directory exists if a path is in the prefix
    outdir = os.path.dirname(os.path.abspath(args.output_prefix))
    os.makedirs(outdir, exist_ok=True)

    hap_path = f"{args.output_prefix}.haplotypes.tsv"
    reads_path = f"{args.output_prefix}.reads.tsv"
    summary_path = f"{args.output_prefix}.summary.json"

    write_haplotypes_tsv(H, hap_path)
    write_reads_tsv(hap_ids, obs, reads_path)

    summary = {
        "num_variants": N,
        "num_reads": R,
        "error_rate": args.error_rate,
        "missing_rate": args.missing_rate,
        "homozygous_rate_param": args.homozygous_rate,
        "hom_alt_prob_param": args.hom_alt_prob,
        "seed": args.seed,
        "observed_missing_fraction_over_all_bases": observed_missing_frac,
        "observed_error_rate_over_non_missing": observed_error_rate,
        "haplotype_homozygous_sites": hom_sites,
        "haplotype_heterozygous_sites": heter_sites,
    }
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("Wrote:")
    print(f"  haplotypes: {hap_path}")
    print(f"  reads:      {reads_path}")
    print(f"  summary:    {summary_path}")
    print("\nKey summary:")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()