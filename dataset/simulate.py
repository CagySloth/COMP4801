import argparse
import numpy as np
from pathlib import Path

from algorithms.io.writer import (
    write_haplotypes_tsv,
    write_reads_sparse_tsv,
    write_haplotypes_npz,
)
from algorithms.io.reads_data import ReadsData


def generate_diploid_data(num_variants, num_reads, read_length, error_rate, missing_rate, allow_monomorphic):
    hap1 = np.random.randint(0, 2, size=num_variants)
    hap2 = hap1.copy()
    switch_indices = np.random.rand(num_variants) < 0.1
    hap2[switch_indices] ^= 1

    if not allow_monomorphic and np.all(hap1 == hap2):
        hap2[0] ^= 1  # ensure difference

    haplotypes = np.stack([hap1, hap2])

    reads = []
    for _ in range(num_reads):
        start = np.random.randint(0, num_variants - read_length + 1)
        hap_idx = np.random.randint(2)
        segment = haplotypes[hap_idx, start : start + read_length].copy()

        noise = np.random.rand(read_length) < error_rate
        segment[noise] ^= 1
        missing = np.random.rand(read_length) < missing_rate
        segment[missing] = -1

        indices = list(range(start, start + read_length))
        values = segment.tolist()
        reads.append({"id": len(reads), "indices": indices, "values": values})

    return haplotypes, reads


def generate_polyploid_data(num_variants, ploidy, num_reads, read_length, error_rate, missing_rate, alpha, beta, allow_monomorphic):
    haplotypes = np.random.binomial(1, np.random.beta(alpha, beta, size=num_variants), size=(ploidy, num_variants))
    
    if not allow_monomorphic:
        while np.all(haplotypes == haplotypes[0]):
            haplotypes[np.random.randint(ploidy), np.random.randint(num_variants)] ^= 1

    reads = []
    for i in range(num_reads):
        start = np.random.randint(0, num_variants - read_length + 1)
        hap = haplotypes[np.random.randint(ploidy)]
        segment = hap[start : start + read_length].copy()

        # Apply errors and missing
        noise = np.random.rand(read_length) < error_rate
        segment[noise] ^= 1
        missing = np.random.rand(read_length) < missing_rate
        segment[missing] = -1

        indices = list(range(start, start + read_length))
        values = segment.tolist()
        reads.append({"id": i, "indices": indices, "values": values})

    return haplotypes, reads


def main():
    parser = argparse.ArgumentParser(description="Simulate haplotype and read data.")
    parser.add_argument("-p", "--ploidy", type=int, required=True)
    parser.add_argument("-n", "--num-variants", type=int, required=True)
    parser.add_argument("-r", "--num-reads", type=int, required=True)
    parser.add_argument("-l", "--read-length", type=int, required=True)
    parser.add_argument("-e", "--error-rate", type=float, default=0.01)
    parser.add_argument("-m", "--missing-rate", type=float, default=0.0)
    parser.add_argument("--maf-alpha", type=float, default=1.0)
    parser.add_argument("--maf-beta", type=float, default=1.0)
    parser.add_argument("--allow-monomorphic", action="store_true")
    parser.add_argument("-o", "--output-prefix", required=True)

    args = parser.parse_args()

    if args.ploidy == 2:
        haps, reads = generate_diploid_data(
            args.num_variants,
            args.num_reads,
            args.read_length,
            args.error_rate,
            args.missing_rate,
            args.allow_monomorphic,
        )
    else:
        haps, reads = generate_polyploid_data(
            args.num_variants,
            args.ploidy,
            args.num_reads,
            args.read_length,
            args.error_rate,
            args.missing_rate,
            args.maf_alpha,
            args.maf_beta,
            args.allow_monomorphic,
        )

    prefix = Path(args.output_prefix)
    write_haplotypes_tsv(f"{prefix}.haplotypes.tsv", haps)
    write_reads_sparse_tsv(f"{prefix}.reads.sparse.tsv", reads)
    reads_np = ReadsData.from_fragments(reads)
    write_haplotypes_npz(f"{prefix}.reads.npz", reads_np.reads)

    print(f"✅ Simulated data written to: {prefix}.[haplotypes.tsv, reads.sparse.tsv, reads.npz]")


if __name__ == "__main__":
    main()