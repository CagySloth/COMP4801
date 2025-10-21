import argparse
import json
from pathlib import Path

import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description="Benchmark phasing accuracy against ground truth.")
    parser.add_argument("--truth", type=Path, required=True, help="Path to ground truth haplotypes TSV")
    parser.add_argument("--pred", type=Path, required=True, help="Path to predicted haplotypes TSV")
    parser.add_argument("--output", type=Path, default=None, help="Optional path to save accuracy report as JSON")
    return parser.parse_args()


def read_haplotypes_tsv(path):
    haplotypes = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) != 2:
                raise ValueError(f"Invalid line: {line.strip()}")
            haplotypes.append([int(x) for x in parts[1]])
    return np.array(haplotypes, dtype=np.uint8)


def best_label_mapping(truth, pred):
    from itertools import permutations

    P = truth.shape[0]
    best_acc = 0.0
    best_perm = None

    for perm in permutations(range(P)):
        aligned = pred[list(perm)]
        acc = np.mean(aligned == truth)
        if acc > best_acc:
            best_acc = acc
            best_perm = perm

    return best_acc, best_perm


def main():
    args = parse_args()

    truth = read_haplotypes_tsv(args.truth)
    pred = read_haplotypes_tsv(args.pred)

    if truth.shape != pred.shape:
        raise ValueError(f"Mismatch in shape: truth {truth.shape}, pred {pred.shape}")

    acc, perm = best_label_mapping(truth, pred)

    print(f"Best accuracy: {acc:.4f}")
    print(f"Best label permutation: {perm}")

    if args.output:
        result = {
            "accuracy": acc,
            "label_permutation": perm,
            "truth_shape": truth.shape,
        }
        with open(args.output, "w") as f:
            json.dump(result, f, indent=2)


if __name__ == "__main__":
    main()