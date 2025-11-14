import numpy as np
import itertools


def global_majority(reads: np.ndarray) -> np.ndarray:
    # Return consensus across reads, ignoring -1 (missing)
    num_variants = reads.shape[1]
    consensus = []
    for j in range(num_variants):
        col = reads[:, j]
        col = col[col != -1]
        if len(col) == 0:
            consensus.append(0)
        else:
            consensus.append(int(np.sum(col) >= len(col) / 2))
    return np.array(consensus)


def compute_mec(reads: np.ndarray, haplotypes: np.ndarray, assignments: np.ndarray) -> int:
    mec = 0
    for i, read in enumerate(reads):
        hap = haplotypes[assignments[i]]
        mismatches = (read != hap) & (read != -1)
        mec += np.sum(mismatches)
    return int(mec)


def hap_truth_accuracy(pred_haps: np.ndarray, true_haps: np.ndarray) -> float:
    if pred_haps.shape != true_haps.shape:
        raise ValueError("Shape mismatch between predicted and true haplotypes")

    perms = [np.array(p) for p in itertools.permutations(range(true_haps.shape[0]))]
    best_score = 0
    for p in perms:
        matches = (pred_haps == true_haps[p])
        score = np.sum(matches)
        best_score = max(best_score, score)
    total = np.prod(pred_haps.shape)
    return best_score / total