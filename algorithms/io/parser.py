import numpy as np
from pathlib import Path
from algorithms.io.reads_data import ReadsData


def parse_sparse_tsv(path):
    fragments = []
    with open(path) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 2:
                continue
            read_id = int(fields[0])
            indices, values = zip(*(field.split(":") for field in fields[1:]))
            fragment = {
                "id": read_id,
                "indices": list(map(int, indices)),
                "values": list(map(int, values)),
            }
            fragments.append(fragment)
    return fragments


def parse_dense_tsv(path):
    reads = []
    with open(path) as f:
        for line in f:
            row = [int(x) for x in line.strip()]
            reads.append(row)
    return np.array(reads, dtype=int)


def load_reads(path: str | Path) -> ReadsData:
    path = Path(path)
    if path.suffix == ".npz":
        data = np.load(path)
        return ReadsData(
            reads=data["reads"],
            positions=data.get("positions", None),
            weights=data.get("weights", None),
        )
    elif path.suffix == ".tsv":
        if is_sparse_tsv(path):
            fragments = parse_sparse_tsv(path)
            return ReadsData.from_fragments(fragments)
        else:
            reads = parse_dense_tsv(path)
            return ReadsData(reads=reads)
    else:
        raise ValueError(f"Unsupported file format: {path}")


def is_sparse_tsv(path: Path) -> bool:
    # Simple heuristic: check for colons (sparse format has indices:values)
    with open(path) as f:
        for line in f:
            if ":" in line:
                return True
            if line.strip() != "":
                return False
    return False


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