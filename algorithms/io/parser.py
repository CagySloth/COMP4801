import numpy as np
import itertools

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
            reads = parse_dense_tsv(path)           # NumPy array of shape (R, N)
            R, N = reads.shape
            positions = np.tile(np.arange(N), (R, 1))
            return ReadsData(reads=reads, positions=positions, num_variants=N)

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