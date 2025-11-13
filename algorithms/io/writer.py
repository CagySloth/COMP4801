import numpy as np
import json
from pathlib import Path
from typing import Optional

def write_haplotypes_tsv(path, haplotypes):
    with open(path, "w") as f:
        for row in haplotypes:
            f.write("".join(str(x) for x in row) + "\n")

def write_assignments_tsv(path, assignments):
    with open(path, "w") as f:
        for a in assignments:
            f.write(str(a) + "\n")

def write_reads_sparse_tsv(path, fragments):
    with open(path, "w") as f:
        for fragment in fragments:
            out = [str(fragment['id'])]
            for i, val in zip(fragment['indices'], fragment['values']):
                out.append(f"{i}:{val}")
            f.write("\t".join(out) + "\n")

def write_haplotypes_npz(path, haplotypes, assignments=None):
    np.savez_compressed(path, haplotypes=haplotypes, assignments=assignments)

def write_summary_json(path, mec: int, time_sec: float, acc: Optional[float] = None):
    summary = {"mec": mec, "time": time_sec}
    if acc is not None:
        summary["accuracy"] = acc
    with open(path, "w") as f:
        json.dump(summary, f, indent=2)
