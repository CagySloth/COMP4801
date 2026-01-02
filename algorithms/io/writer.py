import numpy as np
import json
from pathlib import Path
from typing import Any


def _ensure_parent_dir(path: str) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)


def write_haplotypes_tsv(path, haplotypes):
    _ensure_parent_dir(path)
    with open(path, "w") as f:
        for i, row in enumerate(haplotypes):
            f.write(f"{i}\t" + "".join(str(x) for x in row) + "\n")


def write_assignments_tsv(path, assignments):
    _ensure_parent_dir(path)
    with open(path, "w") as f:
        for a in assignments:
            f.write(str(a) + "\n")


def write_reads_sparse_tsv(path, fragments):
    _ensure_parent_dir(path)
    with open(path, "w") as f:
        for fragment in fragments:
            out = [str(fragment["id"])]
            for i, val in zip(fragment["indices"], fragment["values"]):
                out.append(f"{i}:{val}")
            f.write("\t".join(out) + "\n")


def write_haplotypes_npz(path, haplotypes, assignments=None):
    _ensure_parent_dir(path)
    np.savez_compressed(path, haplotypes=haplotypes, assignments=assignments)


def write_summary_json(summary_dict: dict[str, Any], path: str):
    """Write a full summary dictionary to a JSON file."""
    _ensure_parent_dir(path)
    with open(path, "w") as f:
        json.dump(summary_dict, f, indent=2)



# import numpy as np
# import json
# from pathlib import Path
# from typing import Optional, Any

# def write_haplotypes_tsv(path, haplotypes):
#     with open(path, "w") as f:
#         for i, row in enumerate(haplotypes):
#             f.write(f"{i}\t" + "".join(str(x) for x in row) + "\n")

# def write_assignments_tsv(path, assignments):
#     with open(path, "w") as f:
#         for a in assignments:
#             f.write(str(a) + "\n")

# def write_reads_sparse_tsv(path, fragments):
#     with open(path, "w") as f:
#         for fragment in fragments:
#             out = [str(fragment['id'])]
#             for i, val in zip(fragment['indices'], fragment['values']):
#                 out.append(f"{i}:{val}")
#             f.write("\t".join(out) + "\n")

# def write_haplotypes_npz(path, haplotypes, assignments=None):
#     np.savez_compressed(path, haplotypes=haplotypes, assignments=assignments)

# def write_summary_json(summary_dict: dict[str, Any], path: str):
#     """Write a full summary dictionary to a JSON file."""
#     with open(path, "w") as f:
#         json.dump(summary_dict, f, indent=2)