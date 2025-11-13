from dataclasses import dataclass
import numpy as np
from typing import Optional


@dataclass
class ReadsData:
    reads: np.ndarray
    positions: Optional[np.ndarray] = None
    weights: Optional[np.ndarray] = None

    def __post_init__(self):
        if not isinstance(self.reads, np.ndarray):
            raise TypeError("reads must be a numpy array")
        if self.reads.ndim != 2:
            raise ValueError("reads must be a 2D numpy array")

        n_reads, n_variants = self.reads.shape
        if self.positions is not None and len(self.positions) != n_variants:
            raise ValueError("positions length must match number of columns in reads")
        if self.weights is not None and len(self.weights) != n_reads:
            raise ValueError("weights length must match number of rows in reads")

    @classmethod
    def from_fragments(cls, fragments: list[dict]) -> "ReadsData":
        # Infer max variant index
        max_idx = max(i for frag in fragments for i in frag["indices"]) + 1
        reads = np.full((len(fragments), max_idx), -1, dtype=int)

        for i, frag in enumerate(fragments):
            for j, val in zip(frag["indices"], frag["values"]):
                reads[i, j] = val

        return cls(reads)