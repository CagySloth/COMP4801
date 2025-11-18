# algorithms/io/reads_data.py
from dataclasses import dataclass
import numpy as np

@dataclass
class ReadsData:
    reads: np.ndarray                # shape (num_reads, read_length)
    positions: np.ndarray            # shape (num_reads, read_length)
    num_variants: int
    hap_truth: Optional[np.ndarray] = None
    
    @property
    def alleles(self):
        return self.reads
    
    @property
    def R(self):
        return self.reads.shape[0]

    @property
    def N(self):
        return self.reads.shape[1]

    @staticmethod
    def from_fragments(fragments: list[dict]) -> "ReadsData":
        """
        Convert list of sparse fragment dicts into dense reads and positions matrices.
        Each fragment must have keys: 'id', 'indices', 'values'
        """
        num_reads = len(fragments)
        max_index = max(idx for f in fragments for idx in f["indices"])
        num_variants = max_index + 1

        reads = np.full((num_reads, num_variants), -1, dtype=int)
        positions = np.full((num_reads, num_variants), -1, dtype=int)

        for i, frag in enumerate(fragments):
            for j, (idx, val) in enumerate(zip(frag["indices"], frag["values"])):
                reads[i, idx] = val
                positions[i, j] = idx

        return ReadsData(reads=reads, positions=positions, num_variants=num_variants)

    @staticmethod
    def from_npz(path: str) -> "ReadsData":
        data = np.load(path)
        reads = data["reads"]
        positions = data.get("positions")  # may be None if not stored
        if positions is None:
            # If no positions saved (e.g. dense data saved without positions), assume dense full coverage
            num_variants = reads.shape[1] if reads.size else 0
            # Each read covers all variant indices 0..N-1
            positions = np.tile(np.arange(num_variants), (reads.shape[0], 1))
        else:
            num_variants = int(np.max(positions[positions >= 0]) + 1) if positions.size else 0
        return ReadsData(reads=reads, positions=positions, num_variants=num_variants)
    
    def to_npz(self, filepath: str):
        # Save reads and positions arrays (and any other relevant data)
        np.savez_compressed(filepath, 
                            reads=self.reads, 
                            positions=self.positions)
    
    def to_sparse_tsv(self, filepath: str):
        with open(filepath, "w") as f:
            for i in range(self.R):
                row = []
                for j in range(self.reads.shape[1]):
                    val = self.reads[i, j]
                    if val != -1:
                        row.append(f"{j}:{val}")
                f.write(f"{i}\t" + "\t".join(row) + "\n")
