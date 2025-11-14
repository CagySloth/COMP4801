from dataclasses import dataclass
import numpy as np

@dataclass
class ReadsData:
    reads: np.ndarray                # shape (num_reads, read_length)
    positions: np.ndarray            # shape (num_reads, read_length)
    num_variants: int

    @staticmethod
    def from_fragments(fragments: list[dict]) -> "ReadsData":
        """
        Convert list of sparse fragment dicts into dense reads and positions matrices.
        Each fragment must have keys: 'id', 'indices', 'values'
        """
        num_reads = len(fragments)
        max_length = max(len(frag["indices"]) for frag in fragments)
        all_positions = np.full((num_reads, max_length), -1, dtype=int)
        all_reads = np.full((num_reads, max_length), -1, dtype=int)

        all_variant_indices = []

        for i, frag in enumerate(fragments):
            indices = frag["indices"]
            values = frag["values"]
            all_positions[i, :len(indices)] = indices
            all_reads[i, :len(values)] = values
            all_variant_indices.extend(indices)

        num_variants = max(all_variant_indices) + 1 if all_variant_indices else 0

        return ReadsData(reads=all_reads, positions=all_positions, num_variants=num_variants)
