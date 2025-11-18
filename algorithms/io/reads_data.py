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
        fragments = []
        for i in range(self.reads.shape[0]):
            # Determine how many entries in this read (positions matrix may have padding -1 values)
            pos_row = self.positions[i]
            val_row = self.reads[i]
            valid_mask = (pos_row != -1)
            count = int(np.sum(valid_mask))
            indices = pos_row[:count].tolist()
            values = val_row[:count].tolist()
            fragments.append({"id": i, "indices": indices, "values": values})
        # Use the writer utility to output TSV
        from algorithms.io.writer import write_reads_sparse_tsv
        write_reads_sparse_tsv(filepath, fragments)