import numpy as np
from dataset.simulate import generate_diploid_data
from algorithms.io.reads_data import ReadsData


def test_readsdata_dense_matrix_semantics():
    haps, reads = generate_diploid_data(
        num_variants=5,
        num_reads=3,
        read_length=3,
        error_rate=0.0,
        missing_rate=0.0,
        allow_monomorphic=True,
    )

    data = ReadsData.from_fragments(reads)

    # shape checks
    assert data.alleles.shape == (3, 5)
    assert data.positions.shape == (3, 5)

    # value domain checks
    assert np.all(np.isin(data.alleles, [-1, 0, 1]))

    # positions should be 0-based or -1
    valid_positions = (data.positions >= 0) | (data.positions == -1)
    assert np.all(valid_positions)

from algorithms.diploid.whatshap_adapter import build_readset_from_readsdata

def test_readset_matches_dense_matrix():
    _, reads = generate_diploid_data(
        num_variants=6,
        num_reads=4,
        read_length=3,
        error_rate=0.0,
        missing_rate=0.0,
        allow_monomorphic=True,
    )

    data = ReadsData.from_fragments(reads)
    readset = build_readset_from_readsdata(data)

    expected = sum(not np.all(data.alleles[i] < 0) for i in range(data.R))
    assert len(readset) == expected

    for r in readset:
        last_pos = -1
        for v in r:
            # allele correctness
            assert v.allele in (0, 1)

            # sorted by position
            assert v.position > last_pos
            last_pos = v.position
