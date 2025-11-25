# tests/test_whatshap_adapter.py

import numpy as np
import pytest

from algorithms.io import ReadsData

try:
    from algorithms.diploid.whatshap_adapter import build_readset_from_readsdata
    from vendor.whcore.py import core
    HAS_WHCORE = True
except ImportError:
    HAS_WHCORE = False


pytestmark = pytest.mark.skipif(
    not HAS_WHCORE,
    reason="WhatsHap core (vendor.whcore) not available / not built",
)


def test_build_readset_from_readsdata_basic():
    """
    Simple sanity check that a dense 0/1 matrix turns into a ReadSet with
    matching positions and alleles.
    """
    # 3 reads, 4 variants; no missing values
    reads = np.array(
        [
            [0, 1, 0, 1],
            [1, 0, 1, 0],
            [0, 0, 1, 1],
        ],
        dtype=int,
    )
    R, N = reads.shape
    # Dense positions: every read covers all variants 0..N-1
    positions = np.tile(np.arange(N), (R, 1))

    data = ReadsData(reads=reads, positions=positions, num_variants=N)

    readset = build_readset_from_readsdata(data)

    # One Read per row
    assert len(readset) == R

    # For each read, every variant should match the source matrix
    for r_idx, read in enumerate(readset):
        for variant in read:
            pos = variant.position
            allele = variant.allele
            assert 0 <= pos < N
            assert allele == reads[r_idx, pos]
