from __future__ import annotations

from itertools import combinations

from dataset.longread.reference import generate_reference


def _overlaps(a0: int, a1: int, b0: int, b1: int) -> bool:
    return max(a0, b0) < min(a1, b1)


def test_generate_reference_duplications_are_valid():
    seq, regions, duplications = generate_reference(
        length=80000,
        seed=0,
        contig="chr1",
        preset="plain",
        homopolymer_count=0,
        homopolymer_min_len=6,
        homopolymer_max_len=20,
        str_count=0,
        str_min_motif=2,
        str_max_motif=6,
        str_min_repeats=3,
        str_max_repeats=10,
        gc_window_count=0,
        gc_window_len=200,
        gc_high=0.75,
        gc_low=0.25,
        dup_segments=5,
        dup_len=3000,
        dup_min_gap=500,
    )

    assert isinstance(seq, str)
    assert len(seq) == 80000
    assert len(duplications) == 5

    intervals = []
    for d in duplications:
        ss, se = int(d["src_start0"]), int(d["src_end0"])
        ds, de = int(d["dst_start0"]), int(d["dst_end0"])

        assert se - ss == 3000
        assert de - ds == 3000

        # dst is an exact copy of src
        assert seq[ss:se] == seq[ds:de]

        intervals.append(("src", ss, se))
        intervals.append(("dst", ds, de))

    for (t1, a0, a1), (t2, b0, b1) in combinations(intervals, 2):
        assert not _overlaps(a0, a1, b0, b1), f"Overlap: {t1}[{a0},{a1}) vs {t2}[{b0},{b1})"