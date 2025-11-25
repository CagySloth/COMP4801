# tests/test_diploid_whatshap.py

import numpy as np
import os
import json
import pytest

try:
    from algorithms.diploid import whatshap_driver
    from vendor.whcore.py import core  # ensure core is importable
    HAS_WHCORE = True
except ImportError:
    HAS_WHCORE = False


pytestmark = pytest.mark.skipif(
    not HAS_WHCORE,
    reason="WhatsHap core or driver not available / not built",
)


def _write_npz(tmp_path, reads_matrix, name="input.npz"):
    path = tmp_path / name
    np.savez_compressed(path, reads=reads_matrix)
    return path


def _load_hap_strings(path):
    # writer.write_haplotypes_tsv writes: "i\t0101..."
    hap_strings = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            assert len(parts) == 2
            hap_strings.append(parts[1])
    return hap_strings


def test_diploid_whatshap_two_reads_two_variants(tmp_path):
    """
    Tiny toy case:
      - true haps could be [0,1] and [1,0]
      - two reads, each equal to one haplotype
    We only assert that the output has 2 haplotypes whose columns are {0,1}.
    """
    reads_matrix = np.array(
        [
            [0, 1],  # read 0
            [1, 0],  # read 1
        ],
        dtype=int,
    )

    npz_file = _write_npz(tmp_path, reads_matrix, "wh_input.npz")
    out_prefix = tmp_path / "wh_out"

    # Fake args object like the CLI will pass
    args = type("Args", (), {})()
    args.input = str(npz_file)
    args.output_prefix = str(out_prefix)
    args.max_coverage = 10
    args.error_rate = 0.01

    whatshap_driver.main(args)

    hap_file = str(out_prefix) + ".haplotypes.tsv"
    assign_file = str(out_prefix) + ".assignments.tsv"
    summary_file = str(out_prefix) + ".summary.json"

    assert os.path.exists(hap_file), "Haplotype TSV not written"
    assert os.path.exists(assign_file), "Assignments TSV not written"
    assert os.path.exists(summary_file), "Summary JSON not written"

    haps = _load_hap_strings(hap_file)
    assert len(haps) == 2, "Diploid phasing should output 2 haplotypes"

    # Each hap string has length 2 and characters in {0,1} (no -1 for this tiny clean case)
    assert all(len(h) == 2 for h in haps)
    assert all(set(h).issubset({"0", "1"}) for h in haps)

    col0 = {haps[0][0], haps[1][0]}
    col1 = {haps[0][1], haps[1][1]}
    # At each position we should see both alleles
    assert col0 == {"0", "1"}
    assert col1 == {"0", "1"}

    # Summary JSON at least has algorithm / R / N
    with open(summary_file) as f:
        summary = json.load(f)
    assert summary.get("algorithm", "").startswith("diploid_whats")
    assert summary.get("R") == 2
    assert summary.get("N") == 2
