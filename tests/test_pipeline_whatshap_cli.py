# tests/test_pipeline_whats_cli.py

import os
import json
import subprocess
import sys
import pytest

try:
    from vendor.whcore.py import core  # noqa: F401
    HAS_WHCORE = True
except ImportError:
    HAS_WHCORE = False


pytestmark = pytest.mark.skipif(
    not HAS_WHCORE,
    reason="WhatsHap core (vendor.whcore) not available / not built",
)


def test_full_pipeline_with_diploid_whats(tmp_path):
    """
    End-to-end:
      - simulate a small diploid dataset
      - run `diploid-whats` via the unified CLI
      - check that WhatsHap-based outputs exist and are parseable.
    """
    outprefix = tmp_path / "whatsmini"
    ploidy = 2
    num_variants = 20
    num_reads = 40

    # 1) Simulate data
    sim_cmd = [
        sys.executable,
        "dataset/simulate.py",
        "--ploidy", str(ploidy),
        "--num-variants", str(num_variants),
        "--num-reads", str(num_reads),
        "--read-length", "10",
        "--error-rate", "0.0",
        "--missing-rate", "0.0",
        "--output-prefix", str(outprefix),
    ]
    subprocess.run(sim_cmd, check=True)

    npz_file = f"{outprefix}.reads.npz"
    assert os.path.exists(npz_file)

    # 2) Phase with diploid-whats
    phase_cmd = [
        sys.executable,
        "-m", "algorithms.cli.phase",
        "diploid-whats",
        "-i", npz_file,
        "--output-prefix", f"{outprefix}.whats",
    ]
    result = subprocess.run(phase_cmd)
    assert result.returncode == 0, "diploid-whats CLI returned non-zero exit status"

    hap_file = f"{outprefix}.whats.haplotypes.tsv"
    assign_file = f"{outprefix}.whats.assignments.tsv"
    summary_file = f"{outprefix}.whats.summary.json"

    assert os.path.exists(hap_file)
    assert os.path.exists(assign_file)
    assert os.path.exists(summary_file)

    with open(summary_file) as f:
        summary = json.load(f)
    assert "algorithm" in summary
    assert summary["algorithm"].startswith("diploid_whats")
