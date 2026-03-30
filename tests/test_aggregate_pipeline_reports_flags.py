from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

from benchmark import aggregate_pipeline_reports


def test_aggregate_includes_snp_only_flags(tmp_path: Path):
    root = tmp_path / "runs"
    root.mkdir()

    pipe = {
        "prefix": "output/demo",
        "seed": 0,
        "params": {
            "ref_preset": "plain",
            "ref_length": 20000,
            "num_snps": 200,
            "het_rate": 0.8,
            "platform": "ont",
            "ont_profile": "q20",
            "num_reads": 200,
            "min_len": 2000,
            "max_len": 6000,
            "dup_segments": 0,
            "phase_snps_only": True,
            "eval_snps_only": True,
        },
        "callset": {"precision": 1.0, "recall": 0.9, "truth_snps": 200, "shared_snps": 180},
        "phasing_runs": {
            "called": {
                "eval": {"phase_accuracy_blockflip": 1.0, "switch_error_rate": 0.0, "num_phase_sets": 1},
                "derived": {"effective_phased_recall": 0.8, "shared_het_recall": 0.9},
            },
            "oracle": {
                "eval": {"phase_accuracy_blockflip": 1.0, "switch_error_rate": 0.0, "num_phase_sets": 1},
                "derived": {"effective_phased_recall": 0.9, "shared_het_recall": 0.95},
            },
        },
    }

    p = root / "demo.pipeline.json"
    p.write_text(json.dumps(pipe, indent=2))

    out_csv = tmp_path / "aggregate.csv"

    argv_old = sys.argv[:]
    try:
        sys.argv = ["prog", "--root", str(root), "--out", str(out_csv)]
        aggregate_pipeline_reports.main()
    finally:
        sys.argv = argv_old

    assert out_csv.exists()

    with out_csv.open(newline="") as f:
        r = csv.DictReader(f)
        rows = list(r)

    assert len(rows) == 1
    row = rows[0]

    assert "phase_snps_only" in row
    assert "eval_snps_only" in row
    assert row["phase_snps_only"] not in ("", None)
    assert row["eval_snps_only"] not in ("", None)