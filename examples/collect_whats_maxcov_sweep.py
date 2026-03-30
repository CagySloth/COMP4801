#!/usr/bin/env python3
import argparse
import csv
import glob
import json
import os
import re
from typing import Any, Dict, Optional


def load_json(path: str) -> Dict[str, Any]:
    with open(path, "r") as f:
        return json.load(f)


def get_first(d: Dict[str, Any], *keys: str, default=None):
    for k in keys:
        if k in d:
            return d[k]
    return default


def parse_maxcov_from_path(p: str) -> Optional[int]:
    m = re.search(r"wh_maxcov(\d+)\.summary\.json$", os.path.basename(p))
    return int(m.group(1)) if m else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--dir",
        default="output/prelim/whats_maxcov_sweep_vcf",
        help="Directory containing wh_maxcov*.summary.json and wh_maxcov*.accuracy.json",
    )
    ap.add_argument(
        "--out",
        default=None,
        help="Output CSV path (default: <dir>/whats_maxcov_sweep.csv)",
    )
    args = ap.parse_args()

    base_dir = args.dir
    out_csv = args.out or os.path.join(base_dir, "whats_maxcov_sweep.csv")

    summary_files = sorted(glob.glob(os.path.join(base_dir, "wh_maxcov*.summary.json")))
    if not summary_files:
        raise SystemExit(f"No summary files found in {base_dir}: expected wh_maxcov*.summary.json")

    rows = []
    for sp in summary_files:
        mc = parse_maxcov_from_path(sp)
        if mc is None:
            continue

        base = sp.replace(".summary.json", "")
        acc_path = base + ".accuracy.json"
        if not os.path.exists(acc_path):
            raise SystemExit(f"Missing accuracy file for {sp}: expected {acc_path}")

        s = load_json(sp)
        a = load_json(acc_path)

        rows.append({
            "max_coverage": mc,
            "solver": get_first(s, "solver", default=None),

            # dataset sizes (if available)
            "N_total_variants": get_first(s, "N", "num_variants", default=None),
            "N_hets_phased": get_first(s, "num_hets", "num_het_sites", "N_hets", default=None),
            "R_total": get_first(s, "num_reads_total", "R", default=None),
            "R_informative": get_first(s, "num_informative_reads", default=None),
            "R_selected": get_first(s, "selected_reads", default=None),

            # timings (tolerate legacy names)
            "time_total_sec": get_first(s, "time_total_sec", "phase_time_total_sec", default=None),
            "time_build_readset_sec": get_first(s, "time_build_readset_sec", "phase_time_build_readset_sec", default=None),
            "time_readselection_sec": get_first(s, "time_readselection_sec", "phase_time_readselection_sec", default=None),
            "time_solve_sec": get_first(s, "time_solve_sec", "phase_time_solve_sec", default=None),

            # accuracy
            "accuracy": a.get("accuracy", None),
        })

    rows.sort(key=lambda r: r["max_coverage"])

    # Choose a stable column order
    fieldnames = [
        "max_coverage", "solver",
        "N_total_variants", "N_hets_phased",
        "R_total", "R_informative", "R_selected",
        "time_total_sec", "time_build_readset_sec", "time_readselection_sec", "time_solve_sec",
        "accuracy",
    ]

    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote: {out_csv}")
    print(f"Rows: {len(rows)}")


if __name__ == "__main__":
    main()
