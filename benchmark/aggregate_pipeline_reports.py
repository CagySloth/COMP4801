import argparse
import json
import csv
from pathlib import Path
from typing import Any, Dict, Optional


def get(d: Dict[str, Any], path: str, default=None):
    cur: Any = d
    for p in path.split("."):
        if not isinstance(cur, dict) or p not in cur:
            return default
        cur = cur[p]
    return cur


def load_json(p: Path) -> Optional[Dict[str, Any]]:
    try:
        with open(p, "r") as f:
            return json.load(f)
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser(description="Aggregate *.pipeline.json into a CSV.")
    ap.add_argument("--root", default="output/exp_baseline_q20", help="Directory to search for pipeline.json files")
    ap.add_argument("--out", default="output/exp_baseline_q20/aggregate.csv", help="Output CSV path")
    args = ap.parse_args()

    root = Path(args.root)
    files = sorted(root.rglob("*.pipeline.json"))
    if not files:
        raise SystemExit(f"No pipeline json files found under {root}")

    rows = []
    for fp in files:
        j = load_json(fp)
        if not j:
            continue

        row = {
            "pipeline_json": str(fp),

            # params
            "seed": get(j, "seed"),
            "platform": get(j, "params.platform"),
            "ont_profile": get(j, "params.ont_profile"),
            "ref_preset": get(j, "params.ref_preset"),
            "num_reads": get(j, "params.num_reads"),
            "min_len": get(j, "params.min_len"),
            "max_len": get(j, "params.max_len"),
            "max_coverage": get(j, "params.max_coverage"),

            # calling quality (only meaningful for called)
            "call_precision": get(j, "callset.precision"),
            "call_recall": get(j, "callset.recall"),
            "called_snps": get(j, "callset.called_snps"),
            "truth_snps": get(j, "callset.truth_snps"),
            "truth_het_snps": get(j, "callset.truth_het_snps"),

            # oracle phasing (WhatsHap-only)
            "oracle_phase_accuracy": get(j, "phasing_runs.oracle.eval.phase_accuracy_blockflip"),
            "oracle_switch_error": get(j, "phasing_runs.oracle.eval.switch_error_rate"),
            "oracle_num_phase_sets": get(j, "phasing_runs.oracle.eval.num_phase_sets"),
            "oracle_phasing_rate_shared_het": get(j, "phasing_runs.oracle.derived.phasing_rate_on_shared_het"),
            "oracle_effective_phased_recall": get(j, "phasing_runs.oracle.derived.effective_phased_recall"),
            "oracle_shared_het_recall": get(j, "phasing_runs.oracle.derived.shared_het_recall"),

            # called phasing (end-to-end)
            "called_phase_accuracy": get(j, "phasing_runs.called.eval.phase_accuracy_blockflip"),
            "called_switch_error": get(j, "phasing_runs.called.eval.switch_error_rate"),
            "called_num_phase_sets": get(j, "phasing_runs.called.eval.num_phase_sets"),
            "called_phasing_rate_shared_het": get(j, "phasing_runs.called.derived.phasing_rate_on_shared_het"),
            "called_effective_phased_recall": get(j, "phasing_runs.called.derived.effective_phased_recall"),
            "called_shared_het_recall": get(j, "phasing_runs.called.derived.shared_het_recall"),
            
            # duplications
            "dup_segments": get(j, "params.dup_segments"),
            "dup_len": get(j, "params.dup_len"),
            "dup_min_gap": get(j, "params.dup_min_gap"),

            "len_model": get(j, "params.len_model"),
            "ln_mean": get(j, "params.ln_mean"),
            "ln_sigma": get(j, "params.ln_sigma"),
            "start_model": get(j, "params.start_model"),
            "dropout_fraction": get(j, "params.dropout_fraction"),
            "dropout_block_len": get(j, "params.dropout_block_len"),
        }
        rows.append(row)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = list(rows[0].keys())
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    print(f"✅ Wrote {len(rows)} rows to: {out}")


if __name__ == "__main__":
    main()