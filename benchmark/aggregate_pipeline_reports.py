# benchmark/aggregate_pipeline_reports.py
"""
Like aggregate_pipeline_reports.py, but includes more params and fields useful for experiments.

Usage:
  python -m benchmark.aggregate_pipeline_reports --root output/exp --out output/exp/aggregate.csv
"""
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
    ap = argparse.ArgumentParser(description="Aggregate *.pipeline.json into a CSV (full).")
    ap.add_argument("--root", required=True, help="Directory to search for pipeline.json files")
    ap.add_argument("--out", required=True, help="Output CSV path")
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

            # top-level
            "prefix": get(j, "prefix"),
            "seed": get(j, "seed"),
            "time_total_sec": get(j, "time_total_sec"),

            # common params
            "ref_preset": get(j, "params.ref_preset"),
            "ref_length": get(j, "params.ref_length"),
            "num_snps": get(j, "params.num_snps"),
            "het_rate": get(j, "params.het_rate"),
            "avoid_regions": get(j, "params.avoid_regions"),

            "platform": get(j, "params.platform"),
            "ont_profile": get(j, "params.ont_profile"),
            "map_preset": get(j, "params.map_preset"),

            # reads
            "num_reads": get(j, "params.num_reads"),
            "min_len": get(j, "params.min_len"),
            "max_len": get(j, "params.max_len"),
            "len_model": get(j, "params.len_model"),
            "ln_mean": get(j, "params.ln_mean"),
            "ln_sigma": get(j, "params.ln_sigma"),
            "start_model": get(j, "params.start_model"),
            "dropout_fraction": get(j, "params.dropout_fraction"),
            "dropout_block_len": get(j, "params.dropout_block_len"),

            # realism knobs
            "dup_segments": get(j, "params.dup_segments"),
            "dup_len": get(j, "params.dup_len"),
            "dup_min_gap": get(j, "params.dup_min_gap"),

            "burst_prob": get(j, "params.burst_prob"),
            "burst_count": get(j, "params.burst_count"),
            "burst_len": get(j, "params.burst_len"),
            "burst_mult": get(j, "params.burst_mult"),

            "num_indels": get(j, "params.num_indels"),
            "indel_min_len": get(j, "params.indel_min_len"),
            "indel_max_len": get(j, "params.indel_max_len"),
            "indel_het_rate": get(j, "params.indel_het_rate"),

            "phase_snps_only": get(j, "params.phase_snps_only"),
            "eval_snps_only": get(j, "params.eval_snps_only"),

            # calling knobs
            "call_min_mapq": get(j, "params.call_min_mapq"),
            "call_min_baseq": get(j, "params.call_min_baseq"),

            # phasing knobs
            "max_coverage": get(j, "params.max_coverage"),
            "min_mapq": get(j, "params.min_mapq"),
            "min_baseq": get(j, "params.min_baseq"),

            # callset quality
            "call_precision": get(j, "callset.precision"),
            "call_recall": get(j, "callset.recall"),
            "truth_snps": get(j, "callset.truth_snps"),
            "truth_het_snps": get(j, "callset.truth_het_snps"),
            "called_snps": get(j, "callset.called_snps"),
            "shared_snps": get(j, "callset.shared_snps"),

            # oracle runtime breakdown
            "oracle_time_vcf_parse_sec": get(j, "phasing_runs.oracle.summary.time_vcf_parse_sec"),
            "oracle_time_build_readset_sec": get(j, "phasing_runs.oracle.summary.time_build_readset_sec"),
            "oracle_time_readselection_sec": get(j, "phasing_runs.oracle.summary.time_readselection_sec"),
            "oracle_time_solve_sec": get(j, "phasing_runs.oracle.summary.time_solve_sec"),
            "oracle_time_ps_sec": get(j, "phasing_runs.oracle.summary.time_ps_sec"),
            "oracle_time_write_sec": get(j, "phasing_runs.oracle.summary.time_write_sec"),
            "oracle_time_total_sec": get(j, "phasing_runs.oracle.summary.time_total_sec"),

            # called runtime breakdown
            "called_time_vcf_parse_sec": get(j, "phasing_runs.called.summary.time_vcf_parse_sec"),
            "called_time_build_readset_sec": get(j, "phasing_runs.called.summary.time_build_readset_sec"),
            "called_time_readselection_sec": get(j, "phasing_runs.called.summary.time_readselection_sec"),
            "called_time_solve_sec": get(j, "phasing_runs.called.summary.time_solve_sec"),
            "called_time_ps_sec": get(j, "phasing_runs.called.summary.time_ps_sec"),
            "called_time_write_sec": get(j, "phasing_runs.called.summary.time_write_sec"),
            "called_time_total_sec": get(j, "phasing_runs.called.summary.time_total_sec"),

            # oracle eval
            "oracle_phase_accuracy": get(j, "phasing_runs.oracle.eval.phase_accuracy_blockflip"),
            "oracle_switch_error": get(j, "phasing_runs.oracle.eval.switch_error_rate"),
            "oracle_num_phase_sets": get(j, "phasing_runs.oracle.eval.num_phase_sets"),
            "oracle_effective_phased_recall": get(j, "phasing_runs.oracle.derived.effective_phased_recall"),
            "oracle_phasing_rate_shared_het": get(j, "phasing_runs.oracle.derived.phasing_rate_on_shared_het"),
            "oracle_shared_het_recall": get(j, "phasing_runs.oracle.derived.shared_het_recall"),

            # called eval
            "called_phase_accuracy": get(j, "phasing_runs.called.eval.phase_accuracy_blockflip"),
            "called_switch_error": get(j, "phasing_runs.called.eval.switch_error_rate"),
            "called_num_phase_sets": get(j, "phasing_runs.called.eval.num_phase_sets"),
            "called_effective_phased_recall": get(j, "phasing_runs.called.derived.effective_phased_recall"),
            "called_phasing_rate_shared_het": get(j, "phasing_runs.called.derived.phasing_rate_on_shared_het"),
            "called_shared_het_recall": get(j, "phasing_runs.called.derived.shared_het_recall"),
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
