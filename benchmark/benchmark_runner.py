import argparse
import subprocess
import time
import json
import csv
from pathlib import Path

def run_simulator(ploidy: int, num_variants: int, num_reads: int, read_length: int,
                  error_rate: float, missing_rate: float, alpha: float, beta: float,
                  allow_mono: bool, seed: int, out_prefix: Path):
    
    cmd = [
        "python", "-m", "dataset.simulate",
        "-p", str(ploidy),
        "-n", str(num_variants),
        "-r", str(num_reads),
        "-l", str(read_length),
        "-e", str(error_rate),
        "-m", str(missing_rate),
        "--maf-alpha", str(alpha),
        "--maf-beta", str(beta),
        "--seed", str(seed),
        "-o", str(out_prefix),
    ]

    if allow_mono:
        cmd.append("--allow-monomorphic")

    subprocess.run(cmd, check=True)

def run_algorithm(algo: str, input_npz: Path, output_prefix: Path, ploidy: int = None, vcf_path: Path = None):
    cmd = [
        "python", "-m", "algorithms.cli.phase",
        algo,
        "-i", str(input_npz),
        "--output-prefix", str(output_prefix),
    ]

    # polyploid subcommands accept --ploidy
    if ploidy is not None and algo.startswith("polyploid-"):
        cmd += ["--ploidy", str(ploidy)]

    # WhatsHap VCF-mode
    if vcf_path is not None and algo == "diploid-whats":
        cmd += ["--vcf", str(vcf_path)]

    start = time.perf_counter()
    try:
        subprocess.run(cmd, check=True)
        return time.perf_counter() - start
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Algorithm {algo} failed: {e}")
        return None

def run_accuracy_benchmark(truth_path: Path, pred_path: Path, json_out: Path):
    cmd = [
        "python", "-m", "benchmark.benchmark_accuracy",
        "--truth", str(truth_path),
        "--pred", str(pred_path),
        "--output", str(json_out)
    ]
    subprocess.run(cmd, check=True)
    
def write_summary_csv(records: list[dict], out_csv: Path) -> None:
    # Expand label_permutation (list) into separate columns for easy plotting
    max_perm = 0
    for r in records:
        perm = r.get("label_permutation")
        if isinstance(perm, list):
            max_perm = max(max_perm, len(perm))

    # Collect all keys except label_permutation (handled separately)
    keys = set()
    for r in records:
        keys.update(r.keys())
    keys.discard("label_permutation")

    # Stable column order: put common ones first, then the rest
    preferred = [
        "run", "algorithm",
        "ploidy", "num_variants", "num_reads", "read_length",
        "error_rate", "missing_rate",
        "runtime_seconds", "accuracy", "note",
        # timing breakdown columns we may add (see next section)
        "phase_solver",
        "phase_time_total_sec", "phase_time_build_readset_sec",
        "phase_time_readselection_sec", "phase_time_solve_sec",
        "phase_selected_reads", "phase_num_phase_sets",
    ]
    cols = [c for c in preferred if c in keys]
    cols += sorted(keys - set(cols))

    # Add label permutation columns
    cols += [f"label_perm_{i}" for i in range(max_perm)]

    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in records:
            row = {k: r.get(k) for k in cols if not k.startswith("label_perm_")}
            perm = r.get("label_permutation")
            if isinstance(perm, list):
                for i, v in enumerate(perm):
                    row[f"label_perm_{i}"] = v
            w.writerow(row)


def parse_args():
    parser = argparse.ArgumentParser(description="Benchmark phasing algorithms with synthetic datasets.")
    parser.add_argument("--algorithms", nargs="+", required=True, help="Algorithms to test (e.g. em, mst)")
    parser.add_argument("--ploidy", type=int, default=2)
    parser.add_argument("--num-variants", type=int, default=1000)
    parser.add_argument("--num-reads", type=int, default=5000)
    parser.add_argument("--read-length", type=int, default=50)
    parser.add_argument("--error-rate", type=float, default=0.01)
    parser.add_argument("--missing-rate", type=float, default=0.05)
    parser.add_argument("--maf-alpha", type=float, default=0.4)
    parser.add_argument("--maf-beta", type=float, default=0.4)
    parser.add_argument("--allow-monomorphic", action="store_true")
    parser.add_argument("--num-runs", type=int, default=1)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--vary", type=str, default=None, help="Parameter to vary (e.g. num_reads)")
    parser.add_argument("--vary-values", nargs="+", default=None, help="Values to sweep for --vary")
    return parser.parse_args()

def add_ext(prefix: Path, ext: str) -> Path:
    # Safe even if prefix name contains dots (e.g., 0.005)
    return prefix.parent / (prefix.name + ext)

def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    summary = []

    vary_field = args.vary
    vary_values = args.vary_values or [None]

    for vary_val in vary_values:
        sweep_args = vars(args).copy()
        if vary_field:
            field_type = type(getattr(args, vary_field))
            sweep_args[vary_field] = (
                vary_val.lower() == "true" if field_type == bool else field_type(vary_val)
            )

        for run_id in range(1, args.num_runs + 1):
            if vary_field:
                sim_prefix = args.outdir / f"sim_{vary_field}{vary_val}_run{run_id}"
            else:
                sim_prefix = args.outdir / f"sim_run{run_id}"
            run_simulator(
                ploidy=sweep_args["ploidy"],
                num_variants=sweep_args["num_variants"],
                num_reads=sweep_args["num_reads"],
                read_length=sweep_args["read_length"],
                error_rate=sweep_args["error_rate"],
                missing_rate=sweep_args["missing_rate"],
                alpha=sweep_args["maf_alpha"],
                beta=sweep_args["maf_beta"],
                allow_mono=sweep_args["allow_monomorphic"],
                seed=42 + run_id,
                out_prefix=sim_prefix
            )

            truth = add_ext(sim_prefix, ".haplotypes.tsv")
            sparse_tsv = add_ext(sim_prefix, ".reads.sparse.tsv")
            reads_npz = add_ext(sim_prefix, ".reads.npz")
            vcf_path = add_ext(sim_prefix, ".vcf")

            for algo in args.algorithms:
                print(f"[INFO] Run {run_id} | {vary_field}={vary_val} | Algorithm: {algo}")
                algo_out = sim_prefix.parent / f"{sim_prefix.name}_{algo}"
                vcf_for_algo = vcf_path if (algo == "diploid-whats" and vcf_path.exists()) else None

                duration = run_algorithm(
                    algo,
                    reads_npz,
                    algo_out,
                    ploidy=sweep_args["ploidy"],
                    vcf_path=vcf_for_algo,
                )
                
                algo_summary_path = add_ext(algo_out, ".summary.json")
                algo_summary = None
                if algo_summary_path.exists():
                    with open(algo_summary_path) as f:
                        algo_summary = json.load(f)


                pred_hap = add_ext(algo_out, ".haplotypes.tsv")
                acc_json = add_ext(algo_out, ".accuracy.json")

                if pred_hap.exists():
                    run_accuracy_benchmark(truth, pred_hap, acc_json)
                    with open(acc_json) as f:
                        acc = json.load(f)
                else:
                    acc = {"error": "haplotypes file missing"}

                record = {
                    "run": run_id,
                    "algorithm": algo,
                    "ploidy": sweep_args["ploidy"],
                    "num_reads": sweep_args["num_reads"],
                    "num_variants": sweep_args["num_variants"],
                    "read_length": sweep_args["read_length"],
                    "error_rate": sweep_args["error_rate"],
                    "missing_rate": sweep_args["missing_rate"],
                    "runtime_seconds": duration,
                    "accuracy": acc.get("accuracy"),
                    "label_permutation": acc.get("label_permutation"),
                    "note": acc.get("error"),
                }
                
                if algo_summary:
                    # Keep these names stable for Excel/R plots
                    record["phase_solver"] = algo_summary.get("solver")
                    record["phase_selected_reads"] = algo_summary.get("selected_reads")
                    record["phase_num_phase_sets"] = algo_summary.get("num_phase_sets")

                    record["phase_time_total_sec"] = algo_summary.get("time_total_sec")
                    record["phase_time_build_readset_sec"] = algo_summary.get("time_build_readset_sec")
                    record["phase_time_readselection_sec"] = algo_summary.get("time_readselection_sec")
                    record["phase_time_solve_sec"] = algo_summary.get("time_solve_sec")

                if vary_field:
                    record[vary_field] = sweep_args[vary_field]

                summary.append(record)

    out_path = args.outdir / "benchmark_summary.json"
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[DONE] Benchmark summary written to {out_path}")
    
    out_csv = args.outdir / "benchmark_summary.csv"
    write_summary_csv(summary, out_csv)
    print(f"[DONE] Benchmark CSV written to {out_csv}")


if __name__ == "__main__":
    main()
