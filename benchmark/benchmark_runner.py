import argparse
import subprocess
import time
import json
from pathlib import Path
from typing import List
import shutil
import uuid

def run_simulator(ploidy: int, num_variants: int, num_reads: int, read_length: int,
                  error_rate: float, missing_rate: float, alpha: float, beta: float,
                  allow_mono: bool, seed: int, out_prefix: Path):
    if ploidy == 2:
        sim_script = Path("dataset/0/simulate_phasing_data.py")
    else:
        sim_script = Path("dataset/1/simulate_phasing_data.py")

    cmd = ["python", str(sim_script),
           "-p", str(ploidy),
           "-n", str(num_variants),
           "-r", str(num_reads),
           "-l", str(read_length),
           "-e", str(error_rate),
           "-m", str(missing_rate),
           "--maf-alpha", str(alpha),
           "--maf-beta", str(beta),
           "-s", str(seed),
           "-o", str(out_prefix)]

    if allow_mono:
        cmd.append("--allow-monomorphic")

    subprocess.run(cmd, check=True)

def run_algorithm(algo: str, input_npz: Path, output_prefix: Path, ploidy: int = None):
    cmd = ["python", "-m", "algorithms.cli.phase", algo, "-i", str(input_npz), "-o", str(output_prefix)]
    if ploidy and "polyploid" in algo:
        cmd += ["-k", str(ploidy)]
    start = time.time()
    try:
        subprocess.run(cmd, check=True)
        duration = time.time() - start
        return duration
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Algorithm {algo} failed with error: {e}")
        return None

def run_accuracy_benchmark(truth_path: Path, pred_path: Path, json_out: Path):
    cmd = [
        "python", "benchmarking/benchmark_accuracy.py",
        "--truth", str(truth_path),
        "--pred", str(pred_path),
        "--output", str(json_out)
    ]
    subprocess.run(cmd, check=True)

def parse_args():
    parser = argparse.ArgumentParser(description="Benchmark phasing algorithm under synthetic simulations.")
    parser.add_argument("--algorithms", nargs="+", required=True, help="List of algorithms to run")
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
    parser.add_argument("--vary", type=str, default=None, help="Parameter to sweep (e.g. num_reads, error_rate)")
    parser.add_argument("--vary-values", nargs="+", default=None, help="List of values to sweep for the chosen vary param")
    return parser.parse_args()

def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    summary = []

    vary_field = args.vary
    vary_values = args.vary_values or [None]

    for vary_val in vary_values:
        sweep_kwargs = vars(args).copy()
        if vary_field:
            field_type = type(getattr(args, vary_field))
            if field_type == bool:
                sweep_kwargs[vary_field] = (vary_val.lower() == "true")
            else:
                sweep_kwargs[vary_field] = field_type(vary_val)

        for run_id in range(1, args.num_runs + 1):
            sim_prefix = args.outdir / f"sim_{vary_field}{vary_val}_run{run_id}"
            run_simulator(
                ploidy=sweep_kwargs["ploidy"],
                num_variants=sweep_kwargs["num_variants"],
                num_reads=sweep_kwargs["num_reads"],
                read_length=sweep_kwargs["read_length"],
                error_rate=sweep_kwargs["error_rate"],
                missing_rate=sweep_kwargs["missing_rate"],
                alpha=sweep_kwargs["maf_alpha"],
                beta=sweep_kwargs["maf_beta"],
                allow_mono=sweep_kwargs["allow_monomorphic"],
                seed=42 + run_id,
                out_prefix=sim_prefix
            )

            truth = sim_prefix.with_suffix(".haplotypes.tsv")
            sparse_tsv = sim_prefix.with_suffix(".reads.sparse.tsv")
            reads_npz = sim_prefix.with_suffix(".reads.npz")

            subprocess.run([
                "python", "algorithms/read_tsv.py",
                "-i", str(sparse_tsv),
                "--to-npz", str(reads_npz)
            ], check=True)

            for algo in args.algorithms:
                print(f"[INFO] Run {run_id}, {vary_field}={vary_val}, algorithm: {algo}")
                algo_out = sim_prefix.parent / f"{sim_prefix.name}_{algo}"
                duration = run_algorithm(algo, reads_npz, algo_out, ploidy=sweep_kwargs["ploidy"])

                pred_hap = algo_out.with_suffix(".haplotypes.tsv")
                acc_json = algo_out.with_suffix(".accuracy.json")

                if pred_hap.exists():
                    run_accuracy_benchmark(truth, pred_hap, acc_json)
                    with open(acc_json) as f:
                        acc = json.load(f)
                else:
                    acc = {"error": "haplotypes file missing"}

                summary.append({
                    "run": run_id,
                    "algorithm": algo,
                    "ploidy": sweep_kwargs["ploidy"],
                    "num_reads": sweep_kwargs["num_reads"],
                    "num_variants": sweep_kwargs["num_variants"],
                    "read_length": sweep_kwargs["read_length"],
                    "error_rate": sweep_kwargs["error_rate"],
                    "missing_rate": sweep_kwargs["missing_rate"],
                    "runtime_seconds": duration,
                    "accuracy": acc.get("accuracy"),
                    "label_permutation": acc.get("label_permutation"),
                    "note": acc.get("error"),
                    vary_field: vary_val
                })

    with open(args.outdir / "benchmark_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[DONE] Benchmark summary written to {args.outdir}/benchmark_summary.json")


if __name__ == "__main__":
    main()
