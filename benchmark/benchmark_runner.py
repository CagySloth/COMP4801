import argparse
import subprocess
import time
import json
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


# def run_algorithm(algo: str, input_npz: Path, output_prefix: Path, ploidy: int = None):
#     cmd = ["python", "-m", "algorithms.cli.phase", algo, "-i", str(input_npz), "-o", str(output_prefix)]
#     if ploidy and "polyploid" in algo.lower():
#         cmd += ["-k", str(ploidy)]
#     start = time.time()
#     try:
#         subprocess.run(cmd, check=True)
#         return time.time() - start
#     except subprocess.CalledProcessError as e:
#         print(f"[ERROR] Algorithm {algo} failed: {e}")
#         return None

def run_accuracy_benchmark(truth_path: Path, pred_path: Path, json_out: Path):
    cmd = [
        "python", "-m", "benchmark.benchmark_accuracy",
        "--truth", str(truth_path),
        "--pred", str(pred_path),
        "--output", str(json_out)
    ]
    subprocess.run(cmd, check=True)

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

            truth = sim_prefix.with_suffix(".haplotypes.tsv")
            sparse_tsv = sim_prefix.with_suffix(".reads.sparse.tsv")
            reads_npz = sim_prefix.with_suffix(".reads.npz")
            vcf_path = sim_prefix.with_suffix(".vcf")

            # subprocess.run([
            #     "python", "-m", "algorithms.cli.convert",
            #     "-i", str(sparse_tsv),
            #     "--to-npz", str(reads_npz)
            # ], check=True)

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

                pred_hap = algo_out.with_suffix(".haplotypes.tsv")
                acc_json = algo_out.with_suffix(".accuracy.json")

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

                if vary_field:
                    record[vary_field] = sweep_args[vary_field]

                summary.append(record)

    out_path = args.outdir / "benchmark_summary.json"
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[DONE] Benchmark summary written to {out_path}")

if __name__ == "__main__":
    main()
