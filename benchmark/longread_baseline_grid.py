import argparse
import subprocess
from pathlib import Path


def run_one(prefix: str, seed: int, num_reads: int, max_cov: int, ont_profile: str) -> None:
    cmd = [
        "python", "-m", "benchmark.longread_pipeline_runner",
        "--prefix", prefix,
        "--seed", str(seed),
        "--vcf-source", "both",
        "--platform", "ont",
        "--ont-profile", ont_profile,
        "--num-reads", str(num_reads),
        "--max-coverage", str(max_cov),
    ]
    print("[RUN]", " ".join(cmd))
    subprocess.check_call(cmd)


def main():
    ap = argparse.ArgumentParser(description="Run a grid of long-read pipeline experiments.")
    ap.add_argument("--outdir", default="output/exp_baseline_q20", help="Directory to store runs")
    ap.add_argument("--ont-profile", default="q20", choices=["q20", "classic"])
    ap.add_argument("--seeds", default="0,1,2")
    ap.add_argument("--num-reads", default="50,100,200,400")
    ap.add_argument("--max-coverage", default="10,15,20")

    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    seeds = [int(x) for x in args.seeds.split(",") if x.strip() != ""]
    reads_list = [int(x) for x in args.num_reads.split(",") if x.strip() != ""]
    max_cov_list = [int(x) for x in args.max_coverage.split(",") if x.strip() != ""]

    for num_reads in reads_list:
        for max_cov in max_cov_list:
            for seed in seeds:
                run_name = f"ont_{args.ont_profile}_r{num_reads}_mc{max_cov}_s{seed}"
                prefix = str(outdir / run_name / run_name)
                Path(prefix).parent.mkdir(parents=True, exist_ok=True)
                run_one(prefix, seed, num_reads, max_cov, args.ont_profile)

    print(f"✅ Completed grid. Results under: {outdir}")


if __name__ == "__main__":
    main()