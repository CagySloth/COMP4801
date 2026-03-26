import argparse
import pandas as pd
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="output/exp_baseline_q20/aggregate.csv")
    ap.add_argument("--outdir", default="output/exp_baseline_q20/plots")
    args = ap.parse_args()

    df = pd.read_csv(args.csv)
    outdir = args.outdir
    import os
    os.makedirs(outdir, exist_ok=True)

    # Group by num_reads and max_coverage, average over seeds
    grp = df.groupby(["num_reads", "max_coverage"], as_index=False).mean(numeric_only=True)

    # Plot 1: Calling recall vs num_reads (for each max_coverage)
    plt.figure()
    for mc in sorted(grp["max_coverage"].unique()):
        sub = grp[grp["max_coverage"] == mc].sort_values("num_reads")
        plt.plot(sub["num_reads"], sub["call_recall"], marker="o", label=f"max_cov={mc}")
    plt.xlabel("num_reads")
    plt.ylabel("calling recall (truth SNPs)")
    plt.legend()
    plt.savefig(f"{outdir}/calling_recall_vs_num_reads.png", dpi=200, bbox_inches="tight")
    plt.close()

    # Plot 2: Effective phased recall (oracle vs called) vs num_reads, mc=15 (default slice)
    mc_pick = 15
    sub = grp[grp["max_coverage"] == mc_pick].sort_values("num_reads")
    plt.figure()
    plt.plot(sub["num_reads"], sub["oracle_effective_phased_recall"], marker="o", label="oracle (phasing-only)")
    plt.plot(sub["num_reads"], sub["called_effective_phased_recall"], marker="o", label="called (end-to-end)")
    plt.xlabel("num_reads")
    plt.ylabel("effective phased recall (vs truth het SNPs)")
    plt.legend()
    plt.savefig(f"{outdir}/effective_phased_recall_vs_num_reads_mc{mc_pick}.png", dpi=200, bbox_inches="tight")
    plt.close()

    # Plot 3: Phase set fragmentation (num_phase_sets) vs num_reads, oracle vs called for mc=15
    plt.figure()
    plt.plot(sub["num_reads"], sub["oracle_num_phase_sets"], marker="o", label="oracle num_phase_sets")
    plt.plot(sub["num_reads"], sub["called_num_phase_sets"], marker="o", label="called num_phase_sets")
    plt.xlabel("num_reads")
    plt.ylabel("num_phase_sets (fragmentation)")
    plt.legend()
    plt.savefig(f"{outdir}/phase_sets_vs_num_reads_mc{mc_pick}.png", dpi=200, bbox_inches="tight")
    plt.close()

    print(f"✅ Wrote plots to: {outdir}")


if __name__ == "__main__":
    main()