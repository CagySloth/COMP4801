import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt


def _mean_std(df, value_col):
    """Return a dataframe with mean and std of value_col grouped by ref_preset, num_reads."""
    g = df.groupby(["ref_preset", "num_reads"])[value_col].agg(["mean", "std"]).reset_index()
    return g


def _plot_errorbar(g, value_col, ylabel, title, outpath):
    """
    g must have columns: ref_preset, num_reads, mean, std
    """
    plt.figure()
    for preset in sorted(g["ref_preset"].unique()):
        sub = g[g["ref_preset"] == preset].sort_values("num_reads")
        plt.errorbar(
            sub["num_reads"],
            sub["mean"],
            yerr=sub["std"].fillna(0.0),
            marker="o",
            capsize=3,
            label=preset,
        )
    plt.xlabel("num_reads")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close()


def main():
    ap = argparse.ArgumentParser(description="Plot realism comparison (plain vs realistic) from aggregate.csv.")
    ap.add_argument("--csv", default="output/exp_refpreset_compare/aggregate.csv")
    ap.add_argument("--outdir", default="output/exp_refpreset_compare/plots")
    ap.add_argument("--max-coverage", type=int, default=None, help="Filter to a specific max_coverage (e.g., 15).")
    args = ap.parse_args()

    df = pd.read_csv(args.csv)

    # Optional filter: max_coverage (if your CSV includes multiple values)
    if args.max_coverage is not None and "max_coverage" in df.columns:
        df = df[df["max_coverage"] == args.max_coverage]

    os.makedirs(args.outdir, exist_ok=True)

    # --- Plot 1: Calling recall ---
    if "call_recall" in df.columns:
        g = _mean_std(df, "call_recall")
        _plot_errorbar(
            g,
            "call_recall",
            ylabel="Variant calling recall (truth SNPs)",
            title="Calling recall vs num_reads (plain vs realistic)",
            outpath=os.path.join(args.outdir, "calling_recall_plain_vs_realistic.png"),
        )

    # --- Plot 2: Effective phased recall (oracle vs called) ---
    for col, label, fname in [
        ("oracle_effective_phased_recall", "Oracle effective phased recall", "oracle_effective_phased_recall.png"),
        ("called_effective_phased_recall", "Called effective phased recall", "called_effective_phased_recall.png"),
    ]:
        if col in df.columns:
            g = _mean_std(df, col)
            _plot_errorbar(
                g,
                col,
                ylabel=label,
                title=f"{label} vs num_reads (plain vs realistic)",
                outpath=os.path.join(args.outdir, fname),
            )

    # --- Plot 3: Phase set fragmentation ---
    # (More phase sets => more fragmentation)
    for col, label, fname in [
        ("oracle_num_phase_sets", "Oracle num_phase_sets", "oracle_num_phase_sets.png"),
        ("called_num_phase_sets", "Called num_phase_sets", "called_num_phase_sets.png"),
    ]:
        if col in df.columns:
            g = _mean_std(df, col)
            _plot_errorbar(
                g,
                col,
                ylabel=label,
                title=f"{label} vs num_reads (plain vs realistic)",
                outpath=os.path.join(args.outdir, fname),
            )

    # --- Plot 4 (optional but great): Shared het recall (called) ---
    if "called_shared_het_recall" in df.columns:
        g = _mean_std(df, "called_shared_het_recall")
        _plot_errorbar(
            g,
            "called_shared_het_recall",
            ylabel="Shared het recall (called VCF ∩ truth het) / truth het",
            title="Called shared het recall vs num_reads (plain vs realistic)",
            outpath=os.path.join(args.outdir, "called_shared_het_recall.png"),
        )

    print(f"✅ Plots written to: {args.outdir}")


if __name__ == "__main__":
    main()