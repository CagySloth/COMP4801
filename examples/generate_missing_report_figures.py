from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

report_dir = Path("report/figures/10_experiments")
report_dir.mkdir(parents=True, exist_ok=True)

def plot_errorbar(df: pd.DataFrame, xcol: str, ycol: str, out_png: Path, ylabel: str, title: str = None):
    if xcol not in df.columns or ycol not in df.columns:
        print(f"Skipping {out_png.name}: missing {xcol} or {ycol}")
        return
    g = (
        df.groupby(xcol, as_index=False)[ycol]
          .agg(["mean", "std"])
          .reset_index()
          .sort_values(xcol)
    )
    plt.figure(figsize=(5.5, 4.2))
    plt.errorbar(g[xcol], g["mean"], yerr=g["std"], marker="o", capsize=4, linewidth=2)
    plt.xlabel(xcol)
    plt.ylabel(ylabel)
    plt.title(title or ylabel)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_png}")

def plot_heatmap(df: pd.DataFrame, xcol: str, ycol: str, metric: str, out_png: Path, title: str):
    if xcol not in df.columns or ycol not in df.columns or metric not in df.columns:
        print(f"Skipping {out_png.name}: missing {xcol}, {ycol}, or {metric}")
        return

    xs = sorted(df[xcol].unique())
    ys = sorted(df[ycol].unique())
    piv = (
        df.groupby([ycol, xcol])[metric]
          .mean()
          .unstack(xcol)
          .reindex(index=ys, columns=xs)
    )

    plt.figure(figsize=(5.5, 4.2))
    im = plt.imshow(piv.values, aspect="auto", origin="lower")
    plt.xticks(range(len(xs)), xs)
    plt.yticks(range(len(ys)), ys)
    plt.xlabel(xcol)
    plt.ylabel(ycol)
    plt.title(title)
    plt.colorbar(im, label=metric)

    for yi, y in enumerate(ys):
        for xi, x in enumerate(xs):
            val = piv.loc[y, x]
            if pd.notna(val):
                plt.text(xi, yi, f"{val:.3f}", ha="center", va="center", fontsize=8)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_png}")

# 10.5.5 — oracle effective phased recall vs burst probability
csv_105 = Path("output/experiments_draft/04_ablation_bursts/aggregate.csv")
if csv_105.exists():
    df = pd.read_csv(csv_105)
    plot_errorbar(
        df, "burst_prob", "oracle_effective_phased_recall",
        report_dir / "fig_10_5_5_oracle_effective_phased_recall.png",
        ylabel="Oracle effective phased recall",
        title="Oracle effective phased recall vs burst probability"
    )
else:
    print("Missing:", csv_105)

# 10.9.2.3 + 10.9.2.4 — qual threshold grid heatmaps
csv_1092 = Path("output/experiments_draft/08_optimize_whatshap/qual_threshold_grid/aggregate.csv")
if csv_1092.exists():
    df = pd.read_csv(csv_1092)
    plot_heatmap(
        df, "min_mapq", "min_baseq", "called_num_phase_sets",
        report_dir / "fig_10_9_2_3_called_num_phase_sets_heatmap.png",
        title="Called num phase sets"
    )
    plot_heatmap(
        df, "min_mapq", "min_baseq", "oracle_effective_phased_recall",
        report_dir / "fig_10_9_2_4_oracle_effective_phased_recall_heatmap.png",
        title="Oracle effective phased recall"
    )
else:
    print("Missing:", csv_1092)

# 10.9.4.3 — oracle effective phased recall vs max_coverage
csv_1094 = Path("output/experiments_draft/08_optimize_whatshap/maxcov_runtime_sweep/aggregate.csv")
if csv_1094.exists():
    df = pd.read_csv(csv_1094)
    plot_errorbar(
        df, "max_coverage", "oracle_effective_phased_recall",
        report_dir / "fig_10_9_4_3_oracle_effective_phased_recall.png",
        ylabel="Oracle effective phased recall",
        title="Oracle effective phased recall vs max_coverage"
    )
else:
    print("Missing:", csv_1094)

# 10.9.6.3 — called num phase sets across default vs optimized
csv_1096 = Path("output/experiments_draft/08_optimize_whatshap/confirm_default_vs_optimized/aggregate.csv")
if csv_1096.exists():
    df = pd.read_csv(csv_1096)

    # try config_label first, otherwise infer from path
    if "config_label" not in df.columns:
        def infer_config(path: str) -> str:
            p = str(path).lower()
            if "optimized" in p:
                return "optimized"
            if "default" in p:
                return "default"
            return "unknown"
        df["config_label"] = df["pipeline_json"].apply(infer_config)

    g = (
        df.groupby("config_label", as_index=False)["called_num_phase_sets"]
          .agg(["mean", "std"])
          .reset_index()
    )

    plt.figure(figsize=(5.2, 4.0))
    plt.bar(g["config_label"], g["mean"], yerr=g["std"], capsize=4)
    plt.xlabel("Configuration")
    plt.ylabel("Called num phase sets")
    plt.title("Called num phase sets across default vs optimized")
    plt.tight_layout()
    out = report_dir / "fig_10_9_6_3_called_num_phase_sets.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")
else:
    print("Missing:", csv_1096)