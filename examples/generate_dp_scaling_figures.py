from pathlib import Path
import re
import pandas as pd
import matplotlib.pyplot as plt

csv_path = Path("output/experiments_draft_profile/dp_scaling/aggregate.csv")
report_dir = Path("report/figures/10_experiments")
report_dir.mkdir(parents=True, exist_ok=True)

if not csv_path.exists():
    raise FileNotFoundError(f"Missing DP-scaling aggregate: {csv_path}")

df = pd.read_csv(csv_path)

# Parse config and num_snps from pipeline_json path if needed
if "cfg" not in df.columns or "num_snps" not in df.columns:
    def parse(path):
        path = str(path)
        m = re.search(r"/(default|balanced)_snps(\d+)_s(\d+)\.pipeline\.json$", path)
        if not m:
            return pd.Series([None, None, None])
        return pd.Series([m.group(1), int(m.group(2)), int(m.group(3))])
    df[["cfg", "num_snps", "seed"]] = df["pipeline_json"].apply(parse)

label_map = {"default": "Default", "balanced": "Optimized"}
df["cfg_label"] = df["cfg"].map(label_map)

def plot_metric(metric: str, out_name: str, ylabel: str, title: str):
    g = (
        df.groupby(["cfg_label", "num_snps"], as_index=False)[metric]
          .agg(["mean", "std"])
          .reset_index()
    )
    plt.figure(figsize=(5.8, 4.2))
    for cfg in ["Default", "Optimized"]:
        sub = g[g["cfg_label"] == cfg].sort_values("num_snps")
        plt.errorbar(
            sub["num_snps"], sub["mean"], yerr=sub["std"],
            marker="o", capsize=4, linewidth=2, label=cfg
        )
    plt.xlabel("num_snps")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    out = report_dir / out_name
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")

plot_metric(
    "called_time_solve_sec",
    "fig_10_9_11_1_called_time_solve_sec.png",
    "Called solve time (s)",
    "Called solve time vs SNP count"
)
plot_metric(
    "oracle_time_solve_sec",
    "fig_10_9_11_2_oracle_time_solve_sec.png",
    "Oracle solve time (s)",
    "Oracle solve time vs SNP count"
)
plot_metric(
    "called_effective_phased_recall",
    "fig_10_9_11_3_called_effective_phased_recall.png",
    "Called effective phased recall",
    "Called effective phased recall vs SNP count"
)
plot_metric(
    "called_num_phase_sets",
    "fig_10_9_11_4_called_num_phase_sets.png",
    "Called num phase sets",
    "Called num phase sets vs SNP count"
)