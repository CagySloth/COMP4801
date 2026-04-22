from pathlib import Path
import pandas as pd

csv_path = Path("output/experiments_draft_profile/08_optimize_whatshap/frontier_config_comparison/aggregate.csv")
df = pd.read_csv(csv_path)

keep = ["default", "balanced", "runtime"]
df = df[df["config_label"].isin(keep)].copy()

cols = [
    "called_time_build_readset_sec",
    "called_time_readselection_sec",
    "called_time_solve_sec",
    "called_time_total_sec",
]

g = df.groupby("config_label", as_index=False)[cols].mean()

for c in [
    "called_time_build_readset_sec",
    "called_time_readselection_sec",
    "called_time_solve_sec",
]:
    g[c + "_pct"] = 100 * g[c] / g["called_time_total_sec"]

label_map = {
    "default": "Default",
    "balanced": "Optimized",
    "runtime": "Speed-oriented",
}
g["Configuration"] = g["config_label"].map(label_map)

out = g[
    [
        "Configuration",
        "called_time_build_readset_sec_pct",
        "called_time_readselection_sec_pct",
        "called_time_solve_sec_pct",
    ]
].round(1)

out.columns = [
    "Configuration",
    "Build ReadSet (%)",
    "Read Selection (%)",
    "Solve (%)",
]

print(out.to_string(index=False))