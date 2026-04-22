from pathlib import Path
import pandas as pd

csv_path = Path("output/experiments_draft_profile/08_optimize_whatshap/frontier_config_comparison/aggregate.csv")
out_csv = Path("output/runtime_breakdown_table.csv")

df = pd.read_csv(csv_path)

required = [
    "config_label",
    "called_time_build_readset_sec",
    "called_time_readselection_sec",
    "called_time_solve_sec",
    "called_time_total_sec",
]
missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}")

# Keep only the three poster/report configurations
keep = ["default", "balanced", "runtime"]
df = df[df["config_label"].isin(keep)].copy()

g = (
    df.groupby("config_label", as_index=False)[
        [
            "called_time_build_readset_sec",
            "called_time_readselection_sec",
            "called_time_solve_sec",
            "called_time_total_sec",
        ]
    ]
    .mean()
)

label_map = {
    "default": "Default",
    "balanced": "Optimized",
    "runtime": "Speed-oriented",
}
g["Configuration"] = g["config_label"].map(label_map)

table = g[
    [
        "Configuration",
        "called_time_build_readset_sec",
        "called_time_readselection_sec",
        "called_time_solve_sec",
        "called_time_total_sec",
    ]
].copy()

table.columns = [
    "Configuration",
    "Build ReadSet (s)",
    "Read Selection (s)",
    "Solve (s)",
    "Called Phasing Total (s)",
]

table = table.round(4)
table.to_csv(out_csv, index=False)

print("\nRuntime breakdown table:\n")
print(table.to_string(index=False))
print(f"\nSaved to: {out_csv}")

# Also print markdown-friendly rows
print("\nMarkdown table:\n")
print("| Configuration | Build ReadSet (s) | Read Selection (s) | Solve (s) | Called Phasing Total (s) |")
print("|---|---:|---:|---:|---:|")
for _, r in table.iterrows():
    print(
        f"| {r['Configuration']} | "
        f"{r['Build ReadSet (s)']:.4f} | "
        f"{r['Read Selection (s)']:.4f} | "
        f"{r['Solve (s)']:.4f} | "
        f"{r['Called Phasing Total (s)']:.4f} |"
    )