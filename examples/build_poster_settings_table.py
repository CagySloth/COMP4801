from pathlib import Path
import re
import pandas as pd

transfer_csv = Path("output/experiments_draft/08_optimize_whatshap/transferability_across_scenarios/aggregate.csv")
accuracy_csv = Path("output/poster_accuracy_row/aggregate.csv")  # optional

out_csv = Path("output/poster_settings_table.csv")

def normalize_transfer(df: pd.DataFrame) -> pd.DataFrame:
    # Prefer explicit labels if present
    if "scenario_label" in df.columns and "config_label" in df.columns:
        x = df.copy()
        x["scenario"] = x["scenario_label"].str.lower()
        x["config"] = x["config_label"].str.lower()
        return x

    # Fallback: parse from pipeline_json path
    x = df.copy()
    def parse(path):
        path = str(path).lower()
        cfg = None
        scen = None
        for c in ["default", "balanced", "runtime"]:
            if c in path:
                cfg = c
                break
        for s in ["baseline", "dropout", "interaction", "hard"]:
            if s in path:
                scen = s
                break
        return pd.Series([cfg, scen])

    x[["config", "scenario"]] = x["pipeline_json"].apply(parse)
    return x

def normalize_accuracy(df: pd.DataFrame) -> pd.DataFrame:
    x = df.copy()
    x["config"] = "accuracy-oriented"

    def parse(path):
        path = str(path).lower()
        if "baseline_acc" in path:
            return "baseline"
        if "hard_acc" in path:
            return "hard"
        return None

    x["scenario"] = x["pipeline_json"].apply(parse)
    return x

# Load main transferability data
df_transfer = pd.read_csv(transfer_csv)
df_transfer = normalize_transfer(df_transfer)

frames = [df_transfer]

# Add accuracy-oriented row if file exists
if accuracy_csv.exists():
    df_acc = pd.read_csv(accuracy_csv)
    df_acc = normalize_accuracy(df_acc)
    frames.append(df_acc)

df = pd.concat(frames, ignore_index=True)

# Keep only poster scenarios
df = df[df["scenario"].isin(["baseline", "hard"])].copy()

# Metrics to keep
metrics = [
    "called_effective_phased_recall",
    "called_num_phase_sets",
    "time_total_sec",
]

# Mean across seeds
g = (
    df.groupby(["config", "scenario"], as_index=False)[metrics]
      .mean()
)

# Friendly labels
config_map = {
    "default": "Default",
    "balanced": "Optimized",
    "runtime": "Speed-oriented",
    "accuracy-oriented": "Accuracy-oriented",
}
scenario_map = {
    "baseline": "Normal",
    "hard": "Hard",
}

g["config"] = g["config"].map(config_map)
g["scenario"] = g["scenario"].map(scenario_map)

# Pivot to wide table
rows = []
for config in ["Default", "Optimized", "Speed-oriented", "Accuracy-oriented"]:
    sub = g[g["config"] == config]
    if sub.empty:
        continue

    row = {"Configuration": config}
    for _, r in sub.iterrows():
        scen = r["scenario"]
        row[f"{scen}: Eff. Recall"] = round(r["called_effective_phased_recall"], 4)
        row[f"{scen}: Phase Sets"] = round(r["called_num_phase_sets"], 1)
        row[f"{scen}: Runtime (s)"] = round(r["time_total_sec"], 4)
    rows.append(row)

table = pd.DataFrame(rows)

# Column order
wanted_cols = [
    "Configuration",
    "Normal: Eff. Recall", "Normal: Phase Sets", "Normal: Runtime (s)",
    "Hard: Eff. Recall", "Hard: Phase Sets", "Hard: Runtime (s)",
]
table = table[[c for c in wanted_cols if c in table.columns]]

out_csv.parent.mkdir(parents=True, exist_ok=True)
table.to_csv(out_csv, index=False)

print("\nPoster settings table:\n")
print(table.to_string(index=False))
print(f"\nSaved to: {out_csv}")