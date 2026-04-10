from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

csv_path = Path("output/poster_interaction_dense/aggregate.csv")
out_path = Path("output/poster_interaction_dense/poster_interaction_dup_vs_dropout.png")

df = pd.read_csv(csv_path)

print("Using dup_segments settings:", sorted(df["dup_segments"].unique()))
print("Using dropout_fraction settings:", sorted(df["dropout_fraction"].unique()))

g = (
    df.groupby(["dup_segments", "dropout_fraction"], as_index=False)
      .agg(
          called_mean=("called_effective_phased_recall", "mean"),
          called_std=("called_effective_phased_recall", "std"),
      )
      .sort_values(["dup_segments", "dropout_fraction"])
)

g0 = g[g["dup_segments"] == 0]
g5 = g[g["dup_segments"] == 5]

plt.figure(figsize=(7, 4.5))

plt.errorbar(
    g0["dropout_fraction"], g0["called_mean"], yerr=g0["called_std"],
    marker="o", capsize=4, linewidth=2, label="dup_segments = 0"
)

plt.errorbar(
    g5["dropout_fraction"], g5["called_mean"], yerr=g5["called_std"],
    marker="s", capsize=4, linewidth=2, label="dup_segments = 5"
)

plt.xlabel("dropout_fraction")
plt.ylabel("Called effective phased recall")
plt.title("Duplication × dropout interaction")
plt.ylim(0, 0.6)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

out_path.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(out_path, dpi=300, bbox_inches="tight")
print(f"Saved to: {out_path}")