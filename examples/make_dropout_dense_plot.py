from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

csv_path = Path("output/poster_dropout_dense/aggregate.csv")
out_path = Path("output/poster_dropout_dense/poster_dropout_called_effective_phased_recall_2.png")

df = pd.read_csv(csv_path)

print("Using dropout_fraction settings:", sorted(df["dropout_fraction"].unique()))

g = (
    df.groupby("dropout_fraction", as_index=False)
      .agg(
          mean_called=("called_effective_phased_recall", "mean"),
          std_called=("called_effective_phased_recall", "std"),
      )
      .sort_values("dropout_fraction")
)

plt.figure(figsize=(7, 4.5))

plt.errorbar(
    g["dropout_fraction"],
    g["mean_called"],
    yerr=g["std_called"],
    marker="o",
    markersize=8,
    capsize=5,
    linewidth=3,
    label="Called effective phased recall"
)

plt.xlabel("Dropout fraction")
plt.ylabel("Called effective phased recall")
plt.title("Coverage dropout effect on end-to-end phased recall")
plt.ylim(0, 1.05)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

out_path.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(out_path, dpi=300, bbox_inches="tight")
print(f"Saved to: {out_path}")