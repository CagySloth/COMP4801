from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

csv_path = Path("output/poster_depth_dense/aggregate.csv")
out_path = Path("output/poster_depth_dense/poster_baseline_oracle_vs_called.png")

df = pd.read_csv(csv_path)

print("Using num_reads settings:", sorted(df["num_reads"].unique()))

g = (
    df.groupby("num_reads", as_index=False)
      .agg(
          oracle_mean=("oracle_effective_phased_recall", "mean"),
          oracle_std=("oracle_effective_phased_recall", "std"),
          called_mean=("called_effective_phased_recall", "mean"),
          called_std=("called_effective_phased_recall", "std"),
      )
      .sort_values("num_reads")
)

plt.figure(figsize=(7, 4.5))

plt.errorbar(
    g["num_reads"], g["oracle_mean"], yerr=g["oracle_std"],
    marker="o", capsize=4, linewidth=2, label="Oracle effective phased recall"
)

plt.errorbar(
    g["num_reads"], g["called_mean"], yerr=g["called_std"],
    marker="s", capsize=4, linewidth=2, label="Called effective phased recall"
)

plt.xlabel("Number of reads (coverage proxy)")
plt.ylabel("Effective phased recall")
plt.title("Baseline depth sweep: oracle vs called")
plt.ylim(0, 1.05)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

out_path.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(out_path, dpi=300, bbox_inches="tight")
print(f"Saved to: {out_path}")