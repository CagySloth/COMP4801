import pandas as pd

df = pd.read_csv(
    "output/experiments_draft_profile/08_optimize_whatshap/frontier_config_comparison/aggregate.csv"
)

cols = [
    "called_time_build_readset_sec",
    "called_time_readselection_sec",
    "called_time_solve_sec",
    "called_time_ps_sec",
    "called_time_write_sec",
    "called_time_total_sec",
]

pd.set_option("display.max_columns", None)
pd.set_option("display.width", 200)

g = df.groupby("config_label")[cols].mean().copy()

# Residual runtime after excluding adaptor/readset construction
g["called_time_non_adaptor_sec"] = (
    g["called_time_total_sec"] - g["called_time_build_readset_sec"]
)

# Avoid divide-by-zero just in case
for c in [
    "called_time_readselection_sec",
    "called_time_solve_sec",
    "called_time_ps_sec",
    "called_time_write_sec",
]:
    g[f"{c}_frac_of_non_adaptor"] = g[c] / g["called_time_non_adaptor_sec"]

print("\nMean called-stage runtimes (including residual):\n")
print(
    g[
        [
            "called_time_build_readset_sec",
            "called_time_readselection_sec",
            "called_time_solve_sec",
            "called_time_ps_sec",
            "called_time_write_sec",
            "called_time_total_sec",
            "called_time_non_adaptor_sec",
        ]
    ].round(4)
)

print("\nFractions of non-adaptor runtime:\n")
print(
    g[
        [
            "called_time_readselection_sec_frac_of_non_adaptor",
            "called_time_solve_sec_frac_of_non_adaptor",
            "called_time_ps_sec_frac_of_non_adaptor",
            "called_time_write_sec_frac_of_non_adaptor",
        ]
    ].round(3)
)