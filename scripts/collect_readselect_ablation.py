import json, csv

rows = []
for mode, path in [
    ("With read selection (max_cov=15)", "output/prelim/readselect/wh_sel.summary.json"),
    ("No selection effect (max_cov=1e6)", "output/prelim/readselect/wh_nosel.summary.json"),
]:
    with open(path) as f:
        s = json.load(f)
    rows.append({
        "mode": mode,
        "max_coverage": s.get("max_coverage"),
        "selected_reads": s.get("selected_reads"),
        "num_informative_reads": s.get("num_informative_reads"),
        "time_total_sec": s.get("time_total_sec"),
        "time_readselection_sec": s.get("time_readselection_sec"),
        "time_solve_sec": s.get("time_solve_sec"),
        "time_build_readset_sec": s.get("time_build_readset_sec"),
    })

out = "output/prelim/readselect/readselect_ablation.csv"
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys())
    w.writeheader()
    w.writerows(rows)

print("Wrote:", out)
