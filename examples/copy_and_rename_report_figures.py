from pathlib import Path
import shutil

src_root = Path("output/experiments_draft")
dst_root = Path("report/figures/10_experiments")
dst_root.mkdir(parents=True, exist_ok=True)

mappings = [
    # 10.2 Baseline: depth scaling
    ("01_depth_curve/plots/call_recall.png", "fig_10_2_1_call_recall.png"),
    ("01_depth_curve/plots/oracle_effective_phased_recall.png", "fig_10_2_2_oracle_effective_phased_recall.png"),
    ("01_depth_curve/plots/called_effective_phased_recall.png", "fig_10_2_3_called_effective_phased_recall.png"),
    ("01_depth_curve/plots/oracle_num_phase_sets.png", "fig_10_2_4_oracle_num_phase_sets.png"),
    ("01_depth_curve/plots/called_num_phase_sets.png", "fig_10_2_5_called_num_phase_sets.png"),
    ("01_depth_curve/plots/called_switch_error.png", "fig_10_2_6_called_switch_error.png"),

    # 10.3 Duplicated regions
    ("03_ablation_duplication/plots/call_recall.png", "fig_10_3_1_call_recall.png"),
    ("03_ablation_duplication/plots/called_effective_phased_recall.png", "fig_10_3_2_called_effective_phased_recall.png"),
    ("03_ablation_duplication/plots/called_switch_error.png", "fig_10_3_3_called_switch_error.png"),
    ("03_ablation_duplication/plots/oracle_effective_phased_recall.png", "fig_10_3_4_oracle_effective_phased_recall.png"),
    ("03_ablation_duplication/plots/called_num_phase_sets.png", "fig_10_3_5_called_num_phase_sets.png"),

    # 10.4 Coverage dropout
    ("02_ablation_dropout/plots/oracle_num_phase_sets.png", "fig_10_4_1_oracle_num_phase_sets.png"),
    ("02_ablation_dropout/plots/called_num_phase_sets.png", "fig_10_4_2_called_num_phase_sets.png"),
    ("02_ablation_dropout/plots/oracle_effective_phased_recall.png", "fig_10_4_3_oracle_effective_phased_recall.png"),
    ("02_ablation_dropout/plots/called_effective_phased_recall.png", "fig_10_4_4_called_effective_phased_recall.png"),
    ("02_ablation_dropout/plots/call_recall.png", "fig_10_4_5_call_recall.png"),

    # 10.5 Correlated error bursts
    ("04_ablation_bursts/plots/called_switch_error.png", "fig_10_5_1_called_switch_error.png"),
    ("04_ablation_bursts/plots/called_effective_phased_recall.png", "fig_10_5_2_called_effective_phased_recall.png"),
    ("04_ablation_bursts/plots/call_recall.png", "fig_10_5_3_call_recall.png"),
    ("04_ablation_bursts/plots/called_num_phase_sets.png", "fig_10_5_4_called_num_phase_sets.png"),
    # 10.5.5 missing in archive

    # 10.6 Read length model
    ("06_ablation_lenmodel/plots/oracle_num_phase_sets.png", "fig_10_6_1_oracle_num_phase_sets.png"),
    ("06_ablation_lenmodel/plots/called_num_phase_sets.png", "fig_10_6_2_called_num_phase_sets.png"),
    ("06_ablation_lenmodel/plots/oracle_effective_phased_recall.png", "fig_10_6_3_oracle_effective_phased_recall.png"),
    ("06_ablation_lenmodel/plots/called_effective_phased_recall.png", "fig_10_6_4_called_effective_phased_recall.png"),
    ("06_ablation_lenmodel/plots/call_recall.png", "fig_10_6_5_call_recall.png"),

    # 10.7 Truth indels
    ("05_ablation_indels/plots/oracle_effective_phased_recall.png", "fig_10_7_1_oracle_effective_phased_recall.png"),
    ("05_ablation_indels/plots/called_effective_phased_recall.png", "fig_10_7_2_called_effective_phased_recall.png"),
    ("05_ablation_indels/plots/call_recall.png", "fig_10_7_3_call_recall.png"),
    ("05_ablation_indels/plots/called_num_phase_sets.png", "fig_10_7_4_called_num_phase_sets.png"),

    # 10.8 Interaction: duplication x dropout
    ("07_interaction_dup_x_dropout/plots/call_recall.png", "fig_10_8_1_call_recall.png"),
    ("07_interaction_dup_x_dropout/plots/called_effective_phased_recall.png", "fig_10_8_2_called_effective_phased_recall.png"),
    ("07_interaction_dup_x_dropout/plots/called_num_phase_sets.png", "fig_10_8_3_called_num_phase_sets.png"),
    ("07_interaction_dup_x_dropout/plots/called_switch_error.png", "fig_10_8_4_called_switch_error.png"),

    # 10.9.1 Sweep A: max_coverage
    ("08_optimize_whatshap/cov_sweep/plots/called_effective_phased_recall.png", "fig_10_9_1_1_called_effective_phased_recall.png"),
    ("08_optimize_whatshap/cov_sweep/plots/called_switch_error.png", "fig_10_9_1_2_called_switch_error.png"),
    ("08_optimize_whatshap/cov_sweep/plots/called_num_phase_sets.png", "fig_10_9_1_3_called_num_phase_sets.png"),
    ("08_optimize_whatshap/cov_sweep/plots/call_recall.png", "fig_10_9_1_4_call_recall.png"),

    # 10.9.2 Sweep B: phasing quality thresholds
    ("08_optimize_whatshap/qual_threshold_grid/plots/called_effective_phased_recall.png", "fig_10_9_2_1_called_effective_phased_recall.png"),
    ("08_optimize_whatshap/qual_threshold_grid/plots/called_switch_error.png", "fig_10_9_2_2_called_switch_error.png"),
    # 10.9.2.3 missing
    # 10.9.2.4 missing

    # 10.9.3 Sweep C: calling thresholds
    ("08_optimize_whatshap/calling_threshold_grid/plots/call_recall_heatmap.png", "fig_10_9_3_1_call_recall_heatmap.png"),
    ("08_optimize_whatshap/calling_threshold_grid/plots/called_effective_phased_recall_heatmap.png", "fig_10_9_3_2_called_effective_phased_recall_heatmap.png"),
    ("08_optimize_whatshap/calling_threshold_grid/plots/called_shared_het_recall_heatmap.png", "fig_10_9_3_3_called_shared_het_recall_heatmap.png"),
    ("08_optimize_whatshap/calling_threshold_grid/plots/call_precision_heatmap.png", "fig_10_9_3_4_call_precision_heatmap.png"),

    # 10.9.4 Sweep D: runtime-focused lower max_coverage
    ("08_optimize_whatshap/maxcov_runtime_sweep/plots/called_effective_phased_recall.png", "fig_10_9_4_1_called_effective_phased_recall.png"),
    ("08_optimize_whatshap/maxcov_runtime_sweep/plots/called_num_phase_sets.png", "fig_10_9_4_2_called_num_phase_sets.png"),
    # 10.9.4.3 missing
    ("08_optimize_whatshap/maxcov_runtime_sweep/plots/time_total_sec.png", "fig_10_9_4_4_time_total_sec.png"),

    # 10.9.5 Sweep E: fine phasing min_baseq
    ("08_optimize_whatshap/minbaseq_fine_sweep/plots/called_effective_phased_recall.png", "fig_10_9_5_1_called_effective_phased_recall.png"),
    ("08_optimize_whatshap/minbaseq_fine_sweep/plots/called_num_phase_sets.png", "fig_10_9_5_2_called_num_phase_sets.png"),
    ("08_optimize_whatshap/minbaseq_fine_sweep/plots/oracle_effective_phased_recall.png", "fig_10_9_5_3_oracle_effective_phased_recall.png"),
    ("08_optimize_whatshap/minbaseq_fine_sweep/plots/called_switch_error.png", "fig_10_9_5_4_called_switch_error.png"),

    # 10.9.6 Default vs optimized comparison
    ("08_optimize_whatshap/confirm_default_vs_optimized/plots/called_effective_phased_recall.png", "fig_10_9_6_1_called_effective_phased_recall.png"),
    ("08_optimize_whatshap/confirm_default_vs_optimized/plots/called_shared_het_recall.png", "fig_10_9_6_2_called_shared_het_recall.png"),
    # 10.9.6.3 missing
    ("08_optimize_whatshap/confirm_default_vs_optimized/plots/time_total_sec.png", "fig_10_9_6_4_time_total_sec.png"),

    # 10.9.7 Sweep F: local joint search
    ("08_optimize_whatshap/local_joint_bq_search/plots/called_effective_phased_recall_heatmap.png", "fig_10_9_7_1_called_effective_phased_recall_heatmap.png"),
    ("08_optimize_whatshap/local_joint_bq_search/plots/called_shared_het_recall_heatmap.png", "fig_10_9_7_2_called_shared_het_recall_heatmap.png"),
    ("08_optimize_whatshap/local_joint_bq_search/plots/call_recall_heatmap.png", "fig_10_9_7_3_call_recall_heatmap.png"),
    ("08_optimize_whatshap/local_joint_bq_search/plots/time_total_sec_heatmap.png", "fig_10_9_7_4_time_total_sec_heatmap.png"),

    # 10.9.8 Sweep G: frontier comparison
    ("08_optimize_whatshap/frontier_config_comparison/plots/called_effective_phased_recall.png", "fig_10_9_8_1_called_effective_phased_recall.png"),
    ("08_optimize_whatshap/frontier_config_comparison/plots/called_shared_het_recall.png", "fig_10_9_8_2_called_shared_het_recall.png"),
    ("08_optimize_whatshap/frontier_config_comparison/plots/called_num_phase_sets.png", "fig_10_9_8_3_called_num_phase_sets.png"),
    ("08_optimize_whatshap/frontier_config_comparison/plots/time_total_sec.png", "fig_10_9_8_4_time_total_sec.png"),

    # 10.9.9 Optimization robustness
    ("08_optimize_whatshap/transferability_across_scenarios/plots/called_effective_phased_recall_heatmap.png", "fig_10_9_9_1_called_effective_phased_recall_heatmap.png"),
    ("08_optimize_whatshap/transferability_across_scenarios/plots/called_shared_het_recall_heatmap.png", "fig_10_9_9_2_called_shared_het_recall_heatmap.png"),
    ("08_optimize_whatshap/transferability_across_scenarios/plots/called_num_phase_sets_heatmap.png", "fig_10_9_9_3_called_num_phase_sets_heatmap.png"),
    ("08_optimize_whatshap/transferability_across_scenarios/plots/time_total_sec_heatmap.png", "fig_10_9_9_4_time_total_sec_heatmap.png"),

    # 10.9.10 Sweep H: call_min_mapq rule-out
    ("08_optimize_whatshap/callmapq_ruleout_sweep/plots/called_effective_phased_recall.png", "fig_10_9_10_1_called_effective_phased_recall.png"),
    ("08_optimize_whatshap/callmapq_ruleout_sweep/plots/call_recall.png", "fig_10_9_10_2_call_recall.png"),
    ("08_optimize_whatshap/callmapq_ruleout_sweep/plots/called_shared_het_recall.png", "fig_10_9_10_3_called_shared_het_recall.png"),
    ("08_optimize_whatshap/callmapq_ruleout_sweep/plots/call_precision.png", "fig_10_9_10_4_call_precision.png"),

    # Optional custom combined baseline figure, if you want to keep it for the report or poster
    ("01_depth_curve/plots/poster_baseline_oracle_vs_called.png", "fig_10_2_combined_oracle_vs_called.png"),
]

missing = []
copied = []

for rel_src, new_name in mappings:
    src = src_root / rel_src
    dst = dst_root / new_name
    if src.exists():
        shutil.copy2(src, dst)
        copied.append((src, dst))
    else:
        missing.append(rel_src)

print(f"Copied {len(copied)} figures to {dst_root}")
if missing:
    print("\nMissing source plots:")
    for m in missing:
        print(" -", m)