# benchmark/experiment_driver.py
"""
Run end-to-end experiments for the long-read (ONT) WhatsHap benchmarking pipeline.

Usage:
  python -m benchmark.experiment_driver --outdir output/experiments --seeds 0 1 2

You can run only some sections:
  python -m benchmark.experiment_driver --outdir output/experiments --seeds 0 1 2 --only depth,dropout,dup
  
By default skips a run if prefix.pipeline.json exists, use --force to re-run anyway:
  python -m benchmark.experiment_driver --outdir output/experiments_all --seeds 0 1 2 --force

Reality check:
  python -m benchmark.experiment_driver \
  --outdir output/experiments_all \
  --seeds 0 1 2 \
  --only reality \
  --real-fastq /path/to/real.fastq.gz \
  --real-bam /path/to/real.bam

This driver calls:
  - python -m benchmark.longread_pipeline_runner
  - python -m benchmark.aggregate_pipeline_reports_full
and generates plots (matplotlib) under each experiment directory.

Notes:
- If indels are enabled (num_indels > 0), we automatically pass --phase-snps-only to keep SNP phasing/eval valid.
- This driver is intentionally simple and sequential (easy to debug / reproduce). Use --force to re-run existing prefixes.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import pandas as pd
import matplotlib.pyplot as plt


def _run(cmd: List[str]) -> None:
    print("[RUN]", " ".join(cmd))
    subprocess.check_call(cmd)


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _prefix_done(prefix: Path) -> bool:
    return prefix.with_suffix(".pipeline.json").exists()


def _call_pipeline(prefix: Path, args: Dict[str, Any], force: bool) -> None:
    """
    Call benchmark.longread_pipeline_runner with args dict.

    args keys should correspond to CLI flags without leading '--'.
    Example: {"seed": 0, "ref-length": 80000, "num-snps": 800, ...}
    """
    if (not force) and _prefix_done(prefix):
        print(f"[SKIP] exists: {prefix}.pipeline.json")
        return

    cmd = [sys.executable, "-m", "benchmark.longread_pipeline_runner", "--prefix", str(prefix)]
    for k, v in args.items():
        if v is None:
            continue
        flag = "--" + k.replace("_", "-")
        if isinstance(v, bool):
            if v:
                cmd.append(flag)
        else:
            cmd.extend([flag, str(v)])
    _run(cmd)


def _aggregate(exp_dir: Path) -> Path:
    out_csv = exp_dir / "aggregate.csv"
    _run([sys.executable, "-m", "benchmark.aggregate_pipeline_reports_full",
          "--root", str(exp_dir), "--out", str(out_csv)])
    return out_csv


def _plot_errorbar(df: pd.DataFrame, xcol: str, metric: str, out_png: Path, ylabel: Optional[str] = None) -> None:
    xs = sorted(df[xcol].unique())
    g = df.groupby(xcol)[metric].agg(["mean", "std"]).reindex(xs)
    m = g["mean"].values
    s = g["std"].fillna(0.0).values

    plt.figure()
    plt.errorbar(xs, m, yerr=s, marker="o", capsize=3)
    plt.xlabel(xcol)
    plt.ylabel(ylabel or metric)
    plt.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close()


def _plots_basic(exp_dir: Path, df: pd.DataFrame, xcol: str, metrics: List[Tuple[str, str]]) -> None:
    pdir = exp_dir / "plots"
    _ensure_dir(pdir)
    for metric, ylabel in metrics:
        if metric not in df.columns:
            continue
        _plot_errorbar(df, xcol, metric, pdir / f"{metric}.png", ylabel=ylabel)


def _plot_heatmap(df: pd.DataFrame, xcol: str, ycol: str, metric: str, out_png: Path, title: Optional[str] = None) -> None:
    if metric not in df.columns or xcol not in df.columns or ycol not in df.columns:
        return
    xs = sorted(df[xcol].unique())
    ys = sorted(df[ycol].unique())
    piv = df.groupby([ycol, xcol])[metric].mean().unstack(xcol).reindex(index=ys, columns=xs)

    plt.figure(figsize=(5.5, 4.2))
    im = plt.imshow(piv.values, aspect="auto", origin="lower")
    plt.xticks(range(len(xs)), xs)
    plt.yticks(range(len(ys)), ys)
    plt.xlabel(xcol)
    plt.ylabel(ycol)
    plt.title(title or metric)
    plt.colorbar(im, label=metric)

    for yi, y in enumerate(ys):
        for xi, x in enumerate(xs):
            val = piv.loc[y, x]
            if pd.notna(val):
                plt.text(xi, yi, f"{val:.3f}", ha="center", va="center", fontsize=8)

    plt.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close()


def _plots_heatmap(exp_dir: Path, df: pd.DataFrame, xcol: str, ycol: str, metrics: List[Tuple[str, str]]) -> None:
    pdir = exp_dir / "plots"
    _ensure_dir(pdir)
    for metric, title in metrics:
        _plot_heatmap(df, xcol, ycol, metric, pdir / f"{metric}_heatmap.png", title=title)


def _hard_optimize_args(base: Dict[str, Any]) -> Dict[str, Any]:
    hard = dict(base)
    hard.update({
        "ref_length": 120000,
        "num_snps": 1200,
        "num_reads": 300,
        "dup_segments": 5,
        "dup_len": 3000,
        "dup_min_gap": 500,
        "start_model": "dropout",
        "dropout_fraction": 0.1,
        "dropout_block_len": 1000,
        "burst_prob": 0.6,
        "burst_count": 2,
        "burst_len": 300,
        "burst_mult": 8.0,
        "num_indels": 120,
        "indel_min_len": 1,
        "indel_max_len": 5,
        "indel_het_rate": 0.5,
    })
    _auto_indel_phase_flags(hard)
    return hard


def _bottleneck_decomposition(df: pd.DataFrame, tag: str, out_csv: Path) -> pd.DataFrame:
    """
    Compute a simple decomposition for CALLED and ORACLE:
      - shared_het_recall: how many truth het SNPs survive calling/intersection
      - phasing_rate_shared_het: of shared hets, what fraction are phased
      - phase_accuracy: of phased hets, what fraction are correct (best-flip)
    We map this into:
      - effective_phased_recall = shared_het_recall * phasing_rate_shared_het * phase_accuracy
    plus implied losses.

    tag: "called" or "oracle"
    """
    sh = f"{tag}_shared_het_recall"
    pr = f"{tag}_phasing_rate_shared_het"
    acc = f"{tag}_phase_accuracy"
    eff = f"{tag}_effective_phased_recall"

    use = df.copy()
    for c in (sh, pr, acc, eff):
        if c not in use.columns:
            use[c] = float("nan")

    # derive "loss" terms (not a perfect partition, but good for narrative)
    use[f"{tag}_loss_calling_or_sites"] = 1.0 - use[sh]
    use[f"{tag}_loss_unphased"] = use[sh] * (1.0 - use[pr])
    use[f"{tag}_loss_switch_or_phase_error"] = use[sh] * use[pr] * (1.0 - use[acc])
    use[f"{tag}_effective_phased_recall_check"] = use[sh] * use[pr] * use[acc]

    cols = [
        "pipeline_json", "seed",
        sh, pr, acc, eff,
        f"{tag}_loss_calling_or_sites",
        f"{tag}_loss_unphased",
        f"{tag}_loss_switch_or_phase_error",
        f"{tag}_effective_phased_recall_check",
    ]
    cols = [c for c in cols if c in use.columns]
    out = use[cols]
    out.to_csv(out_csv, index=False)
    return out


def _default_base_args() -> Dict[str, Any]:
    return {
        "seed": 0,
        "ref_length": 80000,
        "num_snps": 800,
        "het_rate": 0.8,
        "avoid_regions": False,
        "num_reads": 200,
        "min_len": 2000,
        "max_len": 6000,
        "platform": "ont",
        "ont_profile": "q20",
        "map_preset": "map-ont",
        "call_min_mapq": 20,
        "call_min_baseq": 15,
        "vcf_source": "both",
        "max_coverage": 15,
        "min_mapq": 20,
        "min_baseq": 20,
        "dup_segments": 0,
        "dup_len": 3000,
        "dup_min_gap": 500,
        "len_model": "uniform",
        "ln_mean": 8.3,
        "ln_sigma": 0.6,
        "start_model": "uniform",
        "dropout_fraction": 0.0,
        "dropout_block_len": 1000,
        "burst_prob": 0.0,
        "burst_count": 1,
        "burst_len": 200,
        "burst_mult": 5.0,
        "num_indels": 0,
        "indel_min_len": 1,
        "indel_max_len": 5,
        "indel_het_rate": 0.5,
        # important: keep SNP-only phasing/eval valid when indels exist
        "phase_snps_only": False,
        "eval_snps_only": False,
        # reference preset
        "ref_preset": "plain",
    }


def _auto_indel_phase_flags(args: Dict[str, Any]) -> None:
    if int(args.get("num_indels", 0) or 0) > 0:
        # Keep SNP phasing stable. This is the validated approach.
        args["phase_snps_only"] = True
        # eval-snps-only can be implicit in your runner; safe to set anyway.
        args["eval_snps_only"] = True


def run_depth(outdir: Path, seeds: List[int], force: bool, base: Dict[str, Any]) -> None:
    exp = outdir / "01_depth_curve"
    _ensure_dir(exp)
    reads_list = [50, 100, 200, 400]

    for r in reads_list:
        for s in seeds:
            args = dict(base)
            args["seed"] = s
            args["num_reads"] = r
            args["dup_segments"] = 0
            args["start_model"] = "uniform"
            args["dropout_fraction"] = 0.0
            args["burst_prob"] = 0.0
            args["num_indels"] = 0
            args["len_model"] = "uniform"
            args["phase_snps_only"] = False
            args["eval_snps_only"] = False
            prefix = exp / f"r{r}_s{s}"
            _call_pipeline(prefix, args, force=force)

    csvp = _aggregate(exp)
    df = pd.read_csv(csvp)
    _plots_basic(exp, df, "num_reads", [
        ("call_recall", "Calling recall"),
        ("oracle_effective_phased_recall", "Oracle effective phased recall"),
        ("called_effective_phased_recall", "Called effective phased recall"),
        ("oracle_num_phase_sets", "Oracle num phase sets"),
        ("called_num_phase_sets", "Called num phase sets"),
        ("called_switch_error", "Called switch error rate"),
    ])
    _bottleneck_decomposition(df, "called", exp / "decomposition_called.csv")
    _bottleneck_decomposition(df, "oracle", exp / "decomposition_oracle.csv")


def run_dropout(outdir: Path, seeds: List[int], force: bool, base: Dict[str, Any]) -> None:
    exp = outdir / "02_ablation_dropout"
    _ensure_dir(exp)
    fracs = [0.0, 0.05, 0.10, 0.20]

    for f in fracs:
        for s in seeds:
            args = dict(base)
            args["seed"] = s
            args["start_model"] = "dropout"
            args["dropout_fraction"] = f
            args["dropout_block_len"] = 1000
            args["dup_segments"] = 0
            args["burst_prob"] = 0.0
            args["num_indels"] = 0
            args["phase_snps_only"] = False
            args["eval_snps_only"] = False
            prefix = exp / f"df{f}_s{s}"
            _call_pipeline(prefix, args, force=force)

    csvp = _aggregate(exp)
    df = pd.read_csv(csvp)
    _plots_basic(exp, df, "dropout_fraction", [
        ("call_recall", "Calling recall"),
        ("oracle_num_phase_sets", "Oracle num phase sets"),
        ("called_num_phase_sets", "Called num phase sets"),
        ("oracle_effective_phased_recall", "Oracle effective phased recall"),
        ("called_effective_phased_recall", "Called effective phased recall"),
    ])
    _bottleneck_decomposition(df, "called", exp / "decomposition_called.csv")


def run_dup(outdir: Path, seeds: List[int], force: bool, base: Dict[str, Any]) -> None:
    exp = outdir / "03_ablation_duplication"
    _ensure_dir(exp)
    dups = [0, 1, 3, 5]

    for d in dups:
        for s in seeds:
            args = dict(base)
            args["seed"] = s
            args["dup_segments"] = d
            args["start_model"] = "uniform"
            args["dropout_fraction"] = 0.0
            args["burst_prob"] = 0.0
            args["num_indels"] = 0
            args["phase_snps_only"] = False
            args["eval_snps_only"] = False
            prefix = exp / f"dup{d}_s{s}"
            _call_pipeline(prefix, args, force=force)

    csvp = _aggregate(exp)
    df = pd.read_csv(csvp)
    _plots_basic(exp, df, "dup_segments", [
        ("call_recall", "Calling recall"),
        ("called_switch_error", "Called switch error rate"),
        ("called_num_phase_sets", "Called num phase sets"),
        ("called_effective_phased_recall", "Called effective phased recall"),
        ("oracle_effective_phased_recall", "Oracle effective phased recall"),
    ])
    _bottleneck_decomposition(df, "called", exp / "decomposition_called.csv")


def run_bursts(outdir: Path, seeds: List[int], force: bool, base: Dict[str, Any]) -> None:
    exp = outdir / "04_ablation_bursts"
    _ensure_dir(exp)
    probs = [0.0, 0.3, 0.6]

    for p in probs:
        for s in seeds:
            args = dict(base)
            args["seed"] = s
            args["burst_prob"] = p
            args["burst_count"] = 2
            args["burst_len"] = 300
            args["burst_mult"] = 8.0
            args["dup_segments"] = 0
            args["start_model"] = "uniform"
            args["dropout_fraction"] = 0.0
            args["num_indels"] = 0
            args["phase_snps_only"] = False
            args["eval_snps_only"] = False
            prefix = exp / f"bp{p}_s{s}"
            _call_pipeline(prefix, args, force=force)

    csvp = _aggregate(exp)
    df = pd.read_csv(csvp)
    _plots_basic(exp, df, "burst_prob", [
        ("call_recall", "Calling recall"),
        ("called_effective_phased_recall", "Called effective phased recall"),
        ("called_switch_error", "Called switch error rate"),
        ("called_num_phase_sets", "Called num phase sets"),
    ])
    _bottleneck_decomposition(df, "called", exp / "decomposition_called.csv")


def run_indels(outdir: Path, seeds: List[int], force: bool, base: Dict[str, Any]) -> None:
    exp = outdir / "05_ablation_indels"
    _ensure_dir(exp)
    indels = [0, 80, 200]

    for k in indels:
        for s in seeds:
            args = dict(base)
            args["seed"] = s
            args["num_indels"] = k
            args["indel_min_len"] = 1
            args["indel_max_len"] = 5
            args["indel_het_rate"] = 0.5
            args["dup_segments"] = 0
            args["start_model"] = "uniform"
            args["dropout_fraction"] = 0.0
            args["burst_prob"] = 0.0
            _auto_indel_phase_flags(args)
            prefix = exp / f"indel{k}_s{s}"
            _call_pipeline(prefix, args, force=force)

    csvp = _aggregate(exp)
    df = pd.read_csv(csvp)
    # If your CSV doesn't include num_indels yet, x-axis won't work; we'll fall back to parsing from filename
    if "num_indels" not in df.columns:
        # derive from pipeline_json path
        df["num_indels"] = df["pipeline_json"].str.extract(r"indel(\d+)_s").astype(float)
    _plots_basic(exp, df, "num_indels", [
        ("call_recall", "Calling recall"),
        ("called_effective_phased_recall", "Called effective phased recall"),
        ("oracle_effective_phased_recall", "Oracle effective phased recall"),
        ("called_num_phase_sets", "Called num phase sets"),
    ])


def run_lenmodel(outdir: Path, seeds: List[int], force: bool, base: Dict[str, Any]) -> None:
    exp = outdir / "06_ablation_lenmodel"
    _ensure_dir(exp)
    models = ["uniform", "lognormal"]

    for m in models:
        for s in seeds:
            args = dict(base)
            args["seed"] = s
            args["len_model"] = m
            if m == "lognormal":
                args["ln_mean"] = 8.3
                args["ln_sigma"] = 0.6
            args["dup_segments"] = 0
            args["start_model"] = "uniform"
            args["burst_prob"] = 0.0
            args["num_indels"] = 0
            args["phase_snps_only"] = False
            args["eval_snps_only"] = False
            prefix = exp / f"{m}_s{s}"
            _call_pipeline(prefix, args, force=force)

    csvp = _aggregate(exp)
    df = pd.read_csv(csvp)
    # make x axis categorical model encoded as 0/1 for plotting
    df["_lenmodel_x"] = df["len_model"].map({"uniform": 0, "lognormal": 1})
    _plots_basic(exp, df, "_lenmodel_x", [
        ("call_recall", "Calling recall"),
        ("oracle_effective_phased_recall", "Oracle effective phased recall"),
        ("called_effective_phased_recall", "Called effective phased recall"),
        ("oracle_num_phase_sets", "Oracle num phase sets"),
        ("called_num_phase_sets", "Called num phase sets"),
    ])


def run_interaction(outdir: Path, seeds: List[int], force: bool, base: Dict[str, Any]) -> None:
    exp = outdir / "07_interaction_dup_x_dropout"
    _ensure_dir(exp)
    dups = [0, 5]
    dfs = [0.0, 0.1]

    for d in dups:
        for f in dfs:
            for s in seeds:
                args = dict(base)
                args["seed"] = s
                args["dup_segments"] = d
                args["start_model"] = "dropout" if f > 0 else "uniform"
                args["dropout_fraction"] = f
                args["dropout_block_len"] = 1000
                args["burst_prob"] = 0.0
                args["num_indels"] = 0
                args["phase_snps_only"] = False
                args["eval_snps_only"] = False
                prefix = exp / f"dup{d}_df{f}_s{s}"
                _call_pipeline(prefix, args, force=force)

    csvp = _aggregate(exp)
    df = pd.read_csv(csvp)
    # plot by a combined categorical label
    df["cond"] = df.apply(lambda r: f"dup{int(r['dup_segments'])}_df{float(r['dropout_fraction']):.2f}", axis=1)
    # order
    order = [f"dup{d}_df{f:.2f}" for d in dups for f in dfs]
    df["cond_i"] = df["cond"].apply(lambda c: order.index(c))
    _plots_basic(exp, df, "cond_i", [
        ("called_effective_phased_recall", "Called effective phased recall"),
        ("called_num_phase_sets", "Called num phase sets"),
        ("called_switch_error", "Called switch error rate"),
        ("call_recall", "Calling recall"),
    ])


def run_optimize(outdir: Path, seeds: List[int], force: bool, base: Dict[str, Any], sweeps: Optional[List[str]] = None) -> None:
    """
    WhatsHap optimization sweeps on a fixed "hard" scenario.

    Available sweeps:
      - cov: existing max_coverage sweep
      - qual: existing phasing min_mapq x min_baseq grid
      - calling: caller-side call_min_mapq x call_min_baseq grid
      - runtime: lower max_coverage sweep to hunt for runtime savings
      - bqfine: fine-grained phasing min_baseq sweep
      - local: local joint search around the current best caller/phasing base-quality settings
      - frontier: representative accuracy/runtime frontier comparison
      - confirm: default-vs-optimized confirmation comparison
    """
    exp = outdir / "08_optimize_whatshap"
    _ensure_dir(exp)

    hard = _hard_optimize_args(base)
    selected = set(sweeps or ["cov", "qual", "calling", "runtime", "bqfine"])

    # 8A) max_coverage sweep
    if "cov" in selected:
        cov_dir = exp / "cov_sweep"
        _ensure_dir(cov_dir)
        for mc in [10, 15, 25, 40]:
            for s in seeds:
                args = dict(hard)
                args["seed"] = s
                args["max_coverage"] = mc
                prefix = cov_dir / f"mc{mc}_s{s}"
                _call_pipeline(prefix, args, force=force)
        csvp = _aggregate(cov_dir)
        df = pd.read_csv(csvp)
        _plots_basic(cov_dir, df, "max_coverage", [
            ("called_effective_phased_recall", "Called effective phased recall"),
            ("called_switch_error", "Called switch error rate"),
            ("called_num_phase_sets", "Called num phase sets"),
            ("call_recall", "Calling recall"),
            ("time_total_sec", "Total runtime (s)"),
        ])
        _bottleneck_decomposition(df, "called", cov_dir / "decomposition_called.csv")

    # 8B) phasing min_mapq x min_baseq grid
    if "qual" in selected:
        qb_dir = exp / "qual_threshold_grid"
        _ensure_dir(qb_dir)
        for mq in [0, 10, 20]:
            for bq in [0, 10, 20]:
                for s in seeds:
                    args = dict(hard)
                    args["seed"] = s
                    args["min_mapq"] = mq
                    args["min_baseq"] = bq
                    prefix = qb_dir / f"mq{mq}_bq{bq}_s{s}"
                    _call_pipeline(prefix, args, force=force)
        csvp2 = _aggregate(qb_dir)
        df2 = pd.read_csv(csvp2)
        _plots_heatmap(qb_dir, df2, "min_mapq", "min_baseq", [
            ("called_effective_phased_recall", "Called effective phased recall"),
            ("called_switch_error", "Called switch error"),
            ("called_num_phase_sets", "Called num phase sets"),
            ("oracle_effective_phased_recall", "Oracle effective phased recall"),
        ])
        _bottleneck_decomposition(df2, "called", qb_dir / "decomposition_called.csv")

    # 8C) NEW: caller-side call_min_mapq x call_min_baseq grid
    if "calling" in selected:
        call_dir = exp / "calling_threshold_grid"
        _ensure_dir(call_dir)
        for cmq in [0, 10, 20]:
            for cbq in [5, 10, 15, 20]:
                for s in seeds:
                    args = dict(hard)
                    args["seed"] = s
                    args["call_min_mapq"] = cmq
                    args["call_min_baseq"] = cbq
                    # Use the best phasing-side settings found so far.
                    args["max_coverage"] = 10
                    args["min_mapq"] = 20
                    args["min_baseq"] = 10
                    prefix = call_dir / f"cmq{cmq}_cbq{cbq}_s{s}"
                    _call_pipeline(prefix, args, force=force)
        csvp3 = _aggregate(call_dir)
        df3 = pd.read_csv(csvp3)
        _plots_heatmap(call_dir, df3, "call_min_mapq", "call_min_baseq", [
            ("call_recall", "Calling recall"),
            ("call_precision", "Calling precision"),
            ("called_shared_het_recall", "Called shared het recall"),
            ("called_effective_phased_recall", "Called effective phased recall"),
        ])
        _bottleneck_decomposition(df3, "called", call_dir / "decomposition_called.csv")

    # 8D) NEW: runtime-oriented lower max_coverage sweep
    if "runtime" in selected:
        rt_dir = exp / "maxcov_runtime_sweep"
        _ensure_dir(rt_dir)
        for mc in [4, 6, 8, 10, 15]:
            for s in seeds:
                args = dict(hard)
                args["seed"] = s
                args["call_min_mapq"] = 20
                args["call_min_baseq"] = 15
                args["min_mapq"] = 20
                args["min_baseq"] = 10
                args["max_coverage"] = mc
                prefix = rt_dir / f"mc{mc}_s{s}"
                _call_pipeline(prefix, args, force=force)
        csvp4 = _aggregate(rt_dir)
        df4 = pd.read_csv(csvp4)
        _plots_basic(rt_dir, df4, "max_coverage", [
            ("called_effective_phased_recall", "Called effective phased recall"),
            ("called_num_phase_sets", "Called num phase sets"),
            ("called_switch_error", "Called switch error rate"),
            ("time_total_sec", "Total runtime (s)"),
        ])
        _bottleneck_decomposition(df4, "called", rt_dir / "decomposition_called.csv")

    # 8E) NEW: fine-grained phasing min_baseq sweep
    if "bqfine" in selected:
        bq_dir = exp / "minbaseq_fine_sweep"
        _ensure_dir(bq_dir)
        for bq in [0, 5, 10, 15, 20]:
            for s in seeds:
                args = dict(hard)
                args["seed"] = s
                args["call_min_mapq"] = 20
                args["call_min_baseq"] = 15
                args["max_coverage"] = 10
                args["min_mapq"] = 20
                args["min_baseq"] = bq
                prefix = bq_dir / f"bq{bq}_s{s}"
                _call_pipeline(prefix, args, force=force)
        csvp5 = _aggregate(bq_dir)
        df5 = pd.read_csv(csvp5)
        _plots_basic(bq_dir, df5, "min_baseq", [
            ("oracle_effective_phased_recall", "Oracle effective phased recall"),
            ("called_effective_phased_recall", "Called effective phased recall"),
            ("called_num_phase_sets", "Called num phase sets"),
            ("called_switch_error", "Called switch error rate"),
        ])
        _bottleneck_decomposition(df5, "called", bq_dir / "decomposition_called.csv")


    # 8F) NEW: local joint caller/phasing base-quality search
    if "local" in selected:
        local_dir = exp / "local_joint_bq_search"
        _ensure_dir(local_dir)
        for cbq in [5, 10, 15]:
            for pbq in [5, 10, 15]:
                for s in seeds:
                    args = dict(hard)
                    args["seed"] = s
                    args["call_min_mapq"] = 20
                    args["call_min_baseq"] = cbq
                    args["max_coverage"] = 8
                    args["min_mapq"] = 20
                    args["min_baseq"] = pbq
                    prefix = local_dir / f"cbq{cbq}_pbq{pbq}_s{s}"
                    _call_pipeline(prefix, args, force=force)
        csvp_local = _aggregate(local_dir)
        df_local = pd.read_csv(csvp_local)
        _plots_heatmap(local_dir, df_local, "call_min_baseq", "min_baseq", [
            ("call_recall", "Calling recall"),
            ("called_shared_het_recall", "Called shared het recall"),
            ("called_effective_phased_recall", "Called effective phased recall"),
            ("time_total_sec", "Total runtime (s)"),
        ])
        _bottleneck_decomposition(df_local, "called", local_dir / "decomposition_called.csv")


    # 8G) NEW: representative accuracy/runtime frontier comparison
    if "frontier" in selected:
        fr_dir = exp / "frontier_config_comparison"
        _ensure_dir(fr_dir)
        configs = [
            ("default", {
                "call_min_mapq": 20,
                "call_min_baseq": 15,
                "max_coverage": 15,
                "min_mapq": 20,
                "min_baseq": 20,
            }),
            ("caller_only", {
                "call_min_mapq": 20,
                "call_min_baseq": 10,
                "max_coverage": 15,
                "min_mapq": 20,
                "min_baseq": 20,
            }),
            ("phasing_only", {
                "call_min_mapq": 20,
                "call_min_baseq": 15,
                "max_coverage": 8,
                "min_mapq": 20,
                "min_baseq": 10,
            }),
            ("balanced", {
                "call_min_mapq": 20,
                "call_min_baseq": 10,
                "max_coverage": 8,
                "min_mapq": 20,
                "min_baseq": 10,
            }),
            ("runtime", {
                "call_min_mapq": 20,
                "call_min_baseq": 10,
                "max_coverage": 6,
                "min_mapq": 20,
                "min_baseq": 10,
            }),
        ]
        for label, overrides in configs:
            for s in seeds:
                args = dict(hard)
                args["seed"] = s
                args.update(overrides)
                prefix = fr_dir / f"{label}_s{s}"
                _call_pipeline(prefix, args, force=force)
        csvp_fr = _aggregate(fr_dir)
        df_fr = pd.read_csv(csvp_fr)
        if "pipeline_json" in df_fr.columns:
            df_fr["config_label"] = df_fr["pipeline_json"].astype(str).str.extract(r"/(default|caller_only|phasing_only|balanced|runtime)_s\d+\.pipeline\.json$")[0]
        elif "prefix" in df_fr.columns:
            df_fr["config_label"] = df_fr["prefix"].astype(str).str.extract(r"/(default|caller_only|phasing_only|balanced|runtime)_s\d+$")[0]
        else:
            df_fr["config_label"] = None
        if df_fr["config_label"].isna().all():
            df_fr["config_label"] = [
                "caller_only" if "caller_only" in str(x) else
                "phasing_only" if "phasing_only" in str(x) else
                "balanced" if "balanced" in str(x) else
                "runtime" if "runtime" in str(x) else
                "default"
                for x in df_fr.get("pipeline_json", df_fr.index)
            ]
        df_fr.to_csv(csvp_fr, index=False)
        _plots_basic(fr_dir, df_fr, "config_label", [
            ("called_effective_phased_recall", "Called effective phased recall"),
            ("called_shared_het_recall", "Called shared het recall"),
            ("called_num_phase_sets", "Called num phase sets"),
            ("time_total_sec", "Total runtime (s)"),
        ])
        _bottleneck_decomposition(df_fr, "called", fr_dir / "decomposition_called.csv")


    # 8H) NEW: default-vs-optimized confirmation comparison
    if "confirm" in selected:
        cf_dir = exp / "confirm_default_vs_optimized"
        _ensure_dir(cf_dir)
        configs = [
            ("default", {
                "call_min_mapq": 20,
                "call_min_baseq": 15,
                "max_coverage": 15,
                "min_mapq": 20,
                "min_baseq": 20,
            }),
            ("optimized", {
                "call_min_mapq": 20,
                "call_min_baseq": 10,
                "max_coverage": 8,
                "min_mapq": 20,
                "min_baseq": 10,
            }),
        ]
        for label, overrides in configs:
            for s in seeds:
                args = dict(hard)
                args["seed"] = s
                args.update(overrides)
                prefix = cf_dir / f"{label}_s{s}"
                _call_pipeline(prefix, args, force=force)
        csvp6 = _aggregate(cf_dir)
        df6 = pd.read_csv(csvp6)
        # Derive a stable categorical column from the output path/prefix naming.
        if "pipeline_json" in df6.columns:
            df6["config_label"] = df6["pipeline_json"].astype(str).str.extract(r"/(default|optimized)_s\d+\.pipeline\.json$")[0]
        elif "prefix" in df6.columns:
            df6["config_label"] = df6["prefix"].astype(str).str.extract(r"/(default|optimized)_s\d+$")[0]
        else:
            df6["config_label"] = None
        if df6["config_label"].isna().all():
            # Fallback in case paths are formatted differently.
            df6["config_label"] = ["optimized" if "optimized" in str(x) else "default" for x in df6.get("pipeline_json", df6.index)]
        df6.to_csv(csvp6, index=False)
        _plots_basic(cf_dir, df6, "config_label", [
            ("called_effective_phased_recall", "Called effective phased recall"),
            ("called_shared_het_recall", "Called shared het recall"),
            ("call_recall", "Calling recall"),
            ("time_total_sec", "Total runtime (s)"),
        ])
        _bottleneck_decomposition(df6, "called", cf_dir / "decomposition_called.csv")


def run_reality_check(outdir: Path, real_fastq: Optional[str], real_bam: Optional[str]) -> None:
    """
    Optional: compare simulated vs a real dataset if you have one.
    If not provided, this section is skipped.
    """
    exp = outdir / "09_reality_check"
    _ensure_dir(exp)

    if not real_fastq and not real_bam:
        (exp / "README.txt").write_text(
            "Reality check skipped (no --real-fastq/--real-bam provided).\n"
            "If you have a real ONT FASTQ/BAM, rerun:\n"
            "  python -m benchmark.experiment_driver --outdir ... --real-fastq path --real-bam path\n"
        )
        print("[INFO] reality check skipped (no real data paths)")
        return

    def read_fastq_lengths(path: str, limit: int = 200000) -> List[int]:
        lens: List[int] = []
        import gzip
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as f:
            while True:
                h = f.readline()
                if not h:
                    break
                seq = f.readline().strip()
                f.readline()
                qual = f.readline().strip()
                lens.append(len(seq))
                if len(lens) >= limit:
                    break
        return lens

    if real_fastq:
        L = read_fastq_lengths(real_fastq, limit=200000)
        import numpy as np
        arr = np.array(L)
        stats = {
            "n": int(arr.size),
            "min": int(arr.min()),
            "mean": float(arr.mean()),
            "median": float(np.median(arr)),
            "max": int(arr.max()),
        }
        (exp / "real_fastq_length_stats.json").write_text(json.dumps(stats, indent=2))
        # histogram
        plt.figure()
        plt.hist(arr, bins=100)
        plt.xlabel("read length")
        plt.ylabel("count")
        plt.savefig(exp / "real_fastq_read_length_hist.png", dpi=200, bbox_inches="tight")
        plt.close()
        print("[INFO] wrote reality check FASTQ stats/plot")

    if real_bam:
        # basic mapping stats via samtools flagstat
        out = subprocess.check_output(["samtools", "flagstat", real_bam], text=True)
        (exp / "real_bam_flagstat.txt").write_text(out)
        print("[INFO] wrote reality check BAM flagstat")


def main():
    ap = argparse.ArgumentParser(description="Run all long-read benchmarking experiment directions.")
    ap.add_argument("--outdir", required=True, help="Root output directory, e.g. output/experiments_final")
    ap.add_argument("--seeds", nargs="+", type=int, default=[0, 1, 2], help="List of seeds to run")
    ap.add_argument("--force", action="store_true", help="Re-run even if prefix.pipeline.json exists")
    ap.add_argument(
        "--only",
        default="",
        help="Comma-separated sections to run. Options: depth,dropout,dup,bursts,indels,lenmodel,interaction,optimize,reality. "
             "If empty, runs all."
    )
    ap.add_argument(
        "--optimize-sweeps",
        default="cov,qual,calling,runtime,bqfine",
        help="Comma-separated optimize sub-sweeps: cov,qual,calling,runtime,bqfine,local,confirm. Used when optimize is selected."
    )
    ap.add_argument("--real-fastq", default=None, help="Optional real FASTQ(.gz) for reality check")
    ap.add_argument("--real-bam", default=None, help="Optional real BAM for reality check")

    args = ap.parse_args()
    outdir = Path(args.outdir)
    _ensure_dir(outdir)

    only = [x.strip() for x in args.only.split(",") if x.strip()]
    run_all = not only

    base = _default_base_args()
    seeds = args.seeds

    if run_all or "depth" in only:
        run_depth(outdir, seeds, args.force, base)
    if run_all or "dropout" in only:
        run_dropout(outdir, seeds, args.force, base)
    if run_all or "dup" in only:
        run_dup(outdir, seeds, args.force, base)
    if run_all or "bursts" in only:
        run_bursts(outdir, seeds, args.force, base)
    if run_all or "indels" in only:
        run_indels(outdir, seeds, args.force, base)
    if run_all or "lenmodel" in only:
        run_lenmodel(outdir, seeds, args.force, base)
    if run_all or "interaction" in only:
        run_interaction(outdir, seeds, args.force, base)
    optimize_sweeps = [x.strip() for x in args.optimize_sweeps.split(",") if x.strip()]

    if run_all or "optimize" in only:
        run_optimize(outdir, seeds, args.force, base, sweeps=optimize_sweeps)
    if run_all or "reality" in only:
        run_reality_check(outdir, args.real_fastq, args.real_bam)

    (outdir / "DONE.txt").write_text("All requested experiment sections completed.\n")
    print("✅ Experiments completed:", outdir)


if __name__ == "__main__":
    main()
