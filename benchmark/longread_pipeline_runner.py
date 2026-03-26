# benchmark/longread_pipeline_runner.py
import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path
import gzip
from typing import Dict, Any, Optional


def _which_or_fail(tool: str) -> None:
    if shutil.which(tool) is None:
        raise RuntimeError(f"Missing system tool '{tool}'. Install it (apt/brew/conda) and retry.")
    

def _bcftools_filter_snps(in_vcfgz: str, out_vcfgz: str) -> None:
    _run(["bcftools", "view", "-v", "snps", "-m2", "-M2", "-Oz", "-o", out_vcfgz, in_vcfgz])
    _run(["bcftools", "index", "-t", out_vcfgz])


def _run(cmd: list[str]) -> None:
    print("[RUN]", " ".join(cmd))
    subprocess.check_call(cmd)


def _piped(cmd1: list[str], cmd2: list[str]) -> None:
    print("[RUN]", " ".join(cmd1), "|", " ".join(cmd2))
    p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
    assert p1.stdout is not None
    p2 = subprocess.Popen(cmd2, stdin=p1.stdout)
    p1.stdout.close()
    rc2 = p2.wait()
    rc1 = p1.wait()
    if rc1 != 0:
        raise subprocess.CalledProcessError(rc1, cmd1)
    if rc2 != 0:
        raise subprocess.CalledProcessError(rc2, cmd2)


def _open_text(path: str, mode: str):
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode)


def _count_vcf_records(path: str) -> int:
    c = 0
    with _open_text(path, "r") as f:
        for line in f:
            if line and not line.startswith("#"):
                c += 1
    return c


def _count_het_snps(path: str) -> int:
    """Count heterozygous SNP records in a single-sample VCF."""
    c = 0
    with _open_text(path, "r") as f:
        sample_col = None
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                sample_col = 9 if len(header) >= 10 else None
                continue
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10 or sample_col is None:
                continue
            ref = fields[3]
            alt = fields[4].split(",")[0]
            if len(ref) != 1 or len(alt) != 1:
                continue

            fmt = fields[8].split(":")
            sval = fields[sample_col].split(":")
            d = dict(zip(fmt, sval))
            gt = d.get("GT", "./.")
            if gt in ("./.", ".|.", "."):
                continue
            sep = "|" if "|" in gt else "/"
            a_str, b_str = gt.split(sep)
            if a_str == "." or b_str == ".":
                continue
            if a_str != b_str:
                c += 1
    return c


def _vcf_site_set(path: str):
    s = set()
    with _open_text(path, "r") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            pos1 = int(fields[1])
            ref = fields[3]
            alt = fields[4].split(",")[0]
            if len(ref) == 1 and len(alt) == 1:
                s.add((chrom, pos1, ref, alt))
    return s


def _write_oracle_vcf(truth_vcfgz: str, oracle_vcf: str) -> None:
    """Create an oracle VCF from truth:
    - same sites and GT
    - GT: 0|1 -> 0/1
    - remove PS if present
    """
    with _open_text(truth_vcfgz, "r") as fin, open(oracle_vcf, "w") as fout:
        for line in fin:
            if line.startswith("##FORMAT=<ID=PS"):
                continue
            if line.startswith("#"):
                fout.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                fout.write(line)
                continue

            fmt = fields[8].split(":")
            samp = fields[9].split(":")
            if len(samp) < len(fmt):
                samp += ["."] * (len(fmt) - len(samp))

            gt_i = None
            ps_i = None
            for i, k in enumerate(fmt):
                if k == "GT":
                    gt_i = i
                elif k == "PS":
                    ps_i = i

            if gt_i is not None:
                samp[gt_i] = samp[gt_i].replace("|", "/")

            if ps_i is not None:
                new_fmt = []
                new_samp = []
                for i in range(len(fmt)):
                    if i == ps_i:
                        continue
                    new_fmt.append(fmt[i])
                    new_samp.append(samp[i] if i < len(samp) else ".")
                fields[8] = ":".join(new_fmt)
                fields[9] = ":".join(new_samp)
            else:
                fields[9] = ":".join(samp[: len(fmt)])

            fout.write("\t".join(fields) + "\n")


def _load_json(path: str) -> Optional[Dict[str, Any]]:
    try:
        with open(path, "r") as f:
            return json.load(f)
    except FileNotFoundError:
        return None


def main(args=None):
    ap = argparse.ArgumentParser(description="Run long-read phasing pipeline end-to-end (Steps 1-7).")

    ap.add_argument("--prefix", required=True, help="Output prefix, e.g. output/exp1/demo")
    ap.add_argument("--seed", type=int, default=0)
    
    ap.add_argument("--ref-preset", choices=["plain", "toy", "realistic"], default="plain",
                    help="Reference complexity preset for dataset.longread.reference")
    
    ap.add_argument("--dup-segments", type=int, default=0)
    ap.add_argument("--dup-len", type=int, default=3000)
    ap.add_argument("--dup-min-gap", type=int, default=500)

    # Reference + truth
    ap.add_argument("--ref-length", type=int, default=20000)
    ap.add_argument("--num-snps", type=int, default=200)
    ap.add_argument("--het-rate", type=float, default=0.8)
    ap.add_argument("--avoid-regions", action="store_true")
    ap.add_argument("--num-indels", type=int, default=0)
    ap.add_argument("--indel-min-len", type=int, default=1)
    ap.add_argument("--indel-max-len", type=int, default=5)
    ap.add_argument("--indel-het-rate", type=float, default=0.5)

    # Reads
    ap.add_argument("--num-reads", type=int, default=200)
    ap.add_argument("--min-len", type=int, default=2000)
    ap.add_argument("--max-len", type=int, default=6000)
    ap.add_argument("--platform", choices=["ont", "perfect"], default="ont")
    ap.add_argument("--ont-profile", choices=["classic", "q20"], default="q20")
    
    ap.add_argument("--len-model", choices=["uniform","lognormal"], default="uniform")
    ap.add_argument("--ln-mean", type=float, default=8.3)
    ap.add_argument("--ln-sigma", type=float, default=0.6)
    ap.add_argument("--start-model", choices=["uniform","dropout"], default="uniform")
    ap.add_argument("--dropout-fraction", type=float, default=0.0)
    ap.add_argument("--dropout-block-len", type=int, default=1000)

    # Align/call
    ap.add_argument("--map-preset", default="map-ont")
    ap.add_argument("--call-min-mapq", type=int, default=20, help="bcftools mpileup -q")
    ap.add_argument("--call-min-baseq", type=int, default=15, help="bcftools mpileup -Q")

    # Which VCF(s) to phase
    ap.add_argument("--vcf-source", choices=["called", "oracle", "both"], default="both")

    # Phase (vendored)
    ap.add_argument("--max-coverage", type=int, default=15)
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-baseq", type=int, default=20)
    
    # Burst errors
    ap.add_argument("--burst-prob", type=float, default=0.0)
    ap.add_argument("--burst-count", type=int, default=1)
    ap.add_argument("--burst-len", type=int, default=200)
    ap.add_argument("--burst-mult", type=float, default=5.0)
    
    ap.add_argument(
        "--phase-snps-only",
        action="store_true",
        help="If set, filter the phasing input VCF to biallelic SNPs only (-v snps -m2 -M2). "
            "Recommended when truth includes indels."
    )
    ap.add_argument(
        "--eval-snps-only",
        action="store_true",
        help="If set, filter truth VCF to biallelic SNPs for evaluation. "
            "Recommended together with --phase-snps-only."
    )

    if args is None:
        args = ap.parse_args()

    want_called = args.vcf_source in ("called", "both")
    want_oracle = args.vcf_source in ("oracle", "both")

    _which_or_fail("minimap2")
    _which_or_fail("samtools")
    _which_or_fail("bcftools")

    prefix = Path(args.prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    py = sys.executable

    # Filenames
    ref_fa = str(prefix) + ".ref.fasta"
    ref_meta = str(prefix) + ".ref.meta.json"

    truth_vcf = str(prefix) + ".truth.vcf"
    truth_vcfgz = truth_vcf + ".gz"
    truth_eval_vcfgz = truth_vcfgz

    oracle_vcf = str(prefix) + ".oracle.vcf"
    oracle_vcfgz = oracle_vcf + ".gz"

    hap1 = str(prefix) + ".hap1.fasta"
    hap2 = str(prefix) + ".hap2.fasta"

    reads_fq = str(prefix) + ".reads.fastq"
    bam = str(prefix) + ".bam"

    called_vcfgz = str(prefix) + ".called.vcf.gz"

    # Phasing outputs (keep .ws for called)
    phase_prefix_called = str(prefix) + ".ws"
    phased_vcf_called = phase_prefix_called + ".phased.vcf"
    eval_json_called = phase_prefix_called + ".eval.json"

    phase_prefix_oracle = str(prefix) + ".ws_oracle"
    phased_vcf_oracle = phase_prefix_oracle + ".phased.vcf"
    eval_json_oracle = phase_prefix_oracle + ".eval.json"

    pipeline_json = str(prefix) + ".pipeline.json"

    t0 = time.perf_counter()

    # Step 1: reference
    _run([
        py, "-m", "dataset.longread.reference",
        "-o", str(prefix),
        "--length", str(args.ref_length),
        "--seed", str(args.seed),
        "--preset", str(args.ref_preset),
        "--dup-segments", str(args.dup_segments),
        "--dup-len", str(args.dup_len),
        "--dup-min-gap", str(args.dup_min_gap),
    ])

    # Step 2: truth
    truth_cmd = [
        py, "-m", "dataset.longread.truth",
        "--ref", ref_fa,
        "-o", str(prefix),
        "--seed", str(args.seed),
        "--num-snps", str(args.num_snps),
        "--het-rate", str(args.het_rate),
        "--phased-truth",
        "--random-phase",
        "--num-indels", str(args.num_indels),
        "--indel-min-len", str(args.indel_min_len),
        "--indel-max-len", str(args.indel_max_len),
        "--indel-het-rate", str(args.indel_het_rate),
    ]
    if args.avoid_regions:
        truth_cmd += ["--ref-meta", ref_meta, "--avoid-regions"]
    _run(truth_cmd)

    _run(["bcftools", "view", "-Oz", "-o", truth_vcfgz, truth_vcf])
    _run(["bcftools", "index", "-t", truth_vcfgz])
    
    # Optional: make SNP-only truth for evaluation (recommended if truth has indels)
    if args.phase_snps_only or args.eval_snps_only:
        truth_eval_vcfgz = str(prefix) + ".truth.snps.vcf.gz"
        _bcftools_filter_snps(truth_vcfgz, truth_eval_vcfgz)
    else:
        truth_eval_vcfgz = truth_vcfgz

    truth_het_snps = _count_het_snps(truth_vcfgz)

    # Oracle VCF
    if want_oracle:
        _write_oracle_vcf(truth_vcfgz, oracle_vcf)
        _run(["bcftools", "view", "-Oz", "-o", oracle_vcfgz, oracle_vcf])
        _run(["bcftools", "index", "-t", oracle_vcfgz])

    # Step 3: reads
    readsim_cmd = [
        py, "-m", "dataset.longread.readsim",
        "--hap1", hap1,
        "--hap2", hap2,
        "-o", str(prefix),
        "--seed", str(args.seed),
        "--num-reads", str(args.num_reads),
        "--min-len", str(args.min_len),
        "--max-len", str(args.max_len),
        "--hap1-frac", "0.5",
        "--platform", str(args.platform),
        "--len-model", str(args.len_model),
        "--ln-mean", str(args.ln_mean),
        "--ln-sigma", str(args.ln_sigma),
        "--start-model", str(args.start_model),
        "--burst-prob", str(args.burst_prob),
        "--burst-count", str(args.burst_count),
        "--burst-len", str(args.burst_len),
        "--burst-mult", str(args.burst_mult),
    ]
    if args.platform == "ont":
        readsim_cmd += ["--ont-profile", str(args.ont_profile)]
    if args.start_model == "dropout":
        readsim_cmd += ["--dropout-fraction", str(args.dropout_fraction),
                        "--dropout-block-len", str(args.dropout_block_len)]    
    _run(readsim_cmd)

    # Step 4: align
    _piped(
        ["minimap2", "-a", "-x", str(args.map_preset), ref_fa, reads_fq],
        ["samtools", "sort", "-o", bam],
    )
    _run(["samtools", "index", bam])

    # Step 5: call (optional)
    if want_called:
        p1 = ["bcftools", "mpileup", "-Ou", "-f", ref_fa, "-q", str(args.call_min_mapq), "-Q", str(args.call_min_baseq), bam]
        p2 = ["bcftools", "call", "-mv", "-Oz", "-o", called_vcfgz]
        _piped(p1, p2)
        _run(["bcftools", "index", "-t", called_vcfgz])

    # Step 6-7: phase + eval
    phase_runs: Dict[str, Dict[str, Any]] = {}

    def _phase_and_eval(tag: str, vcf_in: str, out_prefix: str, out_vcf: str, out_eval: str) -> None:
        vcf_for_phasing = vcf_in
        if args.phase_snps_only:
            vcf_for_phasing = str(prefix) + f".{tag}.snps.vcf.gz"
            _bcftools_filter_snps(vcf_in, vcf_for_phasing)
            
        _run([
            py, "-m", "algorithms.cli.phase", "diploid-whats-bam",
            "--bam", bam,
            "--vcf", vcf_for_phasing,
            "--output-prefix", out_prefix,
            "--output-vcf", out_vcf,
            "--max-coverage", str(args.max_coverage),
            "--min-mapq", str(args.min_mapq),
            "--min-baseq", str(args.min_baseq),
        ])
        _run([py, "-m", "benchmark.vcf_phase_eval", "--truth", truth_eval_vcfgz, "--pred", out_vcf, "--out", out_eval])

        phase_runs[tag] = {
            "vcf_input": vcf_in,
            "vcf_input_for_phasing": vcf_for_phasing,
            "phased_vcf": out_vcf,
            "eval_json": out_eval,
            "summary_json": out_prefix + ".summary.json",
            "eval": _load_json(out_eval),
        }

        ev = phase_runs[tag]["eval"]
        if isinstance(ev, dict):
            phase_sets = ev.get("phase_sets", {}) or {}
            correct = sum(int(st.get("best_match", 0)) for st in phase_sets.values() if isinstance(st, dict))

            shared_het = int(ev.get("shared_het_records", 0) or 0)
            phased_het = int(ev.get("pred_phased_het_records", 0) or 0)

            phase_runs[tag]["derived"] = {
                "truth_het_snps": int(truth_het_snps),
                "shared_het_snps": shared_het,
                "phased_het_snps": phased_het,
                "correct_phased_het_snps_bestflip": int(correct),
                "phasing_rate_on_shared_het": (phased_het / shared_het) if shared_het else None,
                "effective_phased_recall": (correct / truth_het_snps) if truth_het_snps else None,
                "shared_het_recall": (shared_het / truth_het_snps) if truth_het_snps else None,
            }

    if want_called:
        _phase_and_eval("called", called_vcfgz, phase_prefix_called, phased_vcf_called, eval_json_called)
    if want_oracle:
        _phase_and_eval("oracle", oracle_vcfgz, phase_prefix_oracle, phased_vcf_oracle, eval_json_oracle)

    # Calling metrics (only for called)
    truth_sites = _vcf_site_set(truth_vcfgz)
    truth_n = len(truth_sites)
    callset: Dict[str, Any] = {
        "truth_snps": truth_n,
        "truth_het_snps": int(truth_het_snps),
        "called_snps": None,
        "shared_snps": None,
        "precision": None,
        "recall": None,
        "notes": "Precision/recall only computed when --vcf-source includes called.",
    }

    if want_called:
        called_sites = _vcf_site_set(called_vcfgz)
        shared_sites = truth_sites & called_sites
        called_n = len(called_sites)
        shared_n = len(shared_sites)
        callset.update({
            "called_snps": called_n,
            "shared_snps": shared_n,
            "precision": (shared_n / called_n) if called_n else None,
            "recall": (shared_n / truth_n) if truth_n else None,
        })

    total_sec = time.perf_counter() - t0

    counts_raw: Dict[str, Any] = {
        "truth_records": _count_vcf_records(truth_vcfgz),
        "called_records": _count_vcf_records(called_vcfgz) if want_called else None,
        "oracle_records": _count_vcf_records(oracle_vcfgz) if want_oracle else None,
        "phasing": {tag: {"phased_records": _count_vcf_records(info["phased_vcf"])} for tag, info in phase_runs.items()},
    }

    report: Dict[str, Any] = {
        "prefix": str(prefix),
        "seed": int(args.seed),
        "vcf_source": str(args.vcf_source),
        "steps": {
            "ref_fasta": ref_fa,
            "truth_vcf_gz": truth_vcfgz,
            "truth_eval_vcf_gz": truth_eval_vcfgz,
            "oracle_vcf_gz": oracle_vcfgz if want_oracle else None,
            "reads_fastq": reads_fq,
            "bam": bam,
            "called_vcf_gz": called_vcfgz if want_called else None,
        },
        "callset": callset,
        "phasing_runs": phase_runs,
        "counts_raw": counts_raw,
        "params": {
            "ref_preset": str(args.ref_preset),
            "ref_length": int(args.ref_length),
            "num_snps": int(args.num_snps),
            "het_rate": float(args.het_rate),
            "avoid_regions": bool(args.avoid_regions),
            "num_reads": int(args.num_reads),
            "min_len": int(args.min_len),
            "max_len": int(args.max_len),
            "platform": str(args.platform),
            "ont_profile": str(args.ont_profile) if args.platform == "ont" else None,
            "map_preset": str(args.map_preset),
            "call_min_mapq": int(args.call_min_mapq),
            "call_min_baseq": int(args.call_min_baseq),
            "max_coverage": int(args.max_coverage),
            "min_mapq": int(args.min_mapq),
            "min_baseq": int(args.min_baseq),
            "dup_segments": int(args.dup_segments),
            "dup_len": int(args.dup_len),
            "dup_min_gap": int(args.dup_min_gap),
            "num_indels": int(args.num_indels),
            "indel_min_len": int(args.indel_min_len),
            "indel_max_len": int(args.indel_max_len),
            "indel_het_rate": float(args.indel_het_rate),

            "len_model": str(args.len_model),
            "ln_mean": float(args.ln_mean) if args.len_model == "lognormal" else None,
            "ln_sigma": float(args.ln_sigma) if args.len_model == "lognormal" else None,
            "start_model": str(args.start_model),
            "dropout_fraction": float(args.dropout_fraction) if args.start_model == "dropout" else 0.0,
            "dropout_block_len": int(args.dropout_block_len) if args.start_model == "dropout" else None,
            "burst_prob": float(args.burst_prob),
            "burst_count": int(args.burst_count),
            "burst_len": int(args.burst_len),
            "burst_mult": float(args.burst_mult),
            "phase_snps_only": bool(args.phase_snps_only),
            "eval_snps_only": bool(args.eval_snps_only or args.phase_snps_only),
        },
        "time_total_sec": total_sec,
    }

    with open(pipeline_json, "w") as f:
        json.dump(report, f, indent=2)

    print(f"✅ Pipeline report: {pipeline_json}")


if __name__ == "__main__":
    main()
