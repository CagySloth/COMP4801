# benchmark/longread_pipeline_runner.py
import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path
import gzip


def _which_or_fail(tool: str) -> None:
    if shutil.which(tool) is None:
        raise RuntimeError(f"Missing system tool '{tool}'. Install via apt and retry.")


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


def main(args=None):
    ap = argparse.ArgumentParser(description="Run long-read phasing pipeline end-to-end (Steps 1-7).")

    ap.add_argument("--prefix", required=True, help="Output prefix, e.g. output/exp1/demo")
    ap.add_argument("--seed", type=int, default=0)

    # Reference + truth
    ap.add_argument("--ref-length", type=int, default=20000)
    ap.add_argument("--num-snps", type=int, default=200)
    ap.add_argument("--het-rate", type=float, default=0.8)
    ap.add_argument("--avoid-regions", action="store_true")

    # Reads
    ap.add_argument("--num-reads", type=int, default=200)
    ap.add_argument("--min-len", type=int, default=2000)
    ap.add_argument("--max-len", type=int, default=6000)

    # Align/call
    ap.add_argument("--map-preset", default="map-ont")

    # Phase (vendored)
    ap.add_argument("--max-coverage", type=int, default=15)
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-baseq", type=int, default=20)

    if args is None:
        args = ap.parse_args()

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
    hap1 = str(prefix) + ".hap1.fasta"
    hap2 = str(prefix) + ".hap2.fasta"

    reads_fq = str(prefix) + ".reads.fastq"
    bam = str(prefix) + ".bam"

    called_vcfgz = str(prefix) + ".called.vcf.gz"

    phase_prefix = str(prefix) + ".ws"
    phased_vcf = phase_prefix + ".phased.vcf"
    eval_json = phase_prefix + ".eval.json"
    pipeline_json = str(prefix) + ".pipeline.json"

    t0 = time.perf_counter()

    # Step 1: reference
    _run([py, "-m", "dataset.longread.reference", "-o", str(prefix), "--length", str(args.ref_length), "--seed", str(args.seed)])

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
    ]
    if args.avoid_regions:
        truth_cmd += ["--ref-meta", ref_meta, "--avoid-regions"]
    _run(truth_cmd)

    # compress truth for convenience
    _run(["bcftools", "view", "-Oz", "-o", truth_vcfgz, truth_vcf])
    _run(["bcftools", "index", "-t", truth_vcfgz])

    # Step 3: reads (ONT-like)
    _run([
        py, "-m", "dataset.longread.readsim",
        "--hap1", hap1,
        "--hap2", hap2,
        "-o", str(prefix),
        "--seed", str(args.seed),
        "--num-reads", str(args.num_reads),
        "--min-len", str(args.min_len),
        "--max-len", str(args.max_len),
        "--hap1-frac", "0.5",
        "--platform", "ont",
        "--ont-profile", "classic",
    ])

    # # Step 3: reads (perfect reads for now)
    # _run([
    #     py, "-m", "dataset.longread.readsim",
    #     "--hap1", hap1,
    #     "--hap2", hap2,
    #     "-o", str(prefix),
    #     "--seed", str(args.seed),
    #     "--num-reads", str(args.num_reads),
    #     "--min-len", str(args.min_len),
    #     "--max-len", str(args.max_len),
    #     "--hap1-frac", "0.5",
    # ])

    # Step 4: align -> BAM
    _piped(
        ["minimap2", "-a", "-x", str(args.map_preset), ref_fa, reads_fq],
        ["samtools", "sort", "-o", bam],
    )
    _run(["samtools", "index", bam])

    # Step 5: call variants
    p1 = ["bcftools", "mpileup", "-Ou", "-f", ref_fa, bam]
    p2 = ["bcftools", "call", "-mv", "-Oz", "-o", called_vcfgz]
    _piped(p1, p2)
    _run(["bcftools", "index", "-t", called_vcfgz])

    # Step 6: phase (vendored core via your CLI)
    _run([
        py, "-m", "algorithms.cli.phase", "diploid-whats-bam",
        "--bam", bam,
        "--vcf", called_vcfgz,
        "--output-prefix", phase_prefix,
        "--output-vcf", phased_vcf,
        "--max-coverage", str(args.max_coverage),
        "--min-mapq", str(args.min_mapq),
        "--min-baseq", str(args.min_baseq),
    ])

    # Step 7: evaluate phase vs truth
    _run([py, "-m", "benchmark.vcf_phase_eval", "--truth", truth_vcfgz, "--pred", phased_vcf, "--out", eval_json])

    # Quick precision/recall for calling
    truth_sites = _vcf_site_set(truth_vcfgz)
    called_sites = _vcf_site_set(called_vcfgz)
    shared_sites = truth_sites & called_sites

    truth_n = len(truth_sites)
    called_n = len(called_sites)
    shared_n = len(shared_sites)

    precision = (shared_n / called_n) if called_n else None
    recall = (shared_n / truth_n) if truth_n else None

    total_sec = time.perf_counter() - t0

    report = {
        "prefix": str(prefix),
        "seed": int(args.seed),
        "steps": {
            "ref_fasta": ref_fa,
            "truth_vcf_gz": truth_vcfgz,
            "reads_fastq": reads_fq,
            "bam": bam,
            "called_vcf_gz": called_vcfgz,
            "phased_vcf": phased_vcf,
            "eval_json": eval_json,
        },
        "callset": {
            "truth_snps": truth_n,
            "called_snps": called_n,
            "shared_snps": shared_n,
            "precision": precision,
            "recall": recall,
        },
        "counts_raw": {
            "truth_records": _count_vcf_records(truth_vcfgz),
            "called_records": _count_vcf_records(called_vcfgz),
            "phased_records": _count_vcf_records(phased_vcf),
        },
        "params": {
            "ref_length": int(args.ref_length),
            "num_snps": int(args.num_snps),
            "het_rate": float(args.het_rate),
            "avoid_regions": bool(args.avoid_regions),
            "num_reads": int(args.num_reads),
            "min_len": int(args.min_len),
            "max_len": int(args.max_len),
            "map_preset": str(args.map_preset),
            "max_coverage": int(args.max_coverage),
            "min_mapq": int(args.min_mapq),
            "min_baseq": int(args.min_baseq),
        },
        "time_total_sec": total_sec,
    }

    with open(pipeline_json, "w") as f:
        json.dump(report, f, indent=2)

    print(f"✅ Pipeline report: {pipeline_json}")


if __name__ == "__main__":
    main()
