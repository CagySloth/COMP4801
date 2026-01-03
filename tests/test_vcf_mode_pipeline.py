# tests/test_vcf_mode_pipeline.py

import subprocess
from pathlib import Path
import numpy as np

from algorithms.eval.metrics import hap_truth_accuracy


def _read_haplotypes_tsv(path: Path) -> np.ndarray:
    """
    Your haplotypes TSV format:
      <row_index>\t<bitstring>\n
    Returns: array shape (2, N) dtype=int
    """
    rows = []
    with path.open("r") as f:
        for line in f:
            parts = line.strip().split("\t")
            assert len(parts) == 2, f"Unexpected hap TSV line: {line!r}"
            bitstring = parts[1]
            rows.append([int(c) for c in bitstring])
    arr = np.array(rows, dtype=int)
    assert arr.ndim == 2 and arr.shape[0] == 2, f"Expected (2,N) haplotypes, got {arr.shape}"
    return arr


def _count_het_and_phased(vcf_path: Path) -> tuple[int, int]:
    """
    Counts heterozygous GTs among {0/1,1/0,0|1,1|0} and how many are phased (contain '|').
    Assumes single-sample VCF (sample in column 10).
    """
    het = 0
    phased = 0

    with vcf_path.open("r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            sample_field = fields[9]
            gt = sample_field.split(":")[0]

            if gt in ("0/1", "1/0", "0|1", "1|0"):
                het += 1
                if "|" in gt:
                    phased += 1

    return het, phased


def test_vcf_mode_end_to_end(tmp_path: Path) -> None:
    """
    End-to-end regression:
      simulate -> produces {prefix}.vcf and {prefix}.reads.npz
      phase (VCF-mode) -> produces {out_prefix}.phased.vcf and {out_prefix}.haplotypes.tsv
      checks: phased VCF exists, has at least one phased het, accuracy is sane.
    """
    prefix = tmp_path / "demo"
    out_prefix = tmp_path / "demo_phased"

    # 1) Simulate (must write demo.vcf, demo.reads.npz, demo.haplotypes.tsv)
    subprocess.check_call(
        [
            "python",
            "-m",
            "dataset.simulate",
            "-p",
            "2",
            "-n",
            "100",
            "-r",
            "50",
            "-l",
            "30",
            "-e",
            "0.01",
            "-m",
            "0.0",
            "-o",
            str(prefix),
        ]
    )

    in_vcf = Path(str(prefix) + ".vcf")
    in_npz = Path(str(prefix) + ".reads.npz")
    truth_tsv = Path(str(prefix) + ".haplotypes.tsv")

    assert in_vcf.exists(), f"Simulator did not write VCF: {in_vcf}"
    assert in_npz.exists(), f"Simulator did not write reads NPZ: {in_npz}"
    assert truth_tsv.exists(), f"Simulator did not write truth haplotypes TSV: {truth_tsv}"

    # 2) Phase in VCF-mode
    subprocess.check_call(
        [
            "python",
            "-m",
            "algorithms.cli.phase",
            "diploid-whats",
            "-i",
            str(in_npz),
            "--vcf",
            str(in_vcf),
            "--output-prefix",
            str(out_prefix),
        ]
    )

    out_vcf = Path(str(out_prefix) + ".phased.vcf")
    pred_tsv = Path(str(out_prefix) + ".haplotypes.tsv")

    assert out_vcf.exists(), f"Phasing did not write phased VCF: {out_vcf}"
    assert pred_tsv.exists(), f"Phasing did not write predicted haplotypes TSV: {pred_tsv}"

    # 3) VCF sanity: input has hets, output phases at least one het
    het_in, phased_in = _count_het_and_phased(in_vcf)
    het_out, phased_out = _count_het_and_phased(out_vcf)

    assert het_in > 0, "Input VCF has no heterozygous sites; cannot test phasing."
    assert phased_in == 0, f"Input VCF unexpectedly already phased: phased_in={phased_in}"
    assert het_out == het_in, "Het count changed between input and output VCF unexpectedly."
    assert phased_out > 0, "Output VCF did not phase any heterozygous sites (no '|' in GT)."

    # 4) Accuracy sanity (full-length haplotypes reconstructed)
    truth = _read_haplotypes_tsv(truth_tsv)
    pred = _read_haplotypes_tsv(pred_tsv)
    assert truth.shape == pred.shape, f"Truth/pred shape mismatch: {truth.shape} vs {pred.shape}"

    acc = hap_truth_accuracy(pred, truth)

    # Conservative threshold: this should be comfortably > 0.8 in VCF-mode
    # If you later increase error/missing rates, you can lower this threshold for that test.
    assert acc > 0.80, f"Accuracy too low for VCF-mode sanity run: {acc}"
