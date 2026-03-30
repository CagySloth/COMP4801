from __future__ import annotations

import argparse
import json
from pathlib import Path

import pytest

from benchmark import vcf_phase_eval


def _write_vcf(path: Path, records: list[tuple[str, int, str, str, str, str]]) -> None:
    """
    records: list of (chrom, pos1, ref, alt, fmt, sample_value)
    fmt must include GT (and PS if you want phase sets).
    """
    lines = [
        "##fileformat=VCFv4.2",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
    ]
    for chrom, pos1, ref, alt, fmt, sample in records:
        lines.append(f"{chrom}\t{pos1}\t.\t{ref}\t{alt}\t.\tPASS\t.\t{fmt}\t{sample}")
    path.write_text("\n".join(lines) + "\n")


def _run_eval_capture(capsys, truth: Path, pred: Path) -> dict:
    args = argparse.Namespace(truth=str(truth), pred=str(pred), sample=None, out=None)
    vcf_phase_eval.main(args)
    out = capsys.readouterr().out.strip()
    return json.loads(out)


def test_parse_gt_basic():
    a, b, phased = vcf_phase_eval._parse_gt("0|1")
    assert (a, b, phased) == (0, 1, True)

    a, b, phased = vcf_phase_eval._parse_gt("1/0")
    assert (a, b, phased) == (1, 0, False)

    a, b, phased = vcf_phase_eval._parse_gt("./.")
    assert a is None and b is None and phased is False


def test_blockflip_and_counts(capsys, tmp_path: Path):
    """
    Truth: 4 het SNPs (all in PS=10), fully phased.
    Pred: 3 phased het (perfectly flipped) + 1 unphased het.
    Expect best flip chosen, accuracy 1.0 on phased het, switch error 0.0.
    """
    truth = tmp_path / "truth.vcf"
    pred = tmp_path / "pred.vcf"

    truth_recs = [
        ("chr1", 100, "A", "C", "GT:PS", "0|1:10"),
        ("chr1", 200, "A", "C", "GT:PS", "0|1:10"),
        ("chr1", 300, "A", "C", "GT:PS", "1|0:10"),
        ("chr1", 400, "A", "C", "GT:PS", "1|0:10"),
    ]
    pred_recs = [
        ("chr1", 100, "A", "C", "GT:PS", "1|0:10"),
        ("chr1", 200, "A", "C", "GT:PS", "1|0:10"),
        ("chr1", 300, "A", "C", "GT:PS", "0|1:10"),
        ("chr1", 400, "A", "C", "GT:PS", "0/1:10"),  # unphased het
    ]

    _write_vcf(truth, truth_recs)
    _write_vcf(pred, pred_recs)

    rep = _run_eval_capture(capsys, truth, pred)

    assert rep["shared_snp_records"] == 4
    assert rep["shared_het_records"] == 4
    assert rep["pred_phased_het_records"] == 3
    assert rep["pred_unphased_het_records"] == 1

    assert rep["phase_accuracy_blockflip"] == pytest.approx(1.0, abs=1e-9)

    # 3 phased het => 2 adjacency comparisons
    assert rep["switch_den"] == 2
    assert rep["switch_error_rate"] == pytest.approx(0.0, abs=1e-9)

    assert rep["num_phase_sets"] == 1
    ps = next(iter(rep["phase_sets"].keys()))
    assert rep["phase_sets"][ps]["best_flip"] is True


def test_switch_error_nonzero(capsys, tmp_path: Path):
    """
    Force a non-zero switch error case.

    Truth:
      100 0|1
      200 0|1
      300 1|0
      400 1|0   (PS=10)

    Pred:
      100 1|0
      200 1|0
      300 1|0
      400 0|1   (PS=10)

    This yields 2 switches / 3 adjacencies after best flip.
    """
    truth = tmp_path / "truth.vcf"
    pred = tmp_path / "pred.vcf"

    truth_recs = [
        ("chr1", 100, "A", "C", "GT:PS", "0|1:10"),
        ("chr1", 200, "A", "C", "GT:PS", "0|1:10"),
        ("chr1", 300, "A", "C", "GT:PS", "1|0:10"),
        ("chr1", 400, "A", "C", "GT:PS", "1|0:10"),
    ]
    pred_recs = [
        ("chr1", 100, "A", "C", "GT:PS", "1|0:10"),
        ("chr1", 200, "A", "C", "GT:PS", "1|0:10"),
        ("chr1", 300, "A", "C", "GT:PS", "1|0:10"),
        ("chr1", 400, "A", "C", "GT:PS", "0|1:10"),
    ]

    _write_vcf(truth, truth_recs)
    _write_vcf(pred, pred_recs)

    rep = _run_eval_capture(capsys, truth, pred)

    assert rep["shared_het_records"] == 4
    assert rep["pred_phased_het_records"] == 4
    assert rep["switch_den"] == 3
    assert rep["switch_error_rate"] == pytest.approx(2 / 3, abs=1e-9)
    assert rep["phase_accuracy_blockflip"] == pytest.approx(3 / 4, abs=1e-9)