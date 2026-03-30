import json
from pathlib import Path

import numpy as np
import pytest

try:
    from whatshap import core  # noqa: F401
    HAS_WHCORE = True
except ImportError:
    HAS_WHCORE = False

pytestmark = pytest.mark.skipif(
    not HAS_WHCORE,
    reason="whatshap.core not importable (vendored core not built/installed).",
)

if not HAS_WHCORE:
    pytest.skip(
        "whatshap.core not importable (vendored core not built/installed).",
        allow_module_level=True,
    )

from algorithms.diploid import whatshap_driver


def _write_minimal_vcf(vcf_path: Path, gts: list[str], sample: str = "SAMPLE") -> None:
    """
    Write a minimal single-sample VCF with GT only.
    gts length == number of variant records.
    """
    assert all(gt in {"0/0", "0/1", "1/1", "./."} for gt in gts)

    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n")
        for i, gt in enumerate(gts):
            pos = i + 1
            f.write(f"chr1\t{pos}\t.\tA\tC\t.\tPASS\t.\tGT\t{gt}\n")


def _parse_sample_fields(line: str) -> dict[str, str]:
    fields = line.rstrip("\n").split("\t")
    fmt_keys = fields[8].split(":")
    sample_vals = fields[9].split(":")
    d = dict(zip(fmt_keys, sample_vals))
    return d


def test_pedigreedptable_phases_two_hets_same_phase_set(tmp_path: Path):
    # 10 variants total; het at indices 5 and 7 only
    gts = ["0/0"] * 10
    gts[5] = "0/1"
    gts[7] = "0/1"

    vcf_in = tmp_path / "in.vcf"
    _write_minimal_vcf(vcf_in, gts)

    # Dense reads matrix (R x N) with missing=-1 outside covered sites.
    # Two reads that connect variant 5 and 7:
    #  - read0 supports allele 0 at both
    #  - read1 supports allele 1 at both
    reads = np.full((2, 10), -1, dtype=int)
    reads[0, 5] = 0
    reads[0, 7] = 0
    reads[1, 5] = 1
    reads[1, 7] = 1

    npz_in = tmp_path / "reads.npz"
    np.savez_compressed(npz_in, reads=reads)

    out_prefix = tmp_path / "out"

    # Build args object (what your CLI effectively passes)
    args = type("Args", (), {})()
    args.input = str(npz_in)
    args.output_prefix = str(out_prefix)
    args.max_coverage = 10
    args.error_rate = 0.0
    args.vcf = str(vcf_in)
    args.sample = None
    args.output_vcf = None
    args.solver = "whatshap"       # <-- force PedigreeDPTable path
    args.recomb_rate = 1.26

    whatshap_driver.main(args)

    vcf_out = Path(str(out_prefix) + ".phased.vcf")
    summary_path = Path(str(out_prefix) + ".summary.json")
    hap_path = Path(str(out_prefix) + ".haplotypes.tsv")

    assert vcf_out.exists()
    assert summary_path.exists()
    assert hap_path.exists()

    # 1) Confirm solver recorded as whatshap
    summary = json.loads(summary_path.read_text())
    assert summary.get("solver") == "whatshap"

    # 2) Parse output VCF: variants 5 and 7 should be phased with PS=6
    records = [ln for ln in vcf_out.read_text().splitlines() if ln and not ln.startswith("#")]
    assert len(records) == 10

    # Homozygous sites remain unphased; het sites become phased with '|'
    het5 = _parse_sample_fields(records[5])
    het7 = _parse_sample_fields(records[7])

    assert "|" in het5["GT"], f"Expected phased GT at index 5, got {het5['GT']}"
    assert "|" in het7["GT"], f"Expected phased GT at index 7, got {het7['GT']}"

    # WhatsHap-style PS: leftmost variant index in component + 1
    # The component is {5,7}, so PS should be 6
    assert het5.get("PS") == "6"
    assert het7.get("PS") == "6"

    # 3) Also sanity-check haplotypes.tsv columns at 5 and 7 contain both alleles across the 2 haplotypes
    # Format: two lines: "0\t<bitstring>" and "1\t<bitstring>"
    rows = [ln.split("\t")[1].strip() for ln in hap_path.read_text().splitlines() if ln.strip()]
    assert len(rows) == 2
    assert len(rows[0]) == 10 and len(rows[1]) == 10

    col5 = {rows[0][5], rows[1][5]}
    col7 = {rows[0][7], rows[1][7]}
    assert col5 == {"0", "1"}
    assert col7 == {"0", "1"}
