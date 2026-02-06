# benchmark/vcf_phase_eval.py
import argparse
import gzip
import json
from collections import defaultdict
from typing import Dict, Tuple, List, Optional


def _open_text(path: str, mode: str):
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode)


def _parse_gt(gt: str) -> Tuple[Optional[int], Optional[int], bool]:
    """
    Returns (a, b, phased) where a/b are allele ints or None.
    """
    if gt in ("./.", ".|.", "."):
        return None, None, False
    phased = "|" in gt
    sep = "|" if phased else "/"
    a_str, b_str = gt.split(sep)
    if a_str == "." or b_str == ".":
        return None, None, phased
    return int(a_str), int(b_str), phased


def read_vcf_minimal(path: str, sample: Optional[str] = None):
    """
    Minimal single-sample VCF reader.
    Returns dict keyed by (chrom, pos1) with:
      ref, alt, gt_str, (a,b,phased), ps_str
    """
    out = {}
    sample_col = None

    with _open_text(path, "r") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                if len(header) < 10:
                    raise ValueError("VCF has no sample columns")
                samples = header[9:]
                if sample is None:
                    sample_col = 9
                else:
                    if sample not in samples:
                        raise ValueError(f"Sample '{sample}' not in VCF. Available: {samples}")
                    sample_col = header.index(sample)
                continue
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            pos1 = int(fields[1])
            ref = fields[3]
            alt = fields[4]

            fmt = fields[8].split(":")
            sval = fields[sample_col].split(":")
            d = dict(zip(fmt, sval))
            gt = d.get("GT", "./.")
            ps = d.get("PS", ".")

            a, b, phased = _parse_gt(gt)
            out[(chrom, pos1)] = {
                "ref": ref,
                "alt": alt,
                "gt": gt,
                "a": a,
                "b": b,
                "phased": phased,
                "ps": ps,
            }

    return out


def main(args=None):
    ap = argparse.ArgumentParser(description="Evaluate phased VCF vs phased truth VCF (block-flip aware).")
    ap.add_argument("--truth", required=True, help="Truth VCF(.gz) with phased GT (0|1/1|0)")
    ap.add_argument("--pred", required=True, help="Predicted phased VCF(.gz)")
    ap.add_argument("--sample", default=None, help="Sample name (default: first sample)")
    ap.add_argument("--out", default=None, help="Write JSON report here (optional)")
    if args is None:
        args = ap.parse_args()

    truth = read_vcf_minimal(args.truth, sample=args.sample)
    pred = read_vcf_minimal(args.pred, sample=args.sample)

    shared = sorted(set(truth.keys()) & set(pred.keys()))
    shared_total = len(shared)

    # Only evaluate biallelic SNPs (REF/ALT length 1) for now
    shared = [k for k in shared if len(truth[k]["ref"]) == 1 and len(truth[k]["alt"]) == 1 and len(pred[k]["ref"]) == 1 and len(pred[k]["alt"]) == 1]
    shared_snp = len(shared)

    # Het sites in truth
    het_sites = []
    for k in shared:
        ta, tb = truth[k]["a"], truth[k]["b"]
        if ta is None or tb is None:
            continue
        if ta != tb:
            het_sites.append(k)

    shared_het = len(het_sites)

    # Partition predicted phased het sites by PS
    blocks = defaultdict(list)  # ps -> list[(chrom,pos1)]
    unphased_het = 0
    phased_het = 0

    for k in het_sites:
        pa, pb, pph = pred[k]["a"], pred[k]["b"], pred[k]["phased"]
        if pa is None or pb is None:
            continue
        if pph:
            phased_het += 1
            ps = pred[k]["ps"]
            if ps in (None, ".", ""):
                ps = "NO_PS"
            blocks[ps].append(k)
        else:
            unphased_het += 1

    # Compute per-block best flip accuracy and switch errors
    correct = 0
    total_phased_het = 0

    switches = 0
    switch_den = 0  # number of adjacencies considered

    block_stats = {}

    for ps, keys in blocks.items():
        keys_sorted = sorted(keys, key=lambda x: x[1])  # sort by POS

        # For biallelic het: represent phase by hap1 allele (0 or 1)
        truth_h1 = []
        pred_h1 = []

        for k in keys_sorted:
            ta, tb = truth[k]["a"], truth[k]["b"]
            pa, pb = pred[k]["a"], pred[k]["b"]
            # truth must be phased for meaningful switch calc; assume it is
            # hap1 allele = left allele in GT string when phased
            # but we already parsed (a,b) in order as appears in GT
            truth_h1.append(int(ta))   # ta is first allele in GT string
            pred_h1.append(int(pa))    # pa is first allele in GT string

        # Flip option: invert pred_h1 (0->1,1->0) for het biallelic
        # This is valid only because we’re in SNP-only biallelic setup.
        match_no_flip = sum(1 for t, p in zip(truth_h1, pred_h1) if t == p)
        match_flip = sum(1 for t, p in zip(truth_h1, pred_h1) if t == (1 - p))

        do_flip = match_flip > match_no_flip
        best_match = match_flip if do_flip else match_no_flip

        correct += best_match
        total_phased_het += len(keys_sorted)

        # Switch errors: compare relative phase between adjacent sites
        if len(keys_sorted) >= 2:
            adj = len(keys_sorted) - 1
            switch_den += adj
            # apply flip if chosen
            pseq = [(1 - p) for p in pred_h1] if do_flip else pred_h1
            for i in range(adj):
                rel_t = truth_h1[i] ^ truth_h1[i + 1]
                rel_p = pseq[i] ^ pseq[i + 1]
                if rel_t != rel_p:
                    switches += 1

        block_stats[ps] = {
            "size": len(keys_sorted),
            "best_flip": bool(do_flip),
            "match_no_flip": int(match_no_flip),
            "match_flip": int(match_flip),
            "best_match": int(best_match),
        }

    phase_accuracy = (correct / total_phased_het) if total_phased_het else None
    switch_error_rate = (switches / switch_den) if switch_den else None

    report = {
        "shared_total_records": shared_total,
        "shared_snp_records": shared_snp,
        "shared_het_records": shared_het,
        "pred_phased_het_records": phased_het,
        "pred_unphased_het_records": unphased_het,
        "num_phase_sets": len(blocks),
        "phase_accuracy_blockflip": phase_accuracy,
        "switches": int(switches),
        "switch_den": int(switch_den),
        "switch_error_rate": switch_error_rate,
        "phase_sets": block_stats,
        "notes": "Evaluation assumes biallelic SNPs. Indels/multiallelic not handled yet.",
    }

    print(json.dumps(report, indent=2))

    if args.out:
        with open(args.out, "w") as f:
            json.dump(report, f, indent=2)


if __name__ == "__main__":
    main()
