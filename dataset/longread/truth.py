# dataset/longread/truth.py
import argparse
import json
from pathlib import Path
from bisect import bisect_left, insort
from typing import List, Tuple, Dict, Any

import numpy as np

DNA = ["A", "C", "G", "T"]


def read_single_fasta(path: Path) -> Tuple[str, str]:
    """
    Returns (contig_name, sequence).
    Assumes single-contig FASTA.
    """
    name = None
    seq_parts = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    raise ValueError("FASTA has more than one contig; only single-contig supported for now.")
                name = line[1:].split()[0]
            else:
                seq_parts.append(line.upper())
    if name is None:
        raise ValueError("FASTA missing header line.")
    seq = "".join(seq_parts)
    return name, seq


def write_fasta(path: Path, name: str, seq: str, width: int = 60) -> None:
    with open(path, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), width):
            f.write(seq[i : i + width] + "\n")


def load_excluded_intervals(meta_path: Path) -> List[Tuple[int, int]]:
    """
    Reads Step-1 meta.json and returns list of [start0, end0) intervals to exclude.
    """
    with open(meta_path, "r") as f:
        meta = json.load(f)
    regions = meta.get("regions", [])
    intervals = []
    for r in regions:
        s = int(r["start0"])
        e = int(r["end0"])
        intervals.append((s, e))
    return intervals


def is_in_any_interval(pos: int, intervals: List[Tuple[int, int]]) -> bool:
    for s, e in intervals:
        if s <= pos < e:
            return True
    return False


def pick_positions(
    rng: np.random.Generator,
    seq: str,
    num_snps: int,
    min_distance: int,
    excluded_intervals: List[Tuple[int, int]] | None,
) -> List[int]:
    """
    Pick SNP positions (0-based) such that:
    - base at pos is A/C/G/T
    - not inside excluded intervals (if provided)
    - min_distance between picked positions
    """
    n = len(seq)
    candidates = []
    for i, b in enumerate(seq):
        if b in ("A", "C", "G", "T"):
            if excluded_intervals and is_in_any_interval(i, excluded_intervals):
                continue
            candidates.append(i)

    if len(candidates) < num_snps:
        raise ValueError(f"Not enough candidate positions to place {num_snps} SNPs (only {len(candidates)} available).")

    rng.shuffle(candidates)

    chosen_sorted: List[int] = []
    for pos in candidates:
        idx = bisect_left(chosen_sorted, pos)
        if idx > 0 and pos - chosen_sorted[idx - 1] <= min_distance:
            continue
        if idx < len(chosen_sorted) and chosen_sorted[idx] - pos <= min_distance:
            continue
        insort(chosen_sorted, pos)
        if len(chosen_sorted) == num_snps:
            break

    if len(chosen_sorted) < num_snps:
        raise RuntimeError(
            f"Could only place {len(chosen_sorted)}/{num_snps} SNPs with min_distance={min_distance}. "
            f"Try reducing num_snps or min_distance, or disable avoid-regions."
        )

    return chosen_sorted


def choose_alt(rng: np.random.Generator, ref: str) -> str:
    alts = [b for b in DNA if b != ref]
    return str(rng.choice(alts))


def main(args=None) -> None:
    parser = argparse.ArgumentParser(
        description="Generate diploid truth (truth VCF + hap1/hap2 FASTA) from a reference FASTA."
    )
    parser.add_argument("--ref", required=True, help="Input reference FASTA (from Step 1), e.g. output/demo.ref.fasta")
    parser.add_argument("-o", "--output-prefix", required=True, help="Output prefix, e.g. output/demo")

    parser.add_argument("--seed", type=int, default=None, help="Random seed")
    parser.add_argument("--sample", type=str, default="SAMPLE", help="Sample name in VCF")

    # SNP-only for Step 2
    parser.add_argument("--num-snps", type=int, default=1000, help="Number of SNPs to inject")
    parser.add_argument("--het-rate", type=float, default=0.8, help="Fraction of injected SNPs that are heterozygous")
    parser.add_argument("--min-distance", type=int, default=5, help="Minimum distance between SNP sites (bp)")

    # Optional: avoid placing variants inside complex regions from reference meta.json
    parser.add_argument("--ref-meta", help="Reference meta JSON from Step 1 (e.g. output/demo.ref.meta.json)")
    parser.add_argument("--avoid-regions", action="store_true", help="Avoid placing SNPs inside meta regions")

    # Truth phasing style
    parser.add_argument(
        "--phased-truth",
        action="store_true",
        help="Write phased genotypes (0|1 / 1|0) in truth.vcf (recommended for evaluation).",
    )
    parser.add_argument(
        "--random-phase",
        action="store_true",
        help="For het SNPs, randomly choose whether hap1 is REF or ALT (otherwise hap1=REF, hap2=ALT).",
    )

    if args is None:
        args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    ref_path = Path(args.ref)
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    contig, ref_seq = read_single_fasta(ref_path)
    L = len(ref_seq)

    excluded = None
    if args.avoid_regions:
        if not args.ref_meta:
            raise ValueError("--avoid-regions requires --ref-meta")
        excluded = load_excluded_intervals(Path(args.ref_meta))

    positions = pick_positions(
        rng=rng,
        seq=ref_seq,
        num_snps=int(args.num_snps),
        min_distance=int(args.min_distance),
        excluded_intervals=excluded,
    )

    # Initialize haplotypes as reference
    hap1 = list(ref_seq)
    hap2 = list(ref_seq)

    records: List[Dict[str, Any]] = []
    het_count = 0
    hom_alt_count = 0

    for pos0 in positions:
        ref_base = ref_seq[pos0]
        alt_base = choose_alt(rng, ref_base)

        is_het = (rng.random() < float(args.het_rate))
        if is_het:
            het_count += 1
            # default: hap1=REF, hap2=ALT
            if args.random_phase and (rng.random() < 0.5):
                hap1[pos0] = alt_base
                hap2[pos0] = ref_base
                gt = "1|0" if args.phased_truth else "0/1"
            else:
                hap1[pos0] = ref_base
                hap2[pos0] = alt_base
                gt = "0|1" if args.phased_truth else "0/1"
        else:
            hom_alt_count += 1
            hap1[pos0] = alt_base
            hap2[pos0] = alt_base
            gt = "1|1" if args.phased_truth else "1/1"

        records.append(
            {
                "chrom": contig,
                "pos1": pos0 + 1,  # VCF is 1-based
                "ref": ref_base,
                "alt": alt_base,
                "gt": gt,
                "pos0": pos0,
            }
        )

    # Write haplotype FASTAs
    hap1_path = Path(str(prefix) + ".hap1.fasta")
    hap2_path = Path(str(prefix) + ".hap2.fasta")
    write_fasta(hap1_path, name=f"{contig}_hap1", seq="".join(hap1))
    write_fasta(hap2_path, name=f"{contig}_hap2", seq="".join(hap2))

    # Write truth VCF
    vcf_path = Path(str(prefix) + ".truth.vcf")
    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##contig=<ID={contig},length={L}>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + str(args.sample) + "\n")
        for r in records:
            f.write(
                f"{r['chrom']}\t{r['pos1']}\t.\t{r['ref']}\t{r['alt']}\t.\tPASS\t.\tGT\t{r['gt']}\n"
            )

    # Optional meta for debugging/benchmark reports
    meta_path = Path(str(prefix) + ".truth.meta.json")
    meta = {
        "ref": str(ref_path),
        "contig": contig,
        "length": L,
        "seed": args.seed,
        "num_snps": int(args.num_snps),
        "het_rate": float(args.het_rate),
        "het_count": int(het_count),
        "hom_alt_count": int(hom_alt_count),
        "min_distance": int(args.min_distance),
        "avoid_regions": bool(args.avoid_regions),
        "ref_meta": str(args.ref_meta) if args.ref_meta else None,
        "coord_note": "records are 0-based internally but VCF pos is 1-based.",
    }
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)

    print(f"✅ Truth VCF written: {vcf_path}")
    print(f"✅ Haplotype FASTAs written: {hap1_path}, {hap2_path}")
    print(f"✅ Truth meta written: {meta_path}")


if __name__ == "__main__":
    main()
