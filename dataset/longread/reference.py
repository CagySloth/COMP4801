import argparse
import json
from pathlib import Path
from typing import List, Tuple, Dict, Any

import numpy as np


DNA = np.array(["A", "C", "G", "T"], dtype="<U1")


def _wrap_fasta(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width)) + "\n"


def _intervals_overlap(a: Tuple[int, int], b: Tuple[int, int]) -> bool:
    # intervals are [start, end) half-open
    return a[0] < b[1] and a[1] > b[0]


def _interval_gap(a: Tuple[int, int], b: Tuple[int, int]) -> int:
    """Distance between two half-open intervals.
    Returns 0 if they overlap or touch, otherwise the number of bases between them.
    """
    if _intervals_overlap(a, b) or a[1] == b[0] or b[1] == a[0]:
        return 0
    if a[1] < b[0]:
        return b[0] - a[1]
    return a[0] - b[1]


def _pick_non_overlapping_start_with_gap(
    rng: np.random.Generator,
    genome_len: int,
    region_len: int,
    occupied: List[Tuple[int, int]],
    anchor: Tuple[int, int],
    min_gap: int,
    max_tries: int = 5000,
) -> int:
    """Pick a non-overlapping interval that is also at least min_gap away from anchor."""
    if region_len > genome_len:
        raise ValueError(f"region_len {region_len} > genome_len {genome_len}")

    for _ in range(max_tries):
        start = int(rng.integers(0, genome_len - region_len + 1))
        cand = (start, start + region_len)
        if any(_intervals_overlap(cand, x) for x in occupied):
            continue
        if _interval_gap(cand, anchor) < int(min_gap):
            continue
        occupied.append(cand)
        return start

    raise RuntimeError(
        f"Could not place a non-overlapping region of length {region_len} with min_gap={min_gap} "
        f"after {max_tries} tries. Try reducing dup counts/lengths or increasing genome length."
    )


def _pick_non_overlapping_start(
    rng: np.random.Generator,
    genome_len: int,
    region_len: int,
    occupied: List[Tuple[int, int]],
    max_tries: int = 2000,
) -> int:
    if region_len > genome_len:
        raise ValueError(f"region_len {region_len} > genome_len {genome_len}")

    for _ in range(max_tries):
        start = int(rng.integers(0, genome_len - region_len + 1))
        cand = (start, start + region_len)
        if all(not _intervals_overlap(cand, x) for x in occupied):
            occupied.append(cand)
            return start

    raise RuntimeError(
        f"Could not place a non-overlapping region of length {region_len} after {max_tries} tries. "
        f"Try reducing counts/lengths or increasing genome length."
    )


def _sample_bases_with_gc(rng: np.random.Generator, length: int, gc_prob: float) -> List[str]:
    # gc_prob = P(G or C)
    p_gc = float(gc_prob) / 2.0
    p_at = (1.0 - float(gc_prob)) / 2.0
    probs = np.array([p_at, p_gc, p_gc, p_at], dtype=float)  # A,C,G,T
    bases = rng.choice(DNA, size=length, replace=True, p=probs)
    return bases.tolist()


def generate_reference(
    length: int,
    seed: int | None,
    contig: str,
    preset: str,
    homopolymer_count: int,
    homopolymer_min_len: int,
    homopolymer_max_len: int,
    str_count: int,
    str_min_motif: int,
    str_max_motif: int,
    str_min_repeats: int,
    str_max_repeats: int,
    gc_window_count: int,
    gc_window_len: int,
    gc_high: float,
    gc_low: float,
    dup_segments: int = 0,
    dup_len: int = 3000,
    dup_min_gap: int = 500,
) -> Tuple[str, List[Dict[str, Any]], List[Dict[str, Any]]]:
    
    rng = np.random.default_rng(seed)
    duplications: List[Dict[str, Any]] = []

    # Base reference: uniform random DNA
    seq = rng.choice(DNA, size=length, replace=True).tolist()

    regions: List[Dict[str, Any]] = []
    occupied: List[Tuple[int, int]] = []

    # Optional presets (auto-tune counts for "realistic" vs "plain")
    if preset == "plain":
        homopolymer_count = 0
        str_count = 0
        gc_window_count = 0
    elif preset == "realistic":
        # Scale counts with length; keep small enough to fit even for moderate lengths.
        # These are starting defaults you can refine later.
        homopolymer_count = max(homopolymer_count, length // 5000)
        str_count = max(str_count, length // 5000)
        gc_window_count = max(gc_window_count, length // 20000)
    elif preset == "toy":
        # small but visible complexity even for short refs
        homopolymer_count = max(homopolymer_count, 3)
        str_count = max(str_count, 3)
        gc_window_count = max(gc_window_count, 1)
    else:
        raise ValueError(f"Unknown preset: {preset}")

    # --- Homopolymers ---
    for _ in range(int(homopolymer_count)):
        L = int(rng.integers(homopolymer_min_len, homopolymer_max_len + 1))
        start = _pick_non_overlapping_start(rng, length, L, occupied)
        base = str(rng.choice(DNA))
        seq[start : start + L] = [base] * L
        regions.append(
            {
                "type": "homopolymer",
                "start0": start,
                "end0": start + L,
                "base": base,
                "length": L,
            }
        )

    # --- STRs (short tandem repeats) ---
    for _ in range(int(str_count)):
        motif_len = int(rng.integers(str_min_motif, str_max_motif + 1))
        repeats = int(rng.integers(str_min_repeats, str_max_repeats + 1))
        motif = "".join(rng.choice(DNA, size=motif_len, replace=True).tolist())
        pattern = motif * repeats
        L = len(pattern)
        if L <= 0 or L > length:
            continue
        start = _pick_non_overlapping_start(rng, length, L, occupied)
        seq[start : start + L] = list(pattern)
        regions.append(
            {
                "type": "str",
                "start0": start,
                "end0": start + L,
                "motif": motif,
                "repeats": repeats,
                "length": L,
            }
        )

    # --- GC windows (high/low) ---
    for i in range(int(gc_window_count)):
        L = int(gc_window_len)
        if L <= 0 or L > length:
            continue
        start = _pick_non_overlapping_start(rng, length, L, occupied)

        if i % 2 == 0:
            gc = float(gc_high)
            label = "gc_high"
        else:
            gc = float(gc_low)
            label = "gc_low"

        window = _sample_bases_with_gc(rng, L, gc_prob=gc)
        seq[start : start + L] = window
        regions.append(
            {
                "type": label,
                "start0": start,
                "end0": start + L,
                "gc_target": gc,
                "length": L,
            }
        )
        
    # --- Duplicated segments (repeat-like regions) ---
    # These create mapping ambiguity, which is a key real-world challenge.
    # We store them separately from "regions" so that --avoid-regions in truth.py
    # can still mean "avoid complex regions" unless the user explicitly wants
    # to avoid duplications.
    dup_segments = int(dup_segments)
    dup_len = int(dup_len)
    dup_min_gap = int(dup_min_gap)
    if dup_segments > 0:
        if dup_len <= 0:
            raise ValueError("dup_len must be > 0")
        if dup_len * 2 > length:
            raise ValueError("dup_len too large for reference length")
        # Quick feasibility check (approximate)
        required = 2 * dup_segments * dup_len
        if required > int(0.9 * length):
            raise ValueError(
                f"duplication request too large for genome: need ~{required}bp of non-overlapping space "
                f"but length={length}. Reduce --dup-segments/--dup-len or increase --length."
            )

        occupied_dups: List[Tuple[int, int]] = []
        for _ in range(dup_segments):
            src_start = _pick_non_overlapping_start(rng, length, dup_len, occupied_dups)
            src_iv = (src_start, src_start + dup_len)
            dst_start = _pick_non_overlapping_start_with_gap(
                rng,
                genome_len=length,
                region_len=dup_len,
                occupied=occupied_dups,
                anchor=src_iv,
                min_gap=dup_min_gap,
            )

            # Copy the sequence (including any complex patterns already inserted).
            seq[dst_start : dst_start + dup_len] = seq[src_start : src_start + dup_len]

            duplications.append(
                {
                    "type": "duplication",
                    "src_start0": src_start,
                    "src_end0": src_start + dup_len,
                    "dst_start0": dst_start,
                    "dst_end0": dst_start + dup_len,
                    "length": dup_len,
                    "min_gap": dup_min_gap,
                }
            )

    ref_seq = "".join(seq)
    return ref_seq, regions, duplications

def write_reference_fasta(path: Path, contig: str, seq: str) -> None:
    with open(path, "w") as f:
        f.write(f">{contig}\n")
        f.write(_wrap_fasta(seq, width=60))


def main(args=None) -> None:
    parser = argparse.ArgumentParser(
        description="Generate a synthetic reference FASTA with optional complex regions (homopolymers/STR/GC windows)."
    )
    parser.add_argument("-o", "--output-prefix", required=True, help="Prefix for outputs (e.g. output/demo)")
    parser.add_argument("--length", type=int, default=100000, help="Reference length in bp")
    parser.add_argument("--contig", type=str, default="chr1", help="Contig name for FASTA header")
    parser.add_argument("--seed", type=int, default=None, help="Random seed")

    parser.add_argument(
        "--preset",
        choices=["plain", "toy", "realistic"],
        default="plain",
        help="plain=no complex regions; toy=small complexity; realistic=scaled complexity.",
    )

    # Advanced knobs (defaults are conservative; preset may override counts)
    parser.add_argument("--homopolymer-count", type=int, default=0)
    parser.add_argument("--homopolymer-min-len", type=int, default=6)
    parser.add_argument("--homopolymer-max-len", type=int, default=20)

    parser.add_argument("--str-count", type=int, default=0)
    parser.add_argument("--str-min-motif", type=int, default=2)
    parser.add_argument("--str-max-motif", type=int, default=6)
    parser.add_argument("--str-min-repeats", type=int, default=3)
    parser.add_argument("--str-max-repeats", type=int, default=10)

    parser.add_argument("--gc-window-count", type=int, default=0)
    parser.add_argument("--gc-window-len", type=int, default=200)
    parser.add_argument("--gc-high", type=float, default=0.75, help="Target GC fraction for GC-high windows")
    parser.add_argument("--gc-low", type=float, default=0.25, help="Target GC fraction for GC-low windows")

    parser.add_argument("--dup-segments", type=int, default=0, help="Number of duplicated segments to create")
    parser.add_argument("--dup-len", type=int, default=3000, help="Length of each duplicated segment (bp)")
    parser.add_argument("--dup-min-gap", type=int, default=500, help="Minimum gap (bp) between the source and destination copies")

    if args is None:
        args = parser.parse_args()

    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    seq, regions, duplications = generate_reference(
        length=int(args.length),
        seed=args.seed,
        contig=str(args.contig),
        preset=str(args.preset),
        homopolymer_count=int(args.homopolymer_count),
        homopolymer_min_len=int(args.homopolymer_min_len),
        homopolymer_max_len=int(args.homopolymer_max_len),
        str_count=int(args.str_count),
        str_min_motif=int(args.str_min_motif),
        str_max_motif=int(args.str_max_motif),
        str_min_repeats=int(args.str_min_repeats),
        str_max_repeats=int(args.str_max_repeats),
        gc_window_count=int(args.gc_window_count),
        gc_window_len=int(args.gc_window_len),
        gc_high=float(args.gc_high),
        gc_low=float(args.gc_low),
        
        dup_segments=int(args.dup_segments),
        dup_len=int(args.dup_len),
        dup_min_gap=int(args.dup_min_gap),
    )

    fasta_path = Path(str(prefix) + ".ref.fasta")
    meta_path = Path(str(prefix) + ".ref.meta.json")

    write_reference_fasta(fasta_path, contig=str(args.contig), seq=seq)

    meta = {
        "contig": str(args.contig),
        "length": int(args.length),
        "seed": args.seed,
        "preset": str(args.preset),
        "regions": regions,
        "duplications": duplications,
        "coord_note": "All region coordinates are 0-based half-open [start0, end0).",
    }
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)

    print(f"✅ Reference FASTA written: {fasta_path}")
    print(f"✅ Reference metadata written: {meta_path}")


if __name__ == "__main__":
    main()
