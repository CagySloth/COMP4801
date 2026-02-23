# dataset/longread/readsim.py
import argparse
from pathlib import Path
import numpy as np

_BASES = np.array(list("ACGT"), dtype="<U1")


def read_single_fasta(path: Path) -> tuple[str, str]:
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
        raise ValueError("FASTA missing header.")
    return name, "".join(seq_parts)


def write_fastq_record(f, name: str, seq: str, qual: str) -> None:
    f.write(f"@{name}\n")
    f.write(seq + "\n")
    f.write("+\n")
    f.write(qual + "\n")


def _phred_char(q: float) -> str:
    # Clamp to a sensible FASTQ range. (Phred+33)
    qi = int(round(q))
    qi = 0 if qi < 0 else qi
    qi = 41 if qi > 41 else qi
    return chr(qi + 33)


def _sample_q(frac: float, rng: np.random.Generator, q_mean: float, q_std: float, q_end_drop: float) -> float:
    # Simple “worse at ends, better in middle” profile
    # frac in [0,1]
    end_penalty = q_end_drop * (abs(frac - 0.5) / 0.5)  # 0 at center, q_end_drop at ends
    mu = q_mean - end_penalty
    return float(rng.normal(mu, q_std))


def _rand_base(rng: np.random.Generator) -> str:
    return str(rng.choice(_BASES))


def _mutate_base(base: str, rng: np.random.Generator) -> str:
    # pick a different base
    choices = [b for b in "ACGT" if b != base]
    return str(rng.choice(choices))


def simulate_ont_read(
    seq: str,
    rng: np.random.Generator,
    *,
    sub_rate: float,
    ins_rate: float,
    del_rate: float,
    homopolymer_len: int,
    homopolymer_factor: float,
    max_ins_per_base: int,
    q_mean: float,
    q_std: float,
    q_end_drop: float,
    q_sub_penalty: float,
    q_indel_penalty: float,
) -> tuple[str, str, int, int, int]:
    """
    A simple ONT-like error model:
      - substitutions + insertions + deletions
      - indels increased in homopolymers (systematic bias)
      - qualities vary across read; lower qualities for error events

    Returns: (mut_seq, qual_str, n_sub, n_ins, n_del)
    """
    out_bases: list[str] = []
    out_quals: list[str] = []

    n_sub = 0
    n_ins = 0
    n_del = 0

    L = len(seq)
    if L == 0:
        return "", "", 0, 0, 0

    i = 0
    while i < L:
        base = seq[i]

        # homopolymer run length starting at i
        run = 1
        while i + run < L and seq[i + run] == base:
            run += 1
        hp = run >= homopolymer_len

        ins_p = ins_rate * (homopolymer_factor if hp else 1.0)
        del_p = del_rate * (homopolymer_factor if hp else 1.0)

        frac = (i / (L - 1)) if L > 1 else 0.5

        # 1) insertions BEFORE this base (geometric-ish with cap)
        k = 0
        p = ins_p
        while k < max_ins_per_base and rng.random() < p:
            b_ins = _rand_base(rng)
            q = _sample_q(frac, rng, q_mean, q_std, q_end_drop) - q_indel_penalty
            out_bases.append(b_ins)
            out_quals.append(_phred_char(q))
            n_ins += 1
            k += 1
            # reduce probability of multiple consecutive inserts
            p *= 0.5

        # 2) deletion of this base
        if rng.random() < del_p:
            n_del += 1
            i += 1
            continue

        # 3) substitution (or match)
        q = _sample_q(frac, rng, q_mean, q_std, q_end_drop)
        if rng.random() < sub_rate:
            base2 = _mutate_base(base, rng)
            q -= q_sub_penalty
            n_sub += 1
        else:
            base2 = base

        out_bases.append(base2)
        out_quals.append(_phred_char(q))
        i += 1

    return "".join(out_bases), "".join(out_quals), n_sub, n_ins, n_del


def main(args=None) -> None:
    parser = argparse.ArgumentParser(description="Simulate long reads from hap1/hap2 FASTA (perfect or ONT-like).")
    parser.add_argument("--hap1", required=True, help="Haplotype1 FASTA (from Step 2)")
    parser.add_argument("--hap2", required=True, help="Haplotype2 FASTA (from Step 2)")
    parser.add_argument("-o", "--output-prefix", required=True, help="Output prefix, e.g. output/demo")
    parser.add_argument("--num-reads", type=int, default=200, help="Number of reads to simulate")
    parser.add_argument("--seed", type=int, default=None)

    parser.add_argument("--platform", choices=["perfect", "ont"], default="ont", help="Read error/quality model")

    # Simple length model (kept)
    parser.add_argument("--min-len", type=int, default=5000)
    parser.add_argument("--max-len", type=int, default=15000)

    # Optional: haplotype mixture
    parser.add_argument("--hap1-frac", type=float, default=0.5, help="Fraction of reads drawn from hap1")

    # Perfect mode: constant quality
    parser.add_argument("--qual-char", type=str, default="I", help="(perfect mode) single FASTQ quality char")

    # ONT presets (you can override with explicit rates below)
    parser.add_argument("--ont-profile", choices=["classic", "q20"], default="classic")

    # ONT error rates per reference base (approx). You can tune.
    parser.add_argument("--sub-rate", type=float, default=None)
    parser.add_argument("--ins-rate", type=float, default=None)
    parser.add_argument("--del-rate", type=float, default=None)

    # ONT: homopolymer bias (indels more likely)
    parser.add_argument("--homopolymer-len", type=int, default=4)
    parser.add_argument("--homopolymer-factor", type=float, default=2.0)
    parser.add_argument("--max-ins-per-base", type=int, default=3)

    # ONT: quality model (Phred)
    parser.add_argument("--q-mean", type=float, default=None)
    parser.add_argument("--q-std", type=float, default=2.5)
    parser.add_argument("--q-end-drop", type=float, default=4.0)
    parser.add_argument("--q-sub-penalty", type=float, default=3.0)
    parser.add_argument("--q-indel-penalty", type=float, default=4.0)

    if args is None:
        args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    contig1, seq1 = read_single_fasta(Path(args.hap1))
    contig2, seq2 = read_single_fasta(Path(args.hap2))
    if len(seq1) != len(seq2):
        raise ValueError("hap1 and hap2 must have same length")

    L = len(seq1)
    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    out_fastq = Path(str(prefix) + ".reads.fastq")
    out_truth = Path(str(prefix) + ".reads.truth.tsv")

    hap1_frac = float(args.hap1_frac)
    if not (0.0 <= hap1_frac <= 1.0):
        raise ValueError("--hap1-frac must be between 0 and 1")

    min_len = int(args.min_len)
    max_len = int(args.max_len)
    if min_len < 1 or max_len < min_len:
        raise ValueError("invalid min/max len")

    if args.platform == "perfect":
        qual_char = str(args.qual_char)
        if len(qual_char) != 1:
            raise ValueError("--qual-char must be a single character like 'I'")

    # ONT defaults:
    # A widely reported ONT-like split is on the order of a few % each for sub/ins/del,
    # with indels prominent and systematic errors in homopolymers. :contentReference[oaicite:1]{index=1}
    if args.ont_profile == "classic":
        sub_rate = 0.025 if args.sub_rate is None else float(args.sub_rate)
        ins_rate = 0.022 if args.ins_rate is None else float(args.ins_rate)
        del_rate = 0.030 if args.del_rate is None else float(args.del_rate)
        q_mean = 12.0 if args.q_mean is None else float(args.q_mean)
    else:  # q20
        # “Q20 chemistry” aims for >99% single-read accuracy. :contentReference[oaicite:2]{index=2}
        sub_rate = 0.004 if args.sub_rate is None else float(args.sub_rate)
        ins_rate = 0.003 if args.ins_rate is None else float(args.ins_rate)
        del_rate = 0.003 if args.del_rate is None else float(args.del_rate)
        q_mean = 20.0 if args.q_mean is None else float(args.q_mean)

    for r in (sub_rate, ins_rate, del_rate):
        if r < 0 or r > 0.5:
            raise ValueError("sub/ins/del rates must be in a sensible range (0..0.5)")

    with open(out_fastq, "w") as fq, open(out_truth, "w") as tf:
        tf.write("read_id\thap\tstart0\tlen_ref\tlen_obs\tsubs\tins\tdels\tcontig\n")

        for i in range(int(args.num_reads)):
            hap = 1 if (rng.random() < hap1_frac) else 2
            hap_seq = seq1 if hap == 1 else seq2

            rlen = int(rng.integers(min_len, max_len + 1))
            if rlen > L:
                rlen = L

            start0 = int(rng.integers(0, L - rlen + 1))
            read_seq_ref = hap_seq[start0 : start0 + rlen]

            read_id = f"r{i}_hap{hap}_{contig1}:{start0}-{start0+rlen}"

            if args.platform == "perfect":
                read_seq = read_seq_ref
                qual = str(args.qual_char) * len(read_seq)
                subs = ins = dels = 0
            else:
                read_seq, qual, subs, ins, dels = simulate_ont_read(
                    read_seq_ref,
                    rng,
                    sub_rate=sub_rate,
                    ins_rate=ins_rate,
                    del_rate=del_rate,
                    homopolymer_len=int(args.homopolymer_len),
                    homopolymer_factor=float(args.homopolymer_factor),
                    max_ins_per_base=int(args.max_ins_per_base),
                    q_mean=q_mean,
                    q_std=float(args.q_std),
                    q_end_drop=float(args.q_end_drop),
                    q_sub_penalty=float(args.q_sub_penalty),
                    q_indel_penalty=float(args.q_indel_penalty),
                )

            write_fastq_record(fq, read_id, read_seq, qual)
            tf.write(f"{read_id}\t{hap}\t{start0}\t{rlen}\t{len(read_seq)}\t{subs}\t{ins}\t{dels}\t{contig1}\n")

    print(f"✅ FASTQ written: {out_fastq}")
    print(f"✅ Read truth TSV written: {out_truth}")
    if args.platform == "ont":
        total_err = sub_rate + ins_rate + del_rate
        print(f"ℹ️  ONT profile={args.ont_profile} sub/ins/del={sub_rate:.4f}/{ins_rate:.4f}/{del_rate:.4f} (total≈{total_err:.4f})")


if __name__ == "__main__":
    main()
