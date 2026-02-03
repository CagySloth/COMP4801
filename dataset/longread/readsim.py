# dataset/longread/readsim.py
import argparse
from pathlib import Path
import numpy as np


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


def main(args=None) -> None:
    parser = argparse.ArgumentParser(description="Simulate perfect (error-free) long reads from hap1/hap2 FASTA.")
    parser.add_argument("--hap1", required=True, help="Haplotype1 FASTA (from Step 2)")
    parser.add_argument("--hap2", required=True, help="Haplotype2 FASTA (from Step 2)")
    parser.add_argument("-o", "--output-prefix", required=True, help="Output prefix, e.g. output/demo")
    parser.add_argument("--num-reads", type=int, default=200, help="Number of reads to simulate")
    parser.add_argument("--seed", type=int, default=None)

    # Simple length model for Step 3A
    parser.add_argument("--min-len", type=int, default=5000)
    parser.add_argument("--max-len", type=int, default=15000)

    # Optional: haplotype mixture
    parser.add_argument("--hap1-frac", type=float, default=0.5, help="Fraction of reads drawn from hap1")

    # Quality (constant for now)
    parser.add_argument("--qual-char", type=str, default="I", help="Single FASTQ quality character to repeat")

    if args is None:
        args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    contig1, seq1 = read_single_fasta(Path(args.hap1))
    contig2, seq2 = read_single_fasta(Path(args.hap2))
    if len(seq1) != len(seq2):
        raise ValueError("hap1 and hap2 must have same length")
    if contig1 != contig2 and contig2:
        # This is not fatal, but usually both come from same reference contig
        pass

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

    qual_char = str(args.qual_char)
    if len(qual_char) != 1:
        raise ValueError("--qual-char must be a single character like 'I'")

    with open(out_fastq, "w") as fq, open(out_truth, "w") as tf:
        tf.write("read_id\thap\tstart0\tlen\tcontig\n")

        for i in range(int(args.num_reads)):
            # choose hap
            hap = 1 if (rng.random() < hap1_frac) else 2
            hap_seq = seq1 if hap == 1 else seq2

            # choose read length uniformly for now
            rlen = int(rng.integers(min_len, max_len + 1))
            if rlen > L:
                rlen = L

            # choose start
            start0 = int(rng.integers(0, L - rlen + 1))
            read_seq = hap_seq[start0 : start0 + rlen]

            read_id = f"r{i}_hap{hap}_{contig1}:{start0}-{start0+rlen}"

            qual = qual_char * len(read_seq)
            write_fastq_record(fq, read_id, read_seq, qual)

            tf.write(f"{read_id}\t{hap}\t{start0}\t{rlen}\t{contig1}\n")

    print(f"✅ FASTQ written: {out_fastq}")
    print(f"✅ Read truth TSV written: {out_truth}")


if __name__ == "__main__":
    main()
