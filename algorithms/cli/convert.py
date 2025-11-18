import argparse
from algorithms.io import parse_sparse_tsv, ReadsData

def main():
    parser = argparse.ArgumentParser(description="Convert between sparse TSV and NPZ formats.")
    parser.add_argument("input", help="Input .tsv or .npz file")
    parser.add_argument("--output", required=True, help="Output filename (.npz or .tsv)")
    args = parser.parse_args()

    if args.input.endswith(".tsv") and args.output.endswith(".npz"):
        fragments = parse_sparse_tsv(args.input)
        reads_data = ReadsData.from_fragments(fragments)
        reads_data.to_npz(args.output)
        print(f"Converted {args.input} → {args.output}")
    elif args.input.endswith(".npz") and args.output.endswith(".tsv"):
        reads = ReadsData.from_npz(args.input)
        reads.to_sparse_tsv(args.output)
        print(f"Converted {args.input} → {args.output}")
    else:
        raise ValueError("File extensions must match conversion direction (.tsv → .npz or .npz → .tsv)")

if __name__ == "__main__":
    main()
