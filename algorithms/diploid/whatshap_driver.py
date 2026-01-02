import argparse
import numpy as np

from algorithms.io.reads_data import ReadsData
from algorithms.io.writer import (
    write_haplotypes_tsv,
    write_assignments_tsv,
    write_summary_json,
)

from algorithms.diploid.whatshap_adapter import build_readset_from_readsdata
from whatshap import core, readselect
# from whatshap.blocks import compute_overall_components  # optional

def main(args=None):
    parser = argparse.ArgumentParser(description="Diploid phasing via WhatsHap core")
    parser.add_argument("-i", "--input", required=True, help="Input NPZ file")
    parser.add_argument("--output-prefix", required=True, help="Prefix for output files")
    parser.add_argument("--max-coverage", type=int, default=15)
    parser.add_argument("--error-rate", type=float, default=0.01)
    # error_rate kept for interface consistency; not used yet

    if args is None:
        args = parser.parse_args()

    # 1) Load reads
    data = ReadsData.from_npz(args.input)
    N = data.N
    R = data.R

    # 2) Build ReadSet
    readset = build_readset_from_readsdata(data)
    readset.sort()

    # 3) Perform WhatsHap read selection
    selected_indices = readselect.readselection(readset, args.max_coverage, None)
    selected_readset = readset.subset(selected_indices)

    # 4) Run HapChatCore MEC solver (single sample, no pedigree)
    hap_core = core.HapChatCore(selected_readset)
    superreads_list, _ = hap_core.get_super_reads()

    # 5) Extract haplotypes from superreads
    hap1 = np.full(N, -1, dtype=int)
    hap2 = np.full(N, -1, dtype=int)

    # Each element of superreads_list is a ReadSet for one block, with two reads: hap0 and hap1
    for block in superreads_list:
        if len(block) < 2:
            continue
        hap_read0 = block[0]
        hap_read1 = block[1]
        for variant in hap_read0:
            hap1[variant.position] = variant.allele
        for variant in hap_read1:
            hap2[variant.position] = variant.allele

    # 6) Fallback: fill any sites HapChat left as -1 using the observed alleles per column
    for pos in range(N):
        if hap1[pos] == -1 and hap2[pos] == -1:
            col = data.alleles[:, pos]
            alleles = sorted({int(a) for a in col if a >= 0})
            if len(alleles) == 1:
                hap1[pos] = hap2[pos] = alleles[0]
            elif len(alleles) >= 2:
                hap1[pos], hap2[pos] = alleles[0], alleles[1]
            else:
                # no information at this position; arbitrarily choose 0/0
                hap1[pos] = hap2[pos] = 0

    H = np.stack([hap1, hap2], axis=0)

    # 7) Placeholder assignments: all selected reads assigned to haplotype 0
    assignments = np.zeros(len(selected_indices), dtype=np.int32)

    # 8) Write outputs in the same format as other algorithms
    prefix = args.output_prefix
    write_haplotypes_tsv(f"{prefix}.haplotypes.tsv", H)
    write_assignments_tsv(f"{prefix}.assignments.tsv", assignments)
    write_summary_json(
        {
            "algorithm": "diploid_whats",
            "R": int(R),
            "N": int(N),
        },
        f"{prefix}.summary.json",
    )


if __name__ == "__main__":
    main()


# import argparse
# import numpy as np

# from algorithms.io.reads_data import ReadsData
# from algorithms.io.writer import (
#     write_haplotypes_tsv,
#     write_assignments_tsv,
#     write_summary_json,
# )

# from algorithms.diploid.whatshap_adapter import build_readset_from_readsdata
# from whatshap import core, readselect
# # from whatshap.blocks import compute_overall_components

# def main(args=None):
#     parser = argparse.ArgumentParser(description="Diploid phasing via WhatsHap core")
#     parser.add_argument("-i", "--input", required=True, help="Input NPZ file")
#     parser.add_argument("--output-prefix", required=True, help="Prefix for output files")
#     parser.add_argument("--max-coverage", type=int, default=15)
#     parser.add_argument("--error-rate", type=float, default=0.01)
#     # ... any other parameters (e.g. distrust-genotypes flag) you want to add later

#     if args is None:
#         args = parser.parse_args()

#     # 1) Load reads
#     data = ReadsData.from_npz(args.input)
#     N = data.N
#     R = data.R

#     # 2) Build ReadSet
#     readset = build_readset_from_readsdata(data)

#     # 3) Perform WhatsHap read selection
#     selected_indices = readselect.readselection(readset, args.max_coverage, None)
#     selected_readset = readset.subset(selected_indices)

#     # 4) Run PedigreeDPTable (MEC) on the whole chromosome
#     #    For now, no pedigree -> recombi costs empty, distrust_genotypes=False, positions=None
    
#     # hap_core = core.HapChatCore(selected_readset)
#     # superreads_list, _ = hap_core.get_super_reads()
    
#     recombcosts = []  # no recombination cost for simple diploid single-sample
#     pedigree = core.Pedigree()  # empty pedigree; check constructor in core.pyx
#     dp = core.PedigreeDPTable(selected_readset, recombcosts, pedigree, False, None)

#     # This actually runs the DP and returns superreads (consensus haplotypes)
#     superreads_list, transmission = dp.get_super_reads()
#     # For a single sample, superreads_list[0] is a ReadSet with two reads (hap1 & hap2)

#     # 5) Extract haplotypes from superreads
#     hap1 = np.full(N, -1, dtype=int)
#     hap2 = np.full(N, -1, dtype=int)

#     # The positions in superreads correspond to variant positions in the original coordinate system
#     superreads = superreads_list[0]
#     # superreads[0] = hap1-read, superreads[1] = hap2-read (by convention)
#     # iterate through variants in each superread and fill hap1/hap2 arrays
#     hap_read0 = superreads[0]
#     hap_read1 = superreads[1]

#     for variant in hap_read0:
#         hap1[variant.position] = variant.allele
#     for variant in hap_read1:
#         hap2[variant.position] = variant.allele

#     # 6) Optionally compute components/blocks via compute_overall_components
#     #    (use accessible_positions = [0..N-1], homozygous_positions=[], etc.)

#     # 7) Assign each selected read to hap1/hap2 (you can use read's 'phase' info if present,
#     #    or compute your own by comparing to H)

#     H = np.stack([hap1, hap2], axis=0)
#     assignments = np.zeros(len(selected_indices), dtype=np.int32)  # placeholder for now

#     # 8) Write outputs in the same format as other algorithms
#     prefix = args.output_prefix
#     write_haplotypes_tsv(f"{prefix}.haplotypes.tsv", H)
#     write_assignments_tsv(f"{prefix}.assignments.tsv", assignments)
#     write_summary_json(
#         {
#             "algorithm": "diploid_whats",
#             "R": int(R),
#             "N": int(N),
#             # Add MEC, number of phased variants, etc. once you compute them
#         },
#         f"{prefix}.summary.json",
#     )
