import os
import argparse
import numpy as np
import whatshap as wh

from algorithms.io.reads_data import ReadsData
from algorithms.io.writer import (
    write_haplotypes_tsv,
    write_assignments_tsv,
    write_summary_json,
)

from algorithms.diploid.whatshap_adapter import build_readset_from_readsdata
from whatshap import core, readselect
# from whatshap.blocks import compute_overall_components  # optional

def _read_vcf_gt_list(vcf_path: str, sample: str | None = None):
    """
    Return:
      gt_list: list[str] GT field as in VCF (e.g. '0/1', '1/1', './.')
      alleles: list[tuple[int|None,int|None]] parsed alleles (None if missing)
      is_het: list[bool]
      sample_col: int column index used
    """
    gt_list = []
    alleles = []
    is_het = []
    sample_col = None

    with open(vcf_path, "r") as f:
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
                        raise ValueError(f"Sample '{sample}' not found in VCF. Available: {samples}")
                    sample_col = header.index(sample)
                continue
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            fmt = fields[8].split(":")
            sval = fields[sample_col].split(":")
            d = dict(zip(fmt, sval))
            gt = d.get("GT", "./.")
            gt_list.append(gt)

            # parse GT (biallelic only)
            if gt in (".", "./.", ".|."):
                alleles.append((None, None))
                is_het.append(False)
                continue

            sep = "|" if "|" in gt else "/"
            a_str, b_str = gt.split(sep)
            if a_str == "." or b_str == ".":
                alleles.append((None, None))
                is_het.append(False)
                continue

            a = int(a_str)
            b = int(b_str)
            alleles.append((a, b))
            is_het.append(a != b)

    return gt_list, alleles, is_het


def _write_phased_vcf(in_vcf: str, out_vcf: str, phased_gt: list[str], ps: list[str], sample: str | None = None):
    """
    Write a phased VCF by taking input records and replacing GT (and adding PS).
    phased_gt[i] is GT string for record i (e.g. '0|1' or '0/1' or '1/1')
    ps[i] is PS value string or '.' for record i
    """
    saw_ps_meta = False
    sample_col = None

    with open(in_vcf, "r") as fin, open(out_vcf, "w") as fout:
        for line in fin:
            if line.startswith("##"):
                if line.startswith("##FORMAT=<ID=PS,"):
                    saw_ps_meta = True
                fout.write(line)
                continue

            if line.startswith("#CHROM"):
                if not saw_ps_meta:
                    fout.write('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">\n')
                header = line.rstrip("\n").split("\t")
                if len(header) < 10:
                    raise ValueError("VCF has no sample columns")
                samples = header[9:]
                if sample is None:
                    sample_col = 9
                else:
                    if sample not in samples:
                        raise ValueError(f"Sample '{sample}' not found in VCF. Available: {samples}")
                    sample_col = header.index(sample)
                fout.write(line)
                continue

            if line.startswith("#"):
                fout.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            rec_i = len(_write_phased_vcf._seen)  # record index
            _write_phased_vcf._seen.append(1)

            fmt = fields[8].split(":")
            sval = fields[sample_col].split(":")
            d = dict(zip(fmt, sval))

            # ensure PS exists in FORMAT
            if "PS" not in fmt:
                fmt.append("PS")
                sval.append(".")
            d = dict(zip(fmt, sval))

            d["GT"] = phased_gt[rec_i]
            d["PS"] = ps[rec_i]

            # rebuild sample value in same order as fmt
            new_sample = ":".join(d.get(k, ".") for k in fmt)
            fields[8] = ":".join(fmt)
            fields[sample_col] = new_sample
            fout.write("\t".join(fields) + "\n")


# state for indexing records while streaming
_write_phased_vcf._seen = []


def main(args=None):
    parser = argparse.ArgumentParser(description="Diploid phasing via WhatsHap core")
    parser.add_argument("-i", "--input", required=True, help="Input NPZ file")
    parser.add_argument("--output-prefix", required=True, help="Prefix for output files")
    parser.add_argument("--max-coverage", type=int, default=15)
    parser.add_argument("--error-rate", type=float, default=0.01)
    parser.add_argument("--vcf", help="Input VCF with genotypes (unphased). If set, phase only het sites and write phased VCF.")
    parser.add_argument("--sample", help="Sample name in VCF (default: first sample).")
    parser.add_argument("--output-vcf", help="Output phased VCF path (default: <output-prefix>.phased.vcf).")

    # error_rate kept for interface consistency; not used yet

    if args is None:
        args = parser.parse_args()

    # 1) Load reads
    data = ReadsData.from_npz(args.input)
    N = data.N
    R = data.R
    
    vcf_path = getattr(args, "vcf", None)
    sample = getattr(args, "sample", None)
    out_vcf = getattr(args, "output_vcf", None) or f"{args.output_prefix}.phased.vcf"

    het_positions = None
    old_to_new = None
    alleles = None
    is_het = None

    data_for_phasing = data

    if vcf_path:
        _, alleles, is_het = _read_vcf_gt_list(vcf_path, sample=sample)
        if len(is_het) != N:
            raise ValueError(f"VCF variant count ({len(is_het)}) != reads matrix N ({N}). They must match in this simulation pipeline.")

        het_positions = [i for i, h in enumerate(is_het) if h]
        old_to_new = {old: new for new, old in enumerate(het_positions)}

        # Build het-only matrix (this matches real WhatsHap: phase only heterozygous variants)
        A_het = data.alleles[:, het_positions] if len(het_positions) > 0 else np.empty((R, 0), dtype=int)
        positions_het = np.tile(np.arange(A_het.shape[1]), (R, 1))
        data_for_phasing = ReadsData(reads=A_het, positions=positions_het, num_variants=A_het.shape[1])


    # 2) Build ReadSet
    readset = build_readset_from_readsdata(data_for_phasing)
    readset.sort()
    
    # --- Filter non-informative reads (must cover >= 2 variants) ---
    # WhatsHap readselection expects reads with at least two variants.
    informative_idx = [i for i, r in enumerate(readset) if len(r) >= 2]
    informative_readset = readset.subset(informative_idx)
    informative_readset.sort()

    # If there aren't enough variants or informative reads, phasing is impossible.
    N_phase = data_for_phasing.N
    if N_phase < 2 or len(informative_readset) == 0:
        selected_indices = []
        superreads_list = []
    else:
        # 3) Perform WhatsHap read selection on informative reads only
        sel_local = readselect.readselection(informative_readset, args.max_coverage, None)
        selected_readset = informative_readset.subset(sel_local)

        # If you want selected_indices relative to the ORIGINAL readset:
        selected_indices = [informative_idx[i] for i in sel_local]

        # 4) Run HapChatCore MEC solver
        hap_core = core.HapChatCore(selected_readset)
        superreads_list, _ = hap_core.get_super_reads()

    # 5) Extract haplotypes from superreads
    N_phase = data_for_phasing.N
    hap1_phase = np.full(N_phase, -1, dtype=int)
    hap2_phase = np.full(N_phase, -1, dtype=int)
    ps_phase = np.full(N_phase, -1, dtype=int)

    for block_id, block in enumerate(superreads_list, start=1):
        if len(block) < 2:
            continue
        hap_read0 = block[0]
        hap_read1 = block[1]
        for variant in hap_read0:
            hap1_phase[variant.position] = variant.allele
            ps_phase[variant.position] = block_id
        for variant in hap_read1:
            hap2_phase[variant.position] = variant.allele
            ps_phase[variant.position] = block_id

    # 6) Fallback in PHASING SPACE (het-only if VCF provided, else full space)
    A_phase = data_for_phasing.alleles  # shape: R x N_phase

    for pos in range(N_phase):
        if hap1_phase[pos] == -1 and hap2_phase[pos] == -1:
            col = A_phase[:, pos]
            obs = sorted({int(a) for a in col if a >= 0})
            if len(obs) == 1:
                hap1_phase[pos] = hap2_phase[pos] = obs[0]
            elif len(obs) >= 2:
                hap1_phase[pos], hap2_phase[pos] = obs[0], obs[1]
            else:
                hap1_phase[pos] = hap2_phase[pos] = 0
        elif hap1_phase[pos] == -1:
            hap1_phase[pos] = 1 - hap2_phase[pos]
        elif hap2_phase[pos] == -1:
            hap2_phase[pos] = 1 - hap1_phase[pos]
            
    # If VCF is provided: reconstruct FULL haplotypes (N sites) using GT, 
    # and write phased VCF (GT with | and PS for phased het sites).
    if vcf_path:
        hap1 = np.full(N, -1, dtype=int)
        hap2 = np.full(N, -1, dtype=int)
        ps_full = np.full(N, -1, dtype=int)

        for pos in range(N):
            a, b = alleles[pos]  # from VCF parsing earlier
            if a is None or b is None:
                hap1[pos] = hap2[pos] = 0
                continue

            if not is_het[pos]:
                # homozygous genotype: fixed, no phasing needed
                hap1[pos] = hap2[pos] = a
            else:
                j = old_to_new[pos]  # map full index -> het-space index
                hap1[pos] = hap1_phase[j]
                hap2[pos] = hap2_phase[j]
                ps_full[pos] = ps_phase[j]

        H = np.stack([hap1, hap2], axis=0)

        phased_gt = []
        ps_out = []
        for pos in range(N):
            a, b = alleles[pos]
            if a is None or b is None:
                phased_gt.append("./.")
                ps_out.append(".")
            elif not is_het[pos]:
                phased_gt.append(f"{a}/{b}")  # keep unphased (doesn't matter for homozygous)
                ps_out.append(".")
            else:
                # if that het site got assigned a phase set, output phased GT
                if ps_full[pos] >= 0:
                    phased_gt.append(f"{hap1[pos]}|{hap2[pos]}")
                    ps_out.append(str(ps_full[pos]))
                else:
                    phased_gt.append("0/1")
                    ps_out.append(".")

        _write_phased_vcf._seen = []
        _write_phased_vcf(vcf_path, out_vcf, phased_gt, ps_out, sample=sample)

    else:
        # no VCF: keep “matrix mode” output only
        H = np.stack([hap1_phase, hap2_phase], axis=0)


    # 7) Placeholder assignments: all selected reads assigned to haplotype 0
    assignments = np.zeros(len(selected_indices), dtype=np.int32)

    # 8) Write outputs in the same format as other algorithms
    prefix = args.output_prefix
    write_haplotypes_tsv(f"{prefix}.haplotypes.tsv", H)
    write_assignments_tsv(f"{prefix}.assignments.tsv", assignments)
    summary = {
        "algorithm": "diploid_whats",
        "R": int(R),
        "N": int(N),
        "max_coverage": int(args.max_coverage),
        "selected_reads": int(len(selected_indices)),
        "whatshap_module": os.path.realpath(wh.__file__),
        "whatshap_core_module": os.path.realpath(core.__file__),
        "whatshap_readselect_module": os.path.realpath(readselect.__file__),
    }

    if vcf_path:
        summary.update({
            "N_total": int(N),
            "N_het": int(len(het_positions)),
            "vcf_input": os.path.realpath(vcf_path),
            "vcf_output": os.path.realpath(out_vcf),
        })

    write_summary_json(summary, f"{prefix}.summary.json")


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
