import os
import argparse
import time
import numpy as np

from algorithms.vendor.whatshap_vendor import import_whatshap_vendor
wh, core, readselect = import_whatshap_vendor()

from algorithms.io.reads_data import ReadsData
from algorithms.io.writer import (
    write_haplotypes_tsv,
    write_assignments_tsv,
    write_summary_json,
)

from algorithms.diploid.whatshap_adapter import build_readset_from_readsdata


def _read_vcf_gt_list(vcf_path: str, sample: str | None = None):
    """
    Return:
      gt_list: list[str] GT field as in VCF (e.g. '0/1', '1/1', './.')
      alleles: list[tuple[int|None,int|None]] parsed alleles (None if missing)
      is_het: list[bool]
    Notes:
      - This is a *minimal* VCF parser tailored to your simulator output.
      - It assumes records are biallelic and GT is diploid.
    """
    gt_list: list[str] = []
    alleles: list[tuple[int | None, int | None]] = []
    is_het: list[bool] = []
    sample_col: int | None = None

    with open(vcf_path, "r") as f:
        for line in f:
            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
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
            sval = fields[sample_col].split(":") if sample_col is not None else fields[9].split(":")
            d = dict(zip(fmt, sval))
            gt = d.get("GT", "./.")
            gt_list.append(gt)

            if gt in ("./.", ".|."):
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
    """Write a phased VCF by replacing GT (and adding PS) for one sample.

    This keeps all other fields untouched.
    """
    assert len(phased_gt) == len(ps)

    sample_col = None
    with open(in_vcf, "r") as fin, open(out_vcf, "w") as fout:
        for line in fin:
            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
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
            # ensure GT and PS are in FORMAT
            if "GT" not in fmt:
                fmt = ["GT"] + fmt
                sval = ["./."] + sval
            if "PS" not in fmt:
                fmt = fmt + ["PS"]
                sval = sval + ["."]
            d = dict(zip(fmt, sval))

            d["GT"] = phased_gt[rec_i]
            d["PS"] = ps[rec_i]

            new_sample = ":".join(d.get(k, ".") for k in fmt)
            fields[8] = ":".join(fmt)
            fields[sample_col] = new_sample
            fout.write("\t".join(fields) + "\n")


_write_phased_vcf._seen = []


class _UnionFind:
    def __init__(self, items: list[int]):
        self.parent = {x: x for x in items}
        self.rank = {x: 0 for x in items}

    def find(self, x: int) -> int:
        p = self.parent[x]
        if p != x:
            self.parent[x] = self.find(p)
        return self.parent[x]

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


def _component_leftmost(positions: list[int], readset: core.ReadSet) -> dict[int, int]:
    """Map each position -> leftmost position in its connected component.

    Component graph: positions are nodes; an edge exists between two positions if at least one read
    covers both positions. This mirrors WhatsHap's PS definition.
    """
    if not positions:
        return {}
    posset = set(positions)
    uf = _UnionFind(positions)

    for read in readset:
        covered = [v.position for v in read if v.position in posset]
        if len(covered) < 2:
            continue
        anchor = covered[0]
        for p in covered[1:]:
            uf.union(anchor, p)

    # compute leftmost per root
    leftmost_by_root: dict[int, int] = {}
    for p in positions:
        r = uf.find(p)
        leftmost_by_root[r] = p if r not in leftmost_by_root else min(leftmost_by_root[r], p)

    return {p: leftmost_by_root[uf.find(p)] for p in positions}


def main(args=None):
    parser = argparse.ArgumentParser(description="Diploid phasing via WhatsHap core")
    parser.add_argument("-i", "--input", required=True, help="Input NPZ file")
    parser.add_argument("--output-prefix", required=True, help="Prefix for output files")
    parser.add_argument("--max-coverage", type=int, default=15)
    parser.add_argument("--error-rate", type=float, default=0.01)

    parser.add_argument("--vcf", help="Input VCF with genotypes (unphased). If set, phase only het sites and write phased VCF.")
    parser.add_argument("--sample", help="Sample name in VCF (default: first sample).")
    parser.add_argument("--output-vcf", help="Output phased VCF path (default: <output-prefix>.phased.vcf).")

    # NEW: solver selection in VCF mode
    parser.add_argument(
        "--solver",
        choices=["whatshap", "hapchat"],
        default=None,
        help="Which solver to use in VCF-mode. 'whatshap' uses PedigreeDPTable (default WhatsHap). "
             "'hapchat' uses HapChatCore MEC solver.",
    )
    parser.add_argument(
        "--recomb-rate",
        type=float,
        default=1.26,
        help="Recombination rate (cM/Mb) used by UniformRecombinationCostComputer (PedigreeDPTable only).",
    )

    if args is None:
        args = parser.parse_args()

    t_total_start = time.perf_counter()

    # 1) Load reads
    t_load_start = time.perf_counter()
    data = ReadsData.from_npz(args.input)
    t_load_sec = time.perf_counter() - t_load_start
    N = data.N
    R = data.R

    vcf_path = getattr(args, "vcf", None)
    sample = getattr(args, "sample", None)
    out_vcf = getattr(args, "output_vcf", None) or f"{args.output_prefix}.phased.vcf"

    # Resolve solver choice:
    solver = getattr(args, "solver", None) or ("whatshap" if vcf_path else "hapchat")

    # Optional VCF parsing (for GT + het filter)
    het_positions = None
    old_to_new = None
    alleles = None
    is_het = None

    data_for_phasing = data

    t_vcf_parse_sec = 0.0
    t_build_readset_sec = 0.0
    t_readselect_sec = 0.0
    t_solve_sec = 0.0

    if vcf_path:
        t_vcf_start = time.perf_counter()
        _, alleles, is_het = _read_vcf_gt_list(vcf_path, sample=sample)
        t_vcf_parse_sec = time.perf_counter() - t_vcf_start
        if len(is_het) != N:
            raise ValueError(
                f"VCF variant count ({len(is_het)}) != reads matrix N ({N}). "
                "They must match in this simulation pipeline."
            )

        het_positions = [i for i, h in enumerate(is_het) if h]
        old_to_new = {old: new for new, old in enumerate(het_positions)}

        # Build het-only matrix (like real WhatsHap: phase only heterozygous variants).
        A_het = data.alleles[:, het_positions] if len(het_positions) > 0 else np.empty((R, 0), dtype=int)

        # IMPORTANT: keep ORIGINAL positions (0-based, aligned to VCF record index)
        # so that recombination costs and PS values match "default WhatsHap" semantics.
        positions_het = np.tile(np.array(het_positions, dtype=int), (R, 1))
        data_for_phasing = ReadsData(reads=A_het, positions=positions_het, num_variants=A_het.shape[1])

    # 2) Build ReadSet
    t_rs_start = time.perf_counter()
    readset = build_readset_from_readsdata(data_for_phasing, use_positions=bool(vcf_path))
    readset.sort()
    t_build_readset_sec = time.perf_counter() - t_rs_start

    # --- Filter non-informative reads (must cover >= 2 variants) ---
    informative_idx = [i for i, r in enumerate(readset) if len(r) >= 2]
    informative_readset = readset.subset(informative_idx)
    informative_readset.sort()

    # If there aren't enough variants or informative reads, phasing is impossible.
    N_phase = data_for_phasing.N
    selected_indices: list[int] = []
    selected_readset = None
    phase_source = None  # ("blocks", list[ReadSet]) or ("pair", ReadSet)

    if N_phase >= 2 and len(informative_readset) > 0:
        # 3) Perform WhatsHap read selection
        t_sel_start = time.perf_counter()
        sel_local = readselect.readselection(informative_readset, args.max_coverage, None)
        t_readselect_sec = time.perf_counter() - t_sel_start

        selected_readset = informative_readset.subset(sel_local)
        selected_indices = [informative_idx[i] for i in sel_local]

        # 4) Solve
        t_solve_start = time.perf_counter()
        if vcf_path and solver == "whatshap":
            # Default WhatsHap path: PedigreeDPTable
            positions = sorted(selected_readset.get_positions())
            genotypes = [core.Genotype([int(alleles[p][0]), int(alleles[p][1])]) for p in positions]

            numeric_sample_ids = core.NumericSampleIds()
            sample_name = sample or "SAMPLE"
            pedigree = core.Pedigree(numeric_sample_ids)
            pedigree.add_individual(sample_name, genotypes)

            try:
                from whatshap.pedigree import UniformRecombinationCostComputer  # type: ignore
                recombcost = UniformRecombinationCostComputer(float(args.recomb_rate)).compute(positions)
            except Exception:
                # If pedigree module isn't vendored yet, we can still benchmark DP without recombination
                recombcost = [0] * len(positions)

            dp = core.PedigreeDPTable(selected_readset, recombcost, pedigree, False, positions)
            sr_by_individual, _tv = dp.get_super_reads()
            superreads = sr_by_individual[0]  # ReadSet with 2 reads
            phase_source = ("pair", superreads)
        else:
            # HapChat MEC solver
            hap_core = core.HapChatCore(selected_readset)
            blocks, _ = hap_core.get_super_reads()
            phase_source = ("blocks", blocks)

        t_solve_sec = time.perf_counter() - t_solve_start

    # --------- 5) Extract haplotypes + PS ---------
    if vcf_path:
        # FULL space output (N sites) because your downstream evaluation expects N-length haplotypes.tsv
        hap1 = np.full(N, -1, dtype=int)
        hap2 = np.full(N, -1, dtype=int)
        ps_full = np.full(N, -1, dtype=int)

        # Fill homozygous sites from GT and give them no PS (WhatsHap doesn't phase them)
        for pos in range(N):
            a, b = alleles[pos]
            if a is None or b is None:
                continue
            if not is_het[pos]:
                hap1[pos] = hap2[pos] = int(a)

        # Fill het sites that were actually solved
        if phase_source is not None:
            kind, payload = phase_source
            if kind == "pair":
                sr = payload
                if len(sr) >= 2:
                    for v in sr[0]:
                        hap1[v.position] = int(v.allele)
                    for v in sr[1]:
                        hap2[v.position] = int(v.allele)
            else:
                blocks = payload
                for block in blocks:
                    if len(block) < 2:
                        continue
                    for v in block[0]:
                        hap1[v.position] = int(v.allele)
                    for v in block[1]:
                        hap2[v.position] = int(v.allele)

        # Compute PS from connectivity of the *selected reads* (WhatsHap definition)
        if selected_readset is not None:
            accessible_positions = sorted(selected_readset.get_positions())
            leftmost = _component_leftmost(accessible_positions, selected_readset)
            for pos in accessible_positions:
                ps_full[pos] = leftmost[pos] + 1  # VCF PS is 1-based

        # Fallback fill to ensure haplotypes.tsv contains only 0/1 (your test/eval requirement).
        # NOTE: This is NOT "default WhatsHap" output; WhatsHap would leave unconnected hets unphased.
        if het_positions is not None and old_to_new is not None:
            A_het = data_for_phasing.alleles  # R x N_het
            for pos in het_positions:
                if hap1[pos] == -1 and hap2[pos] == -1:
                    j = old_to_new[pos]
                    col = A_het[:, j]
                    obs = sorted({int(a) for a in col if a >= 0})
                    if len(obs) == 1:
                        hap1[pos] = hap2[pos] = obs[0]
                    elif len(obs) >= 2:
                        hap1[pos], hap2[pos] = obs[0], obs[1]
                    else:
                        hap1[pos] = hap2[pos] = 0
                elif hap1[pos] == -1:
                    hap1[pos] = 1 - hap2[pos]
                elif hap2[pos] == -1:
                    hap2[pos] = 1 - hap1[pos]

        H = np.stack([hap1, hap2], axis=0)

        phased_gt: list[str] = []
        ps_out: list[str] = []
        het_phased = 0
        het_unphased = 0

        for pos in range(N):
            a, b = alleles[pos]
            if a is None or b is None:
                phased_gt.append("./.")
                ps_out.append(".")
                continue

            if not is_het[pos]:
                phased_gt.append(f"{int(a)}/{int(b)}")
                ps_out.append(".")
                continue

            if ps_full[pos] >= 0:
                phased_gt.append(f"{int(hap1[pos])}|{int(hap2[pos])}")
                ps_out.append(str(int(ps_full[pos])))
                het_phased += 1
            else:
                phased_gt.append("0/1")
                ps_out.append(".")
                het_unphased += 1

        _write_phased_vcf._seen = []
        _write_phased_vcf(vcf_path, out_vcf, phased_gt, ps_out, sample=sample)

        num_phase_sets = len(set(int(x) for x in ps_full if x >= 0))

    else:
        # MATRIX MODE: treat columns as positions 0..N-1
        hap1 = np.full(N, -1, dtype=int)
        hap2 = np.full(N, -1, dtype=int)
        ps = np.full(N, -1, dtype=int)

        if phase_source is not None:
            kind, payload = phase_source
            if kind == "pair":
                sr = payload
                if len(sr) >= 2:
                    for v in sr[0]:
                        hap1[v.position] = int(v.allele)
                    for v in sr[1]:
                        hap2[v.position] = int(v.allele)
            else:
                blocks = payload
                for block in blocks:
                    if len(block) < 2:
                        continue
                    for v in block[0]:
                        hap1[v.position] = int(v.allele)
                    for v in block[1]:
                        hap2[v.position] = int(v.allele)

        if selected_readset is not None:
            positions = sorted(selected_readset.get_positions())
            leftmost = _component_leftmost(positions, selected_readset)
            for pos in positions:
                ps[pos] = leftmost[pos] + 1

        # Fallback fill (keep your existing behavior for TSV)
        A_phase = data_for_phasing.alleles
        for pos in range(N):
            if hap1[pos] == -1 and hap2[pos] == -1:
                col = A_phase[:, pos]
                obs = sorted({int(a) for a in col if a >= 0})
                if len(obs) == 1:
                    hap1[pos] = hap2[pos] = obs[0]
                elif len(obs) >= 2:
                    hap1[pos], hap2[pos] = obs[0], obs[1]
                else:
                    hap1[pos] = hap2[pos] = 0
            elif hap1[pos] == -1:
                hap1[pos] = 1 - hap2[pos]
            elif hap2[pos] == -1:
                hap2[pos] = 1 - hap1[pos]

        H = np.stack([hap1, hap2], axis=0)
        num_phase_sets = len(set(int(x) for x in ps if x >= 0))

    # 7) Placeholder assignments: all selected reads assigned to haplotype 0
    assignments = np.zeros(len(selected_indices), dtype=np.int32)

    # 8) Write outputs
    prefix = args.output_prefix
    write_haplotypes_tsv(f"{prefix}.haplotypes.tsv", H)
    write_assignments_tsv(f"{prefix}.assignments.tsv", assignments)

    t_total_sec = time.perf_counter() - t_total_start
    summary = {
        "algorithm": "diploid_whats",
        "solver": solver,
        "R": int(R),
        "N": int(N),
        "max_coverage": int(args.max_coverage),

        "num_reads_total": int(len(readset)),
        "num_informative_reads": int(len(informative_readset)),
        "selected_reads": int(len(selected_indices)),
        "num_phase_sets": int(num_phase_sets),

        "time_total_sec": float(t_total_sec),
        "time_load_sec": float(t_load_sec),
        "time_vcf_parse_sec": float(t_vcf_parse_sec),
        "time_build_readset_sec": float(t_build_readset_sec),
        "time_readselection_sec": float(t_readselect_sec),
        "time_solve_sec": float(t_solve_sec),

        "whatshap_module": os.path.realpath(wh.__file__),
        "whatshap_core_module": os.path.realpath(core.__file__),
        "whatshap_readselect_module": os.path.realpath(readselect.__file__),
    }

    if vcf_path:
        summary.update({
            "het_total": int(sum(is_het)),
            "vcf_input": os.path.realpath(vcf_path),
            "vcf_output": os.path.realpath(out_vcf),
        })

    write_summary_json(summary, f"{prefix}.summary.json")


if __name__ == "__main__":
    main()
