import argparse
import gzip
import os
import time
from collections import defaultdict
from typing import Dict, List, Tuple

import numpy as np

from algorithms.io.writer import write_haplotypes_tsv, write_assignments_tsv, write_summary_json
from algorithms.vendor.whatshap_vendor import import_whatshap_vendor

wh, core, readselect = import_whatshap_vendor()


def _open_text(path: str, mode: str):
    # mode: "r" or "w"
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode)


def _read_vcf_minimal(vcf_path: str, sample: str | None = None):
    """
    Minimal VCF parser that supports .vcf and .vcf.gz.
    Returns:
      records: list of (chrom, pos0, ref, alt, gt_str, a, b, is_het)
    Notes:
      - assumes single-sample VCF
      - assumes biallelic SNPs (REF/ALT length 1) for now
    """
    records = []
    sample_col = None

    with _open_text(vcf_path, "r") as f:
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
                        raise ValueError(f"Sample '{sample}' not found. Available: {samples}")
                    sample_col = header.index(sample)
                continue
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            pos1 = int(fields[1])
            pos0 = pos1 - 1
            ref = fields[3].upper()
            alt = fields[4].upper()

            # keep SNPs only (for now)
            if len(ref) != 1 or len(alt) != 1:
                continue

            fmt = fields[8].split(":")
            sval = fields[sample_col].split(":")
            d = dict(zip(fmt, sval))
            gt = d.get("GT", "./.")

            if gt in ("./.", ".|."):
                a = b = None
                is_het = False
            else:
                sep = "|" if "|" in gt else "/"
                a_str, b_str = gt.split(sep)
                if a_str == "." or b_str == ".":
                    a = b = None
                    is_het = False
                else:
                    a = int(a_str)
                    b = int(b_str)
                    is_het = (a != b)

            records.append((chrom, pos0, ref, alt, gt, a, b, is_het))

    return records


def _write_phased_vcf(in_vcf: str, out_vcf: str, phased_gt: List[str], ps: List[str], sample: str | None = None):
    """
    Stream input VCF and replace GT + add PS for one sample.
    Supports .vcf and .vcf.gz.
    """
    assert len(phased_gt) == len(ps)
    saw_ps_meta = False
    sample_col = None
    rec_i = 0

    with _open_text(in_vcf, "r") as fin, _open_text(out_vcf, "w") as fout:
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
                        raise ValueError(f"Sample '{sample}' not found. Available: {samples}")
                    sample_col = header.index(sample)
                fout.write(line)
                continue

            if line.startswith("#"):
                fout.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            fmt = fields[8].split(":")
            sval = fields[sample_col].split(":")

            # ensure GT exists
            if "GT" not in fmt:
                fmt = ["GT"] + fmt
                sval = ["./."] + sval

            # ensure PS exists
            if "PS" not in fmt:
                fmt.append("PS")
                sval.append(".")

            d = dict(zip(fmt, sval))
            if rec_i < len(phased_gt):
                d["GT"] = phased_gt[rec_i]
                d["PS"] = ps[rec_i]
            else:
                # If this record was not part of phasing (e.g., filtered out earlier),
                # keep original GT and write PS as '.'
                d["GT"] = d.get("GT", "./.")
                d["PS"] = "."

            fields[8] = ":".join(fmt)
            fields[sample_col] = ":".join(d.get(k, ".") for k in fmt)
            fout.write("\t".join(fields) + "\n")
            rec_i += 1


class _UnionFind:
    def __init__(self, items):
        self.parent = {x: x for x in items}
        self.rank = {x: 0 for x in items}

    def find(self, x):
        p = self.parent[x]
        if p != x:
            self.parent[x] = self.find(p)
        return self.parent[x]

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


def _component_leftmost(positions, readset: core.ReadSet):
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

    leftmost_by_root = {}
    for p in positions:
        r = uf.find(p)
        leftmost_by_root[r] = p if r not in leftmost_by_root else min(leftmost_by_root[r], p)

    return {p: leftmost_by_root[uf.find(p)] for p in positions}


def _build_readset_from_bam_vcf(
    bam_path: str,
    het_variants: List[Tuple[str, int, str, str]],
    *,
    min_mapq: int = 20,
    min_baseq: int = 20,
    default_qual: int = 40,
) -> core.ReadSet:
    """
    Build WhatsHap ReadSet by extracting allele observations at heterozygous SNP sites.
    Variant positions used inside ReadSet are 0-based reference coords (pos0).
    """
    try:
        import pysam
    except ImportError as e:
        raise RuntimeError("pysam is required for BAM parsing. Run: pip install pysam") from e

    bam = pysam.AlignmentFile(bam_path, "rb")

    reads: Dict[str, core.Read] = {}
    # (optional) avoid duplicate add for same (read,pos)
    seen = set()

    for chrom, pos0, ref, alt in het_variants:
        for col in bam.pileup(chrom, pos0, pos0 + 1, truncate=True, stepper="all"):
            if col.reference_pos != pos0:
                continue
            for pr in col.pileups:
                if pr.is_del or pr.is_refskip or pr.query_position is None:
                    continue
                aln = pr.alignment
                if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                    continue
                if aln.mapping_quality < min_mapq:
                    continue

                qpos = pr.query_position
                base = aln.query_sequence[qpos].upper()
                qual = default_qual
                if aln.query_qualities is not None:
                    qual = int(aln.query_qualities[qpos])
                if qual < min_baseq:
                    continue

                if base == ref:
                    allele = 0
                elif base == alt:
                    allele = 1
                else:
                    continue

                key = (aln.query_name, pos0)
                if key in seen:
                    continue
                seen.add(key)

                r = reads.get(aln.query_name)
                if r is None:
                    r = core.Read(
                        name=aln.query_name,
                        mapq=int(aln.mapping_quality),
                        source_id=0,
                        sample_id=0,
                    )
                    reads[aln.query_name] = r

                r.add_variant(int(pos0), int(allele), int(qual))

    readset = core.ReadSet()
    for r in reads.values():
        r.sort()
        if len(r) >= 2:  # WhatsHap readselection requires >=2 variants covered
            readset.add(r)

    readset.sort()
    return readset


def main(args=None):
    parser = argparse.ArgumentParser(description="WhatsHap-core phasing from BAM+VCF (vendored core)")
    parser.add_argument("--bam", required=True, help="Input BAM (sorted, indexed)")
    parser.add_argument("--vcf", required=True, help="Input VCF/VCF.GZ with genotypes (called, unphased)")
    parser.add_argument("--output-prefix", required=True, help="Prefix for outputs")
    parser.add_argument("--sample", default=None, help="Sample name in VCF (default: first sample)")
    parser.add_argument("--output-vcf", default=None, help="Output phased VCF path (default: <prefix>.phased.vcf)")
    parser.add_argument("--max-coverage", type=int, default=15)
    parser.add_argument("--solver", choices=["whatshap", "hapchat"], default="whatshap")
    parser.add_argument("--recomb-rate", type=float, default=1.26)

    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-baseq", type=int, default=20)

    if args is None:
        args = parser.parse_args()

    t0 = time.perf_counter()

    out_vcf = args.output_vcf or f"{args.output_prefix}.phased.vcf"

    # 1) Parse VCF, keep het SNPs
    records = _read_vcf_minimal(args.vcf, sample=args.sample)
    het_variants = [(c, p0, r, a) for (c, p0, r, a, gt, A, B, is_het) in records if is_het and A is not None and B is not None]
    gt_by_pos = {p0: (A, B) for (c, p0, r, a, gt, A, B, is_het) in records if (A is not None and B is not None)}

    # 2) Build ReadSet from BAM+VCF het sites
    readset = _build_readset_from_bam_vcf(
        args.bam,
        het_variants,
        min_mapq=int(args.min_mapq),
        min_baseq=int(args.min_baseq),
    )

    # If not enough data, still emit outputs but unphased
    selected_indices = []
    selected_readset = None
    phase_source = None

    if len(readset) > 0 and len(het_variants) >= 2:
        # 3) read selection
        sel = readselect.readselection(readset, int(args.max_coverage), None)
        selected_readset = readset.subset(sel)
        selected_indices = list(sel)

        # 4) Solve
        if args.solver == "whatshap":
            positions = sorted(selected_readset.get_positions())
            genotypes = [core.Genotype([int(gt_by_pos[p][0]), int(gt_by_pos[p][1])]) for p in positions]

            numeric_sample_ids = core.NumericSampleIds()
            sample_name = args.sample or "SAMPLE"
            pedigree = core.Pedigree(numeric_sample_ids)
            pedigree.add_individual(sample_name, genotypes)

            try:
                from whatshap.pedigree import UniformRecombinationCostComputer  # type: ignore
                recombcost = UniformRecombinationCostComputer(float(args.recomb_rate)).compute(positions)
            except Exception:
                recombcost = [0] * len(positions)

            dp = core.PedigreeDPTable(selected_readset, recombcost, pedigree, False, positions)
            sr_by_individual, _ = dp.get_super_reads()
            superreads = sr_by_individual[0]  # ReadSet with 2 reads
            phase_source = ("pair", superreads)
        else:
            hc = core.HapChatCore(selected_readset)
            blocks, _ = hc.get_super_reads()
            phase_source = ("blocks", blocks)

    # 5) Extract hap alleles for phased het sites
    hap1_by_pos = {}
    hap2_by_pos = {}

    if phase_source is not None:
        kind, payload = phase_source
        if kind == "pair":
            sr = payload
            if len(sr) >= 2:
                for v in sr[0]:
                    hap1_by_pos[int(v.position)] = int(v.allele)
                for v in sr[1]:
                    hap2_by_pos[int(v.position)] = int(v.allele)
        else:
            for block in payload:
                if len(block) < 2:
                    continue
                for v in block[0]:
                    hap1_by_pos[int(v.position)] = int(v.allele)
                for v in block[1]:
                    hap2_by_pos[int(v.position)] = int(v.allele)

    # 6) Compute PS from connectivity on selected reads
    ps_by_pos = {}
    if selected_readset is not None:
        pos_list = sorted(selected_readset.get_positions())
        leftmost = _component_leftmost(pos_list, selected_readset)
        ps_by_pos = {p: int(leftmost[p]) + 1 for p in pos_list}  # PS is 1-based coordinate

    # 7) Build phased GT + PS per VCF record order + write outputs
    phased_gt = []
    ps_out = []

    # also emit haplotypes.tsv in "variant-index order"
    M = len(records)
    H1 = np.zeros(M, dtype=int)
    H2 = np.zeros(M, dtype=int)

    for i, (chrom, pos0, ref, alt, gt, A, B, is_het) in enumerate(records):
        if A is None or B is None:
            phased_gt.append("./.")
            ps_out.append(".")
            H1[i] = 0
            H2[i] = 0
            continue

        if not is_het:
            phased_gt.append(f"{int(A)}/{int(B)}")
            ps_out.append(".")
            H1[i] = int(A)
            H2[i] = int(B)
            continue

        if pos0 in hap1_by_pos and pos0 in hap2_by_pos and pos0 in ps_by_pos:
            phased_gt.append(f"{hap1_by_pos[pos0]}|{hap2_by_pos[pos0]}")
            ps_out.append(str(ps_by_pos[pos0]))
            H1[i] = hap1_by_pos[pos0]
            H2[i] = hap2_by_pos[pos0]
        else:
            # unphased het (match WhatsHap semantics)
            phased_gt.append("0/1")
            ps_out.append(".")
            # fill TSV deterministically (compat with your eval scripts)
            H1[i] = 0
            H2[i] = 1

    _write_phased_vcf(args.vcf, out_vcf, phased_gt, ps_out, sample=args.sample)

    H = np.stack([H1, H2], axis=0)
    write_haplotypes_tsv(f"{args.output_prefix}.haplotypes.tsv", H)

    assignments = np.zeros(len(selected_indices), dtype=np.int32)
    write_assignments_tsv(f"{args.output_prefix}.assignments.tsv", assignments)

    summary = {
        "algorithm": "diploid_whats_bam",
        "solver": args.solver,
        "vcf_input": os.path.realpath(args.vcf),
        "bam_input": os.path.realpath(args.bam),
        "vcf_output": os.path.realpath(out_vcf),
        "max_coverage": int(args.max_coverage),
        "min_mapq": int(args.min_mapq),
        "min_baseq": int(args.min_baseq),

        "variants_total": int(len(records)),
        "het_variants": int(len(het_variants)),

        "num_reads_total": int(len(readset)),
        "selected_reads": int(len(selected_indices)),
        "num_phase_sets": int(len(set(ps_by_pos.values()))) if ps_by_pos else 0,

        "time_total_sec": float(time.perf_counter() - t0),

        "whatshap_module": os.path.realpath(wh.__file__),
        "whatshap_core_module": os.path.realpath(core.__file__),
        "whatshap_readselect_module": os.path.realpath(readselect.__file__),
    }
    write_summary_json(summary, f"{args.output_prefix}.summary.json")
