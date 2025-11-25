import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from algorithms.io.writer import (
    write_haplotypes_tsv,
    write_reads_sparse_tsv,
    write_haplotypes_npz,
)
from algorithms.io.reads_data import ReadsData


# -------------------------------
#  Technology presets / error model
# -------------------------------

@dataclass
class TechModel:
    name: str
    error_rate_simple: float      # per-site flip rate in "easy" regions
    error_rate_hard: float        # per-site flip rate in "hard" regions
    missing_rate_simple: float    # per-site missing (-1) in "easy" regions
    missing_rate_hard: float      # per-site missing (-1) in "hard" regions
    recall_simple: float          # probability variant is "called" in easy regions
    recall_hard: float            # probability variant is "called" in hard regions


TECH_MODELS = {
    "hifi": TechModel(
        name="hifi",
        error_rate_simple=0.003,
        error_rate_hard=0.02,
        missing_rate_simple=0.01,
        missing_rate_hard=0.05,
        recall_simple=0.99,
        recall_hard=0.93,
    ),
    "ont": TechModel(
        name="ont",
        error_rate_simple=0.01,
        error_rate_hard=0.08,
        missing_rate_simple=0.02,
        missing_rate_hard=0.08,
        recall_simple=0.97,
        recall_hard=0.90,
    ),
}


# -------------------------------
#  Core simulator
# -------------------------------

class LongReadReadsetSimulator:
    """
    Simulate:
      - underlying haplotypes (truth)
      - per-variant "context" (easy vs hard)
      - reads at the variant level, with tech-specific error/missing/recall.

    Output 'reads' is a list of dicts:
      {"id": int, "indices": [int...], "values": [0/1/-1...]}

    That matches what your current pipeline expects.
    """

    def __init__(
        self,
        ploidy: int,
        num_variants: int,
        num_reads: int,
        mean_read_length: int,
        read_length_sd: float,
        tech_model: TechModel,
        maf_alpha: float = 1.0,
        maf_beta: float = 1.0,
        hard_region_fraction: float = 0.1,
        allow_monomorphic: bool = False,
        seed: int | None = None,
    ):
        self.ploidy = int(ploidy)
        self.num_variants = int(num_variants)
        self.num_reads = int(num_reads)
        self.mean_read_length = int(mean_read_length)
        self.read_length_sd = float(read_length_sd)
        self.tech = tech_model
        self.maf_alpha = float(maf_alpha)
        self.maf_beta = float(maf_beta)
        self.hard_region_fraction = float(hard_region_fraction)
        self.allow_monomorphic = bool(allow_monomorphic)
        self.rng = np.random.default_rng(seed)

    # ---- public entry point ----

    def simulate(self):
        """
        Returns:
            haplotypes: (ploidy, num_variants) int8 array of 0/1
            reads: list[{"id", "indices", "values"}]
            contexts: (num_variants,) int8 array; 0 = easy, 1 = hard
        """
        contexts = self._simulate_variant_contexts()
        haplotypes = self._simulate_haplotypes()
        reads = self._simulate_reads(haplotypes, contexts)
        return haplotypes, reads, contexts

    # ---- internal steps ----

    def _simulate_variant_contexts(self):
        """
        0 = "simple" context
        1 = "hard" context (homopolymer/STR/segdup-like)
        """
        hard = int(round(self.num_variants * self.hard_region_fraction))
        ctx = np.zeros(self.num_variants, dtype=np.int8)
        if hard > 0:
            hard_indices = self.rng.choice(self.num_variants, size=hard, replace=False)
            ctx[hard_indices] = 1
        return ctx

    def _simulate_haplotypes(self):
        if self.ploidy == 2:
            return self._simulate_diploid_haplotypes()
        else:
            return self._simulate_polyploid_haplotypes()

    def _simulate_diploid_haplotypes(self):
        hap1 = self.rng.integers(0, 2, size=self.num_variants, dtype=np.int8)
        hap2 = hap1.copy()

        # Local switches to mimic hap-block structure
        switch_indices = self.rng.random(self.num_variants) < 0.1
        hap2[switch_indices] ^= 1

        if not self.allow_monomorphic and np.all(hap1 == hap2):
            hap2[0] ^= 1

        return np.stack([hap1, hap2], axis=0)

    def _simulate_polyploid_haplotypes(self):
        # per-site allele freq from Beta; then ploidy Bernoulli draws
        p = self.rng.beta(self.maf_alpha, self.maf_beta, size=self.num_variants)
        haplotypes = self.rng.binomial(
            1, p, size=(self.ploidy, self.num_variants)
        ).astype(np.int8)

        if not self.allow_monomorphic:
            # ensure at least one differing site
            while np.all(haplotypes == haplotypes[0]):
                i = self.rng.integers(self.ploidy)
                j = self.rng.integers(self.num_variants)
                haplotypes[i, j] ^= 1

        return haplotypes

    def _draw_read_length(self) -> int:
        if self.read_length_sd <= 0:
            L = self.mean_read_length
        else:
            L = int(round(self.rng.normal(self.mean_read_length, self.read_length_sd)))

        if L < 1:
            L = 1
        if L > self.num_variants:
            L = self.num_variants
        return L

    def _simulate_reads(self, haplotypes: np.ndarray, contexts: np.ndarray):
        reads = []
        for read_id in range(self.num_reads):
            L = self._draw_read_length()
            start = self.rng.integers(0, self.num_variants - L + 1)
            hap_idx = self.rng.integers(0, self.ploidy)

            # truth alleles for this read
            segment = haplotypes[hap_idx, start : start + L].copy()
            ctx_segment = contexts[start : start + L]

            # random uniforms for each type of noise
            u_error = self.rng.random(L)
            u_missing = self.rng.random(L)
            u_recall = self.rng.random(L)

            for i in range(L):
                hard = ctx_segment[i] == 1

                # 1) Variant calling "recall": drop entire site as uncalled
                recall_thr = self.tech.recall_hard if hard else self.tech.recall_simple
                if u_recall[i] > recall_thr:
                    segment[i] = -1
                    continue

                # 2) Sequencing/genotyping error: flip allele
                err_thr = self.tech.error_rate_hard if hard else self.tech.error_rate_simple
                if segment[i] != -1 and u_error[i] < err_thr and segment[i] in (0, 1):
                    segment[i] ^= 1

                # 3) Read-level missingness (low qual / mapping issues)
                miss_thr = self.tech.missing_rate_hard if hard else self.tech.missing_rate_simple
                if u_missing[i] < miss_thr:
                    segment[i] = -1

            indices = list(range(start, start + L))
            values = segment.astype(int).tolist()
            reads.append(
                {
                    "id": int(read_id),
                    "indices": indices,
                    "values": values,
                }
            )

        return reads


# -------------------------------
#  CLI wrapper (same I/O as current script)
# -------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Simulate long-read-style haplotype and read data (WhatsHap-like ReadSet)."
    )
    parser.add_argument("-p", "--ploidy", type=int, required=True)
    parser.add_argument("-n", "--num-variants", type=int, required=True)
    parser.add_argument("-r", "--num-reads", type=int, required=True)
    parser.add_argument(
        "-l",
        "--read-length",
        type=int,
        required=True,
        help="Mean read length in #variants.",
    )
    parser.add_argument(
        "--read-length-sd",
        type=float,
        default=0.0,
        help="SD of read length in #variants (0 = fixed).",
    )

    parser.add_argument(
        "-t",
        "--tech",
        choices=sorted(TECH_MODELS.keys()),
        default="hifi",
        help="Technology preset controlling error/missing/recall rates.",
    )
    parser.add_argument(
        "-e",
        "--error-rate",
        type=float,
        default=None,
        help="Override SIMPLE-region error rate (per-allele).",
    )
    parser.add_argument(
        "-m",
        "--missing-rate",
        type=float,
        default=None,
        help="Override SIMPLE-region missing rate.",
    )

    parser.add_argument("--maf-alpha", type=float, default=1.0)
    parser.add_argument("--maf-beta", type=float, default=1.0)
    parser.add_argument(
        "--hard-region-fraction",
        type=float,
        default=0.1,
        help="Fraction of variants marked as 'hard' context.",
    )
    parser.add_argument("--allow-monomorphic", action="store_true")
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("-o", "--output-prefix", required=True)

    args = parser.parse_args()

    tech_model = TECH_MODELS[args.tech]

    # Optional overrides so you can still do simple "global error rate" experiments
    if args.error_rate is not None:
        tech_model = TechModel(
            name=tech_model.name,
            error_rate_simple=args.error_rate,
            error_rate_hard=tech_model.error_rate_hard,
            missing_rate_simple=tech_model.missing_rate_simple,
            missing_rate_hard=tech_model.missing_rate_hard,
            recall_simple=tech_model.recall_simple,
            recall_hard=tech_model.recall_hard,
        )
    if args.missing_rate is not None:
        tech_model = TechModel(
            name=tech_model.name,
            error_rate_simple=tech_model.error_rate_simple,
            error_rate_hard=tech_model.error_rate_hard,
            missing_rate_simple=args.missing_rate,
            missing_rate_hard=tech_model.missing_rate_hard,
            recall_simple=tech_model.recall_simple,
            recall_hard=tech_model.recall_hard,
        )

    sim = LongReadReadsetSimulator(
        ploidy=args.ploidy,
        num_variants=args.num_variants,
        num_reads=args.num_reads,
        mean_read_length=args.read_length,
        read_length_sd=args.read_length_sd,
        tech_model=tech_model,
        maf_alpha=args.maf_alpha,
        maf_beta=args.maf_beta,
        hard_region_fraction=args.hard_region_fraction,
        allow_monomorphic=args.allow_monomorphic,
        seed=args.seed,
    )

    haplotypes, reads, contexts = sim.simulate()

    prefix = Path(args.output_prefix)
    write_haplotypes_tsv(f"{prefix}.haplotypes.tsv", haplotypes)
    write_reads_sparse_tsv(f"{prefix}.reads.sparse.tsv", reads)
    reads_np = ReadsData.from_fragments(reads)
    write_haplotypes_npz(f"{prefix}.reads.npz", reads_np.reads)

    # Optional extra: per-variant context labelling (0 = easy, 1 = hard)
    np.savetxt(f"{prefix}.variant_contexts.txt", contexts, fmt="%d")

    print(
        f"✅ Simulated long-read ReadSet written to: "
        f"{prefix}.[haplotypes.tsv, reads.sparse.tsv, reads.npz, variant_contexts.txt]"
    )


if __name__ == "__main__":
    main()
