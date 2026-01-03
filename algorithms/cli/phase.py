import argparse
import sys
from algorithms import diploid, polyploid

def diploid_em_main(args):
    from algorithms.diploid import em
    em.main(args)

def diploid_mst_main(args):
    from algorithms.diploid import mst
    mst.main(args)

def polyploid_em_main(args):
    from algorithms.polyploid import em
    em.main(args)

def polyploid_spectral_main(args):
    from algorithms.polyploid import spectral
    spectral.main(args)
    
def diploid_whats_main(args):
    from algorithms.diploid import whatshap_driver
    whatshap_driver.main(args)

def main():
    parser = argparse.ArgumentParser(description="Unified CLI for phasing algorithms")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # diploid-em
    dip_em = subparsers.add_parser("diploid-em", help="Diploid phasing using hard EM")
    dip_em.add_argument("-i", "--input", required=True, help="Input NPZ file")
    dip_em.add_argument("--output-prefix", required=True, help="Prefix for output files")
    dip_em.add_argument("--max-iters", type=int, default=20)
    dip_em.add_argument("--tol-iters", type=int, default=5)
    dip_em.add_argument("--seed", type=int, default=1)
    dip_em.set_defaults(func=diploid_em_main)

    # diploid-mst
    dip_mst = subparsers.add_parser("diploid-mst", help="Diploid phasing using MST")
    dip_mst.add_argument("-i", "--input", required=True, help="Input NPZ file")
    dip_mst.add_argument("--output-prefix", required=True, help="Prefix for output files")
    dip_mst.add_argument("--max-iters", type=int, default=20)
    dip_mst.add_argument("--tol-iters", type=int, default=5)
    dip_mst.add_argument("--seed", type=int, default=1)
    dip_mst.add_argument("--min-overlap", type=int, default=1)
    dip_mst.add_argument("--min-het-minor", type=int, default=1)
    dip_mst.set_defaults(func=diploid_mst_main)

    # polyploid-em
    poly_em = subparsers.add_parser("polyploid-em", help="Polyploid phasing using hard EM")
    poly_em.add_argument("-i", "--input", required=True, help="Input NPZ file")
    poly_em.add_argument("--ploidy", type=int, required=True, help="Ploidy (number of haplotypes)")
    poly_em.add_argument("--output-prefix", required=True, help="Prefix for output files")
    poly_em.add_argument("--max-iters", type=int, default=20)
    poly_em.add_argument("--tol-iters", type=int, default=5)
    poly_em.add_argument("--seed", type=int, default=1)
    poly_em.set_defaults(func=polyploid_em_main)

    # polyploid-spectral
    poly_spec = subparsers.add_parser("polyploid-spectral", help="Polyploid phasing using spectral clustering")
    poly_spec.add_argument("-i", "--input", required=True, help="Input NPZ file")
    poly_spec.add_argument("--output-prefix", required=True, help="Prefix for output files")
    poly_spec.add_argument("--max-iters", type=int, default=20)
    poly_spec.add_argument("--tol-iters", type=int, default=5)
    poly_spec.add_argument("--seed", type=int, default=1)
    poly_spec.add_argument("--min-overlap", type=int, default=1)
    poly_spec.add_argument("--ploidy", type=int, default=2)
    poly_spec.set_defaults(func=polyploid_spectral_main)

    # diploid-whatshap
    dip_wh = subparsers.add_parser("diploid-whats", help="Diploid phasing using WhatsHap core")
    dip_wh.add_argument("-i", "--input", required=True, help="Input NPZ file")
    dip_wh.add_argument("--output-prefix", required=True, help="Prefix for output files")
    dip_wh.add_argument("--max-coverage", type=int, default=15)
    dip_wh.add_argument("--error-rate", type=float, default=0.01)
    dip_wh.add_argument("--vcf", help="Input VCF with (unphased) genotypes. If set, only heterozygous sites are phased and homozygous sites are fixed by GT.")
    dip_wh.add_argument("--sample", help="Sample name in VCF (default: first sample).")
    dip_wh.add_argument("--output-vcf", help="Output phased VCF path (default: <output-prefix>.phased.vcf).")
    dip_wh.set_defaults(func=diploid_whats_main)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()