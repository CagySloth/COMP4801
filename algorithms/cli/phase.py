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

def main():
    parser = argparse.ArgumentParser(description="Unified CLI for phasing algorithms")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # diploid-em
    dip_em = subparsers.add_parser("diploid-em", help="Diploid phasing using hard EM")
    dip_em.add_argument("-i", "--input", required=True, help="Input reads .npz")
    dip_em.add_argument("-o", "--output", required=True, help="Output prefix")
    dip_em.add_argument("--max-iters", type=int, default=30)
    dip_em.add_argument("--tol-iters", type=int, default=2)
    dip_em.add_argument("-s", "--seed", type=int, default=None)
    dip_em.set_defaults(func=diploid_em_main)

    # diploid-mst
    dip_mst = subparsers.add_parser("diploid-mst", help="Diploid phasing using MST")
    dip_mst.add_argument("-i", "--input", required=True)
    dip_mst.add_argument("-o", "--output", required=True)
    dip_mst.add_argument("--min-overlap", type=int, default=3)
    dip_mst.add_argument("--min-het-minor", type=int, default=1)
    dip_mst.set_defaults(func=diploid_mst_main)

    # polyploid-em
    poly_em = subparsers.add_parser("polyploid-em", help="Polyploid phasing using hard EM")
    poly_em.add_argument("-i", "--input", required=True)
    poly_em.add_argument("-k", "--ploidy", type=int, required=True)
    poly_em.add_argument("-o", "--output", required=True)
    poly_em.add_argument("--max-iters", type=int, default=40)
    poly_em.add_argument("--tol-iters", type=int, default=3)
    poly_em.add_argument("-s", "--seed", type=int, default=None)
    poly_em.set_defaults(func=polyploid_em_main)

    # polyploid-spectral
    poly_spec = subparsers.add_parser("polyploid-spectral", help="Polyploid phasing using spectral clustering")
    poly_spec.add_argument("-i", "--input", required=True)
    poly_spec.add_argument("-k", "--ploidy", type=int, required=True)
    poly_spec.add_argument("-o", "--output", required=True)
    poly_spec.add_argument("--min-overlap", type=int, default=3)
    poly_spec.add_argument("-s", "--seed", type=int, default=None)
    poly_spec.set_defaults(func=polyploid_spectral_main)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()