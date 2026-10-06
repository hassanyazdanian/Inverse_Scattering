"""Locate a bracketed transmission eigenvalue of a homogeneous sphere."""
import argparse
from scipy.optimize import brentq
from sphere_scattering import transmission_determinant


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=float, default=40)
    parser.add_argument('--radius', type=float, default=1)
    parser.add_argument('--order', type=int, default=0)
    parser.add_argument('--bracket', nargs=2, type=float, default=(0.70, 0.73))
    args = parser.parse_args(argv)
    if args.n == 1:
        parser.error('n=1 has no contrast; the determinant vanishes identically')
    root = brentq(lambda k: transmission_determinant(k, args.n, args.radius, args.order),
                  *args.bracket, xtol=1e-13)
    print(f'n={args.n:g}, order={args.order}, k={root:.15g}')
    return root


if __name__ == '__main__':
    main()
