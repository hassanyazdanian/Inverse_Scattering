"""Scan regularized far-field solution norms for a penetrable sphere."""
import argparse
import numpy as np
from geometry import sphere_quadrature
from regularization import tikhonov, weighted_far_field_svd
from lsm_common import add_common_arguments, prepare_mesh, assemble_far_field


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_arguments(parser, default_n=4)
    parser.add_argument('--k-min', type=float, default=2.8)
    parser.add_argument('--k-max', type=float, default=3.5)
    parser.add_argument('--k-count', type=int, default=51)
    args = parser.parse_args(argv)
    if not (0 < args.k_min < args.k_max and np.isfinite(args.k_max) and args.k_count >= 2):
        parser.error('require 0 < k-min < k-max and k-count >= 2')
    directions, weights = sphere_quadrature(args.n_theta, args.n_phi)
    mesh = prepare_mesh(args)
    ks = np.linspace(args.k_min, args.k_max, args.k_count)
    indicator = np.empty(len(ks))
    # Sampling point is the sphere center, so the RHS phase must translate with it.
    for i, k in enumerate(ks):
        F = assemble_far_field(k, directions, args, mesh)
        U, s, Vh = weighted_far_field_svd(F, weights, weights)
        rhs = np.sqrt(weights)*np.exp(-1j*k*(directions @ args.center))
        indicator[i] = np.linalg.norm(tikhonov(U, s, Vh, rhs, args.regularization))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    np.savez(args.output_dir/'ite_indicator.npz', k=ks, indicator=indicator,
             n=args.n, radius=args.radius, center=args.center, backend=args.backend,
             regularization=args.regularization, n_theta=args.n_theta, n_phi=args.n_phi)
    import matplotlib
    if not args.show:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.semilogy(ks, np.maximum(indicator, np.finfo(float).tiny))
    ax.set(xlabel='Wavenumber k', ylabel='Regularized solution norm',
           title=f'Sphere spectral indicator: n={args.n:g}, {args.backend}')
    ax.grid(True); fig.tight_layout()
    fig.savefig(args.output_dir/'ite_indicator.png', dpi=160)
    if args.show:
        plt.show()
    plt.close(fig)
    print(f'Saved {args.output_dir / "ite_indicator.npz"}')
    return indicator


if __name__ == '__main__':
    main()
