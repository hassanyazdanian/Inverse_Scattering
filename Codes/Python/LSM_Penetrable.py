"""Linear sampling reconstruction of the z=0 cross-section of a sphere."""
import argparse
import numpy as np
from geometry import sphere_quadrature, sphere_parameters
from regularization import tikhonov, weighted_far_field_svd
from lsm_common import add_common_arguments, prepare_mesh, assemble_far_field


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_arguments(parser)
    parser.add_argument('--k', type=float, default=6)
    parser.add_argument('--grid-size', type=int, default=101)
    parser.add_argument('--extent', type=float, default=3)
    args = parser.parse_args(argv)
    sphere_parameters(args.k, args.n, args.radius, args.center)
    if args.grid_size < 3 or not np.isfinite(args.extent) or args.extent <= 0:
        parser.error('grid-size >= 3 and positive finite extent required')
    directions, weights = sphere_quadrature(args.n_theta, args.n_phi)
    mesh = prepare_mesh(args)
    F = assemble_far_field(args.k, directions, args, mesh)
    U, s, Vh = weighted_far_field_svd(F, weights, weights)
    axis = np.linspace(-args.extent, args.extent, args.grid_size)
    X, Y = np.meshgrid(axis, axis)
    points = np.column_stack((X.ravel(), Y.ravel(), np.zeros(X.size)))
    norms = np.empty(len(points))
    for start in range(0, len(points), 2048):
        sample = points[start:start+2048]
        rhs = np.sqrt(weights)[:, None]*np.exp(-1j*args.k*(directions @ sample.T))
        norms[start:start+len(sample)] = np.linalg.norm(
            tikhonov(U, s, Vh, rhs, args.regularization), axis=0)
    indicator = np.log10(np.maximum(norms, np.finfo(float).tiny)).reshape(X.shape)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    np.savez(args.output_dir/'lsm_shape.npz', x=axis, y=axis, indicator=indicator,
             k=args.k, n=args.n, radius=args.radius, center=args.center,
             backend=args.backend, regularization=args.regularization,
             n_theta=args.n_theta, n_phi=args.n_phi)
    import matplotlib
    if not args.show:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    picture = ax.pcolormesh(X, Y, indicator, shading='auto', cmap='viridis')
    fig.colorbar(picture, ax=ax, label='log10 regularized solution norm')
    if abs(args.center[2]) < args.radius:
        from matplotlib.patches import Circle
        circle = Circle(args.center[:2], np.sqrt(args.radius**2-args.center[2]**2),
                        fill=False, edgecolor='white', linestyle='--', label='True boundary')
        ax.add_patch(circle); ax.legend(loc='upper right')
    ax.set(xlabel='x', ylabel='y', aspect='equal', title=f'z=0 slice: n={args.n:g}, {args.backend}')
    fig.tight_layout(); fig.savefig(args.output_dir/'lsm_shape.png', dpi=160)
    if args.show:
        plt.show()
    plt.close(fig)
    print(f'Saved {args.output_dir / "lsm_shape.npz"}')
    return indicator


if __name__ == '__main__':
    main()
