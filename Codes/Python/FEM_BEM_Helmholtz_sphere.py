"""Plot a total-field slice using the optional legacy DOLFIN/Bempp backend."""
import argparse
from pathlib import Path
import numpy as np
from far_field_pattern import sphere_mesh, FEMBEMSolver


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=float, default=0.25)
    parser.add_argument('--k', type=float, default=6)
    parser.add_argument('--radius', type=float, default=1)
    parser.add_argument('--center', nargs=3, type=float, default=(0, 0, 0))
    parser.add_argument('--mesh-size', type=float, default=0.2)
    parser.add_argument('--grid-size', type=int, default=101)
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).resolve().parents[2]/'outputs')
    args = parser.parse_args(argv)
    if args.grid_size < 3:
        parser.error('grid-size must be at least 3')
    mesh = sphere_mesh(args.center, args.radius, args.mesh_size)
    solver = FEMBEMSolver(args.k, args.n, mesh)
    axis = np.linspace(-2, 2, args.grid_size)
    X, Y = np.meshgrid(axis, axis)
    points = np.column_stack((X.ravel(), Y.ravel(), np.zeros(X.size)))
    field = solver.total_field(np.ones(3)/np.sqrt(3), points, args.center, args.radius).reshape(X.shape)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    np.savez(args.output_dir/'sphere_total_field.npz', x=axis, y=axis, field=field,
             k=args.k, n=args.n, radius=args.radius, center=args.center, mesh_size=args.mesh_size,
             residual=solver.last_residual)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    picture = ax.pcolormesh(X, Y, field.real, shading='auto', cmap='RdBu_r')
    fig.colorbar(picture, ax=ax, label='Re(total field)')
    ax.set(xlabel='x', ylabel='y', aspect='equal', title=f'FEM-BEM z=0 slice: n={args.n:g}')
    fig.tight_layout(); fig.savefig(args.output_dir/'sphere_total_field.png', dpi=160)
    plt.close(fig)
    print(f'Saved {args.output_dir / "sphere_total_field.npz"}')


if __name__ == '__main__':
    main()
