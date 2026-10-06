"""Shared command-line settings and optional forward-solver selection."""
from pathlib import Path
import numpy as np
from geometry import sphere_parameters
from sphere_scattering import far_field_matrix


def add_common_arguments(parser, default_n=0.25):
    parser.add_argument('--backend', choices=('analytic', 'fem-bem'), default='analytic')
    parser.add_argument('--n', type=float, default=default_n, help='coefficient n=c0^2/c^2')
    parser.add_argument('--radius', type=float, default=1)
    parser.add_argument('--center', nargs=3, type=float, default=(0, 0, 0))
    parser.add_argument('--n-theta', type=int, default=6)
    parser.add_argument('--n-phi', type=int, default=12)
    parser.add_argument('--mesh-size', type=float, default=0.2)
    parser.add_argument('--regularization', type=float, default=1e-6)
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).resolve().parents[2]/'outputs')
    parser.add_argument('--show', action='store_true')


def prepare_mesh(args):
    sphere_parameters(1, args.n, args.radius, args.center)
    if not np.isfinite(args.regularization) or args.regularization <= 0:
        raise ValueError('regularization must be positive and finite')
    if args.backend == 'fem-bem':
        from far_field_pattern import sphere_mesh
        return sphere_mesh(args.center, args.radius, args.mesh_size)
    return None


def assemble_far_field(k, directions, args, mesh=None):
    if args.backend == 'analytic':
        return far_field_matrix(k, args.n, directions, directions, args.radius, args.center)
    from far_field_pattern import FEMBEMSolver
    solver = FEMBEMSolver(k, args.n, mesh)
    return np.column_stack([solver.far_field(d, directions) for d in directions])
