"""Legacy DOLFIN/Bempp FEM-BEM forward solver, cleaned from the 2024 archive.

Author of the archived workflow: Hassan Yazdanian.
The coupling follows the Bempp Helmholtz/FEniCS example; see THIRD_PARTY_NOTICES.
Optional dependencies are loaded only when this backend is selected.
"""
from importlib import import_module
from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres
from geometry import sphere_parameters, unit_directions


def _backend():
    try:
        dolfin = import_module('dolfin')
        try:
            api = import_module('bempp_cl.api')
            prefix = 'bempp_cl.api'
        except ModuleNotFoundError as exc:
            if exc.name not in ('bempp_cl', 'bempp_cl.api'):
                raise
            api = import_module('bempp.api')
            prefix = 'bempp.api'
        fenics = import_module(prefix+'.external.fenics')
        blocked = import_module(prefix+'.assembly.blocked_operator')
        discrete = import_module(prefix+'.assembly.discrete_boundary_operator')
    except ImportError as exc:
        raise ImportError('FEM-BEM requires legacy DOLFIN and Bempp-cl. '
                          'Use environment-fem.yml; FEniCSx is not a drop-in replacement.') from exc
    return dolfin, api, fenics, blocked, discrete


def sphere_mesh(center=(0, 0, 0), radius=1, mesh_size=0.2):
    center = sphere_parameters(1, 1, radius, center)
    if not np.isfinite(mesh_size) or mesh_size <= 0:
        raise ValueError('mesh_size must be positive and finite')
    dolfin, *_ = _backend()
    import pygalmesh
    import meshio
    raw = pygalmesh.generate_mesh(pygalmesh.Ball(center, radius),
                                 max_cell_circumradius=mesh_size, verbose=False)
    tetra = raw.get_cells_type('tetra')
    if not len(tetra):
        raise RuntimeError('Mesher produced no tetrahedra')
    with TemporaryDirectory(prefix='inverse-scattering-') as tmp:
        path = Path(tmp)/'sphere.xml'
        meshio.write(path, meshio.Mesh(raw.points, [('tetra', tetra)]))
        return dolfin.Mesh(str(path))


class FEMBEMSolver:
    """Assemble once per k; solve for multiple incident waves on one mesh."""

    def __init__(self, k, n, mesh, tolerance=1e-8):
        sphere_parameters(k, n, 1)
        if not np.isfinite(tolerance) or tolerance <= 0:
            raise ValueError('tolerance must be positive and finite')
        dolfin, api, fenics, blocked, discrete = _backend()
        self.dolfin, self.api = dolfin, api
        self.k, self.n, self.tolerance = k, n, tolerance
        self.space = dolfin.FunctionSpace(mesh, 'CG', 1)
        self.trace_space, self.trace = fenics.fenics_to_bempp_trace_data(self.space)
        self.bem_space = api.function_space(self.trace_space.grid, 'DP', 0)
        identity = api.operators.boundary.sparse.identity(
            self.trace_space, self.bem_space, self.bem_space)
        mass = api.operators.boundary.sparse.identity(
            self.bem_space, self.bem_space, self.trace_space)
        double = api.operators.boundary.helmholtz.double_layer(
            self.trace_space, self.bem_space, self.bem_space, k, device_interface='numba')
        single = api.operators.boundary.helmholtz.single_layer(
            self.bem_space, self.bem_space, self.bem_space, k, device_interface='numba')
        u, v = dolfin.TrialFunction(self.space), dolfin.TestFunction(self.space)
        fem = fenics.FenicsOperator((dolfin.inner(dolfin.grad(u), dolfin.grad(v))
                                    - k**2*n*u*v)*dolfin.dx).weak_form()
        trace_op = LinearOperator(self.trace.shape, matvec=lambda x: self.trace @ x,
                                  matmat=lambda x: self.trace @ x, dtype=float)
        blocks = np.empty((2, 2), dtype=object)
        blocks[0, 0] = fem
        blocks[0, 1] = -self.trace.T @ mass.weak_form().to_sparse()
        blocks[1, 0] = (0.5*identity-double).weak_form()*trace_op
        blocks[1, 1] = single.weak_form()
        self.operator = blocked.BlockedDiscreteOperator(blocks)
        inverse = discrete.InverseSparseDiscreteBoundaryOperator
        p1 = inverse(fem.to_sparse().tocsc())
        p2 = inverse(api.operators.boundary.sparse.identity(
            self.bem_space, self.bem_space, self.bem_space).weak_form())
        split = self.space.dim()
        self.preconditioner = LinearOperator(self.operator.shape, dtype=complex,
            matvec=lambda x: np.concatenate((p1 @ x[:split], p2 @ x[split:])))

    def solve(self, incident):
        direction = unit_directions(incident)[0].copy()
        k, api = self.k, self.api

        @api.complex_callable
        def incident_trace(x, normal, domain_index, result):
            result[0] = np.exp(1j*k*np.dot(x, direction))

        wave = api.GridFunction(self.bem_space, fun=incident_trace)
        split = self.space.dim()
        rhs = np.concatenate((np.zeros(split), wave.projections(self.bem_space)))
        solution, info = gmres(self.operator, rhs, M=self.preconditioner,
                               rtol=self.tolerance, atol=0, restart=100, maxiter=1000)
        residual = np.linalg.norm(self.operator @ solution-rhs)/np.linalg.norm(rhs)
        if info != 0 or not np.isfinite(residual) or residual > 20*self.tolerance:
            raise RuntimeError(f'GMRES failed: info={info}, relative residual={residual:.3e}')
        dirichlet = api.GridFunction(self.trace_space, coefficients=self.trace @ solution[:split])
        neumann = api.GridFunction(self.bem_space, coefficients=solution[split:])
        self.last_residual = residual
        return solution[:split], dirichlet, neumann

    def far_field(self, incident, observations):
        observations = unit_directions(observations)
        _, dirichlet, neumann = self.solve(incident)
        operators = self.api.operators.far_field.helmholtz
        double = operators.double_layer(self.trace_space, observations.T, self.k,
                                        device_interface='numba')
        single = operators.single_layer(self.bem_space, observations.T, self.k,
                                        device_interface='numba')
        return np.asarray(double.evaluate(dirichlet)-single.evaluate(neumann)).ravel()

    def total_field(self, incident, points, center=(0, 0, 0), radius=1):
        """Evaluate total field; boundary points are excluded to avoid singular kernels."""
        center = sphere_parameters(self.k, self.n, radius, center)
        points = np.asarray(points, dtype=float)
        direction = unit_directions(incident)[0]
        if points.ndim != 2 or points.shape[1] != 3 or not np.all(np.isfinite(points)):
            raise ValueError('points must have shape (M, 3) and be finite')
        coefficients, dirichlet, neumann = self.solve(direction)
        distance = np.linalg.norm(points-center, axis=1)
        outside = distance > radius*(1+1e-8)
        inside = distance < radius*(1-1e-8)
        result = np.full(len(points), np.nan+1j*np.nan)
        if np.any(outside):
            operators = self.api.operators.potential.helmholtz
            double = operators.double_layer(self.trace_space, points[outside].T, self.k,
                                            device_interface='numba')
            single = operators.single_layer(self.bem_space, points[outside].T, self.k,
                                            device_interface='numba')
            result[outside] = (double.evaluate(dirichlet)-single.evaluate(neumann)).ravel()
            result[outside] += np.exp(1j*self.k*(points[outside] @ direction))
        if np.any(inside):
            real, imag = self.dolfin.Function(self.space), self.dolfin.Function(self.space)
            # Legacy DOLFIN's vector setter needs contiguous real-valued buffers.
            real.vector()[:] = np.ascontiguousarray(coefficients.real)
            imag.vector()[:] = np.ascontiguousarray(coefficients.imag)
            # Curved sphere and tetrahedral boundary differ near the surface.
            real.set_allow_extrapolation(True); imag.set_allow_extrapolation(True)
            result[inside] = [real(self.dolfin.Point(*p))+1j*imag(self.dolfin.Point(*p))
                              for p in points[inside]]
        return result


def far_field_pattern(d, k, x_hat, nn, mesh):
    """Compatibility with the archived function; x_hat has shape (3, M)."""
    return FEMBEMSolver(k, nn, mesh).far_field(d, np.asarray(x_hat).T)
