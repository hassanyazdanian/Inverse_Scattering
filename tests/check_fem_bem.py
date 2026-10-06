"""Optional integration/convergence check; run with environment-fem.yml.

This is deliberately separate from the lightweight default unittest suite.
"""
from pathlib import Path
import sys
import json
import numpy as np
from scipy.special import spherical_jn, spherical_yn, eval_legendre

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'Codes/Python'))
from geometry import sphere_quadrature
from sphere_scattering import far_field_matrix, scattering_coefficients
from far_field_pattern import sphere_mesh, FEMBEMSolver


def exact_field(points, direction, k, n):
    """Independent partial-wave total field for a centered unit sphere."""
    points = np.asarray(points)
    r = np.linalg.norm(points, axis=1)
    cosine = (points @ direction)/r
    result = np.zeros(len(r), complex)
    t = scattering_coefficients(k, n)
    inside = r < 1
    for m, coefficient in enumerate(t):
        radial = np.empty(len(r), complex)
        exterior_r = k*r[~inside]
        radial[~inside] = spherical_jn(m, exterior_r)+coefficient*(
            spherical_jn(m, exterior_r)+1j*spherical_yn(m, exterior_r))
        boundary = spherical_jn(m, k)+coefficient*(spherical_jn(m, k)+1j*spherical_yn(m, k))
        interior_coefficient = boundary/spherical_jn(m, k*np.sqrt(n))
        radial[inside] = interior_coefficient*spherical_jn(m, k*np.sqrt(n)*r[inside])
        result += (1j**m)*(2*m+1)*radial*eval_legendre(m, cosine)
    return result


def main():
    k, n = 1.0, 2.25
    observations, _ = sphere_quadrature(3, 6)
    incident = np.array([1., 0., 0.])
    exact = far_field_matrix(k, n, observations, incident)[:, 0]
    records = []
    for h in (0.45, 0.30, 0.20):
        mesh = sphere_mesh(mesh_size=h)
        solver = FEMBEMSolver(k, n, mesh)
        computed = solver.far_field(incident, observations)
        record = dict(mesh_size=h, vertices=mesh.num_vertices(), cells=mesh.num_cells(),
                      relative_far_field_error=float(np.linalg.norm(computed-exact)/np.linalg.norm(exact)),
                      gmres_residual=float(solver.last_residual))
        records.append(record)
        print(json.dumps(record), flush=True)
    assert records[-1]['relative_far_field_error'] < records[0]['relative_far_field_error']
    assert records[-1]['relative_far_field_error'] < 0.08
    assert max(row['gmres_residual'] for row in records) < 2e-7
    # Reuse the assembled operator for another incident direction.
    second = np.array([0., 0., 1.])
    reference = far_field_matrix(k, n, observations, second)[:, 0]
    assert np.linalg.norm(solver.far_field(second, observations)-reference)/np.linalg.norm(reference) < 0.08
    points = np.array([[0.2, 0.3, 0.1], [-0.3, 0.1, 0.2], [1.5, 0.2, 0.3], [-1.5, 0.3, 0.2]])
    reference = exact_field(points, incident, k, n)
    computed = solver.total_field(incident, points)
    field_error = float(np.linalg.norm(computed-reference)/np.linalg.norm(reference))
    assert field_error < 0.08
    report = dict(k=k, n=n, convergence=records, relative_total_field_error=field_error)
    output = ROOT/'outputs'; output.mkdir(exist_ok=True)
    (output/'fem_validation.json').write_text(json.dumps(report, indent=2)+'\n')
    print('PASS optional FEM-BEM far-field and total-field checks', field_error, flush=True)


if __name__ == '__main__':
    main()
