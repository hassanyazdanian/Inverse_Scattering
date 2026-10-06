"""Independent numerical invariants and command-line smoke tests."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import svd
from scipy.optimize import brentq
from scipy.special import spherical_jn, spherical_yn

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'Codes/Python'))
from geometry import sphere_quadrature, unit_directions
from regularization import tikhonov, weighted_far_field_svd
from sphere_scattering import scattering_coefficients, far_field_matrix, transmission_determinant


class RegularizationTests(unittest.TestCase):
    def test_complex_normal_equations(self):
        rng = np.random.default_rng(487)
        for shape in ((9, 5), (5, 9), (6, 6)):
            A = rng.normal(size=shape)+1j*rng.normal(size=shape)
            for nrhs in (1, 3):
                rhs = rng.normal(size=(shape[0], nrhs))+1j*rng.normal(size=(shape[0], nrhs))
                if nrhs == 1:
                    rhs = rhs[:, 0]
                expected = np.linalg.solve(A.conj().T@A+0.07**2*np.eye(shape[1]), A.conj().T@rhs)
                actual = tikhonov(*svd(A, full_matrices=False), rhs, 0.07)
                assert_allclose(actual, expected, rtol=2e-10, atol=2e-11)

    def test_weighted_equation(self):
        F = np.array([[1+2j, 3], [4j, 2-1j], [0.5, 2]])
        wo, wi = np.array([0.3, 0.8, 0.4]), np.array([0.2, 0.9])
        U, s, Vh = weighted_far_field_svd(F, wo, wi)
        assert_allclose((U*s)@Vh, np.sqrt(wo)[:, None]*F*np.sqrt(wi)[None, :], atol=1e-14)

    def test_invalid_regularization_and_weights(self):
        for lam in (0, -1, np.nan, np.inf):
            with self.assertRaises(ValueError):
                tikhonov(np.eye(2), np.ones(2), np.eye(2), np.ones(2), lam)
        with self.assertRaises(ValueError):
            weighted_far_field_svd(np.eye(2), [1, 0], [1, 1])


class SphereTests(unittest.TestCase):
    def setUp(self):
        self.directions, self.weights = sphere_quadrature(8, 16)

    def test_area_and_moments(self):
        assert_allclose(np.linalg.norm(self.directions, axis=1), 1, atol=1e-14)
        assert_allclose(self.weights.sum(), 4*np.pi, atol=1e-13)
        assert_allclose(self.weights@self.directions, 0, atol=1e-14)
        assert_allclose(self.directions.T@(self.weights[:, None]*self.directions),
                        np.eye(3)*4*np.pi/3, atol=1e-13)
        self.assertEqual(len(np.unique(self.directions, axis=0)), len(self.directions))

    def test_invalid_geometry(self):
        for args in ((1, 12), (6, 2), (3.5, 6)):
            with self.assertRaises(ValueError):
                sphere_quadrature(*args)
        with self.assertRaises(ValueError):
            unit_directions([[1, 1, 1]])

    def test_no_contrast(self):
        F = far_field_matrix(2, 1, self.directions, self.directions)
        assert_allclose(F, 0, atol=1e-14)

    def test_partial_wave_boundary_conditions(self):
        k, n = 1.3, 2.25
        t = scattering_coefficients(k, n, max_order=7)
        m = np.arange(len(t))
        j = spherical_jn(m, k); jp = spherical_jn(m, k, derivative=True)
        h = j+1j*spherical_yn(m, k)
        hp = jp+1j*spherical_yn(m, k, derivative=True)
        ji = spherical_jn(m, k*np.sqrt(n))
        jip = np.sqrt(n)*spherical_jn(m, k*np.sqrt(n), derivative=True)
        interior = (j+t*h)/ji
        assert_allclose(interior*jip, jp+t*hp, atol=1e-13)

    def test_lossless_partial_wave_identity(self):
        t = scattering_coefficients(6, 0.25)
        assert_allclose(t.real+abs(t)**2, 0, atol=1e-13)

    def test_translation_and_reciprocity(self):
        obs, inc = self.directions[:9], self.directions[30:41]
        c = np.array([0.3, -0.2, 0.4]); k = 2.3
        F = far_field_matrix(k, 3, obs, inc)
        shifted = far_field_matrix(k, 3, obs, inc, center=c)
        phase = np.exp(1j*k*((inc@c)[None, :]-(obs@c)[:, None]))
        assert_allclose(shifted, F*phase, atol=1e-14)
        assert_allclose(shifted, far_field_matrix(k, 3, -inc, -obs, center=c).T, atol=1e-14)

    def test_rotation_invariance(self):
        rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        d = self.directions[:12]
        assert_allclose(far_field_matrix(3, 4, d, d),
                        far_field_matrix(3, 4, d@rotation, d@rotation), atol=1e-14)

    def test_series_convergence(self):
        d = self.directions[::8]
        assert_allclose(far_field_matrix(6, 0.25, d, d),
                        far_field_matrix(6, 0.25, d, d, max_order=50), atol=1e-12)

    def test_transmission_root_and_radius_scaling(self):
        root = brentq(lambda k: transmission_determinant(k, 40), 0.70, 0.73)
        self.assertAlmostEqual(root, 0.716709033375646, places=12)
        half = brentq(lambda k: transmission_determinant(k, 40, radius=2), 0.35, 0.365)
        self.assertAlmostEqual(half, root/2, places=12)

    def test_invalid_physical_parameters(self):
        for k, n, r in ((0, 2, 1), (1, -1, 1), (1, 2, 0), (np.nan, 2, 1)):
            with self.assertRaises(ValueError):
                scattering_coefficients(k, n, r)


class CommandLineTests(unittest.TestCase):
    def run_script(self, name, *args):
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', MPLBACKEND='Agg')
        return subprocess.run([sys.executable, str(ROOT/'Codes/Python'/name), *args],
                              check=True, capture_output=True, text=True, timeout=90, env=env)

    def test_all_entry_points_help_without_optional_dependencies(self):
        for script in ('LSM_EigenValue.py', 'LSM_Penetrable.py',
                       'Spherical_domain_det_root.py', 'FEM_BEM_Helmholtz_sphere.py'):
            self.assertIn('usage:', self.run_script(script, '--help').stdout)

    def test_root_cli(self):
        self.assertIn('0.716709033375', self.run_script('Spherical_domain_det_root.py').stdout)

    def test_eigenvalue_scan_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.run_script('LSM_EigenValue.py', '--k-count', '5', '--n-theta', '3',
                            '--n-phi', '6', '--output-dir', tmp)
            result = np.load(Path(tmp)/'ite_indicator.npz')
            self.assertEqual(result['indicator'].shape, (5,))
            self.assertTrue(np.all(np.isfinite(result['indicator'])))
            self.assertGreater((Path(tmp)/'ite_indicator.png').stat().st_size, 1000)

    def test_shape_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.run_script('LSM_Penetrable.py', '--grid-size', '15', '--n-theta', '3',
                            '--n-phi', '6', '--output-dir', tmp)
            result = np.load(Path(tmp)/'lsm_shape.npz')
            self.assertEqual(result['indicator'].shape, (15, 15))
            self.assertTrue(np.all(np.isfinite(result['indicator'])))
            self.assertGreater((Path(tmp)/'lsm_shape.png').stat().st_size, 1000)


if __name__ == '__main__':
    unittest.main()
