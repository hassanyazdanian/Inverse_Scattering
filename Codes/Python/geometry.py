"""Area quadrature and direction validation on the unit sphere."""
import numpy as np
from numpy.polynomial.legendre import leggauss


def sphere_quadrature(n_theta=6, n_phi=12):
    """Return (M, 3) directions and weights summing to 4*pi, without poles."""
    if int(n_theta) != n_theta or int(n_phi) != n_phi or n_theta < 2 or n_phi < 3:
        raise ValueError("n_theta >= 2 and n_phi >= 3 must be integers")
    z, wz = leggauss(n_theta)
    phi = np.arange(n_phi) * (2 * np.pi / n_phi)
    zz = np.repeat(z, n_phi)
    pp = np.tile(phi, n_theta)
    radius = np.sqrt(1 - zz**2)
    directions = np.column_stack((radius*np.cos(pp), radius*np.sin(pp), zz))
    return directions, np.repeat(wz, n_phi) * (2*np.pi/n_phi)


def unit_directions(directions):
    directions = np.atleast_2d(np.asarray(directions, dtype=float))
    if (directions.ndim != 2 or directions.shape[1] != 3
            or not np.all(np.isfinite(directions))
            or not np.allclose(np.linalg.norm(directions, axis=1), 1, atol=1e-10, rtol=0)):
        raise ValueError("directions must have shape (M, 3) and unit length")
    return directions


def sphere_parameters(k, n, radius, center=(0, 0, 0)):
    if not all(np.isfinite(v) and v > 0 for v in (k, n, radius)):
        raise ValueError("k, coefficient n, and radius must be positive and finite")
    center = np.asarray(center, dtype=float)
    if center.shape != (3,) or not np.all(np.isfinite(center)):
        raise ValueError("center must be a finite three-vector")
    return center
