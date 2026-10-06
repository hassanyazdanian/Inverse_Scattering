"""Analytical penetrable sphere benchmark, n = c0**2 / c**2.

The scattered field is exp(i*k*r)/r * u_infinity + O(r**-2).
This independent reference supports the archived FEM/BEM and LSM workflows.
"""
import numpy as np
from scipy.special import spherical_jn, spherical_yn, eval_legendre
from geometry import sphere_parameters, unit_directions


def scattering_coefficients(k, n, radius=1.0, max_order=None):
    sphere_parameters(k, n, radius)
    if max_order is None:
        size = max(k*radius, k*radius*np.sqrt(n))
        max_order = max(12, int(np.ceil(size + 4*size**(1/3) + 8)))
    if int(max_order) != max_order or max_order < 0:
        raise ValueError("max_order must be a nonnegative integer")
    m = np.arange(max_order+1)
    x, y = k*radius, k*radius*np.sqrt(n)
    j, jp = spherical_jn(m, x), spherical_jn(m, x, derivative=True)
    ji, jip = spherical_jn(m, y), np.sqrt(n)*spherical_jn(m, y, derivative=True)
    h = j + 1j*spherical_yn(m, x)
    hp = jp + 1j*spherical_yn(m, x, derivative=True)
    return (j*jip - jp*ji) / (ji*hp - jip*h)


def far_field_matrix(k, n, observations, incidents, radius=1.0,
                     center=(0, 0, 0), max_order=None):
    """Matrix with observation rows and incident-direction columns."""
    center = sphere_parameters(k, n, radius, center)
    obs, inc = unit_directions(observations), unit_directions(incidents)
    cos_angle = np.clip(obs @ inc.T, -1, 1)
    coefficients = scattering_coefficients(k, n, radius, max_order)
    result = np.zeros(cos_angle.shape, dtype=complex)
    for m, value in enumerate(coefficients):
        result += (2*m+1)*value*eval_legendre(m, cos_angle)
    phase = np.exp(1j*k*((inc @ center)[None, :] - (obs @ center)[:, None]))
    return (-1j/k)*result*phase


def transmission_determinant(k, n, radius=1.0, order=0):
    """A root for one spherical harmonic order; not a global minimum search."""
    sphere_parameters(k, n, radius)
    if int(order) != order or order < 0:
        raise ValueError("order must be a nonnegative integer")
    x, y = np.asarray(k)*radius, np.asarray(k)*radius*np.sqrt(n)
    return (np.sqrt(n)*spherical_jn(order, x)*spherical_jn(order, y, derivative=True)
            - spherical_jn(order, y)*spherical_jn(order, x, derivative=True))
