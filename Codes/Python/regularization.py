"""Complex Tikhonov regularization with area-weighted far-field operators."""
import numpy as np
from scipy.linalg import svd


def tikhonov(U, s, Vh, rhs, regularization):
    """Solve using SciPy's A = U @ diag(s) @ Vh, including complex adjoints."""
    if not np.isfinite(regularization) or regularization <= 0:
        raise ValueError("regularization must be positive and finite")
    q = len(s)
    projected = U[:, :q].conj().T @ rhs
    factors = s / (s**2 + regularization**2)
    if projected.ndim == 2:
        factors = factors[:, None]
    return Vh[:q, :].conj().T @ (factors * projected)


def weighted_far_field_svd(F, observation_weights, incident_weights):
    """SVD of sqrt(W_obs) F sqrt(W_inc); solution norm is the L2 norm."""
    F = np.asarray(F)
    wo, wi = np.asarray(observation_weights), np.asarray(incident_weights)
    if (F.ndim != 2 or wo.shape != (F.shape[0],) or wi.shape != (F.shape[1],)
            or not np.all(np.isfinite(F)) or not np.all(np.isfinite(wo))
            or not np.all(np.isfinite(wi)) or np.any(wo <= 0) or np.any(wi <= 0)):
        raise ValueError("finite matrix and matching positive quadrature weights required")
    return svd(np.sqrt(wo)[:, None] * F * np.sqrt(wi)[None, :], full_matrices=False)
