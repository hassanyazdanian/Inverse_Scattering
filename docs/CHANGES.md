# Cleanup and numerical changes

Prepared against GitHub commit `1f1b6d713879131dc1bbe7cdf7a3f94c0c69e4f2`.
This is a functional research-code cleanup, not just a formatting change.

## Restored and integrated

- Restored the missing disk basis from the supplied `Archive.zip`.
- Converted `basis_functions_3D_1.m` from a plotting script into a callable,
  side-effect-free reference function returning exactly `N` modes.
- Consolidated the archived Python sphere workflows into reusable geometry,
  forward solver, regularization, and command-line modules.
- Recovered disk far-field/LSM workflows with an optional μ-diff backend.
  Analytical disk and sphere solvers were added as independent benchmarks and
  lightweight defaults; these analytical Python utilities are new additions.

## Numerical corrections

1. The Galerkin solver now assembles all mode couplings with normalized basis
   functions, exact angular orthogonality, and Gauss-Legendre radial quadrature.
   It replaces hand-selected matrix entries and repeated symbolic integration
   in the main workflows. Optional full symbolic assembly remains for reference.
2. Three-dimensional angular factors consistently use Legendre polynomials.
   The archived `_1`/`_2` scripts mixed Legendre and cosine factors in later modes.
3. The radial interval covers the full unit domain. Gaussian nodes avoid direct
   evaluation of removable origin singularities; a hole `r<0.02` is no longer
   cut out of the ball. Numeric spherical Bessel helpers have finite origin limits.
4. Eigenvalue filtering accepts small relative imaginary roundoff, rejects
   nonpositive/nonfinite values before taking square roots, and reports polynomial
   residuals. Tiny residuals certify the discrete equation, not continuum accuracy.
5. Python Tikhonov regularization uses **conjugate transposes**. SciPy returns
   `Vh`, not MATLAB's `V`; the archived `.T` implementation was incorrect for
   complex far-field matrices.
6. Sphere directions use Gauss-Legendre quadrature in `cos(theta)` and uniform
   azimuth. There are no repeated poles, and weights sum to `4*pi`. Both LSM
   implementations use the corresponding weighted discrete operator.
7. The coefficient convention is consistently `k_inside = k*sqrt(n)` and
   `Delta(u) + k^2*n*u = 0`. Conflicting `n`/`n^2` labels from the archive were
   resolved. The Python shape default represents `n=0.25`.
8. Shape reconstruction evaluates the requested `z=0` plane directly in chunks.
   The old Python workflow allocated and evaluated a full volume before plotting
   a slice. Translation phases and plotted boundary centers are retained correctly.
9. The FEM-BEM solver assembles once for a fixed wavenumber, supports the
   `bempp_cl.api` namespace, checks GMRES status and residuals, and creates meshes
   in a temporary directory without shipping stale mesh files.

## Portability and packaging

Windows export paths, mandatory plotting side effects, `clear all`, obsolete
exploratory utilities, and unused third-party plotting code were removed.
Generated output goes to `outputs/` and is excluded from git. Default workflows
no longer need the missing `csvd`, `discrep`, or external `tikhonov` functions.
The cleaned LSM examples use explicit fixed regularization; automatic Morozov
parameter selection from the exploratory scripts was not reimplemented.

The repository now includes installation instructions, separate optional FEM
dependencies, numerical tests, CI, source provenance, and scoped third-party
license notices. The presentation was already present and is unchanged.

## Interpretation of older results

The saved 2D sweep differs from the cleaned solver by up to **13.122649%** on the
original candidate grid. It remains a historical comparison fixture, not a
golden accuracy target. Selected comparisons against independent characteristic
roots are in [VALIDATION.md](VALIDATION.md). No single cause is asserted for
every historical discrepancy, and the previous published figures are not claimed
to be exact reproductions of this revised implementation.
