# Validation report

Validated on Linux on 2026-10-06. All reported numbers below were measured on
the cleaned implementation. This report distinguishes discrete solver residuals,
errors against analytical solutions, and tests that remain unperformed.

## Executed checks

| Check | Outcome |
| --- | --- |
| Python standard suite | 17 tests passed on Python 3.12.14 and Python 3.11.16 |
| MATLAB-compatible numerical suite | Passed under GNU Octave 10.3.0 |
| Both MATLAB/Octave main scripts | Executed; recovered `n=3` and `n=40`; exported data and plots |
| Analytical sphere CLI examples | Root, spectral scan, and shape slice ran and produced finite data and PNG files |
| Disk LSM examples | Spectral scan and shape smoke tests passed |
| μ-diff optional connection | Analytical/multiple-scattering solver comparison passed for one disk |
| FEM-BEM optional backend | Three meshes, two incident directions, and interior/exterior total-field checks passed |
| Presentation | Attached and repository PDFs are byte-identical; 41 pages |

The Python suite checks complex rectangular/multiple-RHS Tikhonov solutions
against normal equations, weighted SVD reconstruction, sphere area and moments,
absence of duplicate directions, zero scattering at `n=1`, partial-wave boundary
conditions and lossless identities, translation/reciprocity/rotation invariance,
series convergence, root/radius scaling, invalid inputs, and all CLI help paths.

The Octave suite checks polynomial quadrature moments through degree 20,
disk/ball determinant roots, Bessel values at the origin, all 40 tabulated
clamped basis derivatives, normalized mass matrices, symmetry/positivity,
80-versus-160 point radial quadrature, polynomial residuals, both reconstruction
sweeps, increasing-basis convergence, invalid basis/zero-contrast guards,
complex Tikhonov solutions, disk scattering identities, and LSM smoke tests.

## Galerkin accuracy

Unit disk, `N=20`. Analytical references below are independently bracketed
order-zero characteristic roots for these examples.

| n | Archived sweep | Cleaned Galerkin k | Analytical k |
| ---: | ---: | ---: | ---: |
| 3 | 4.187932090867 | 4.167392441300 | 4.159235819878 |
| 5.5 | 2.356579253927 | 2.047333620204 | 2.046295617688 |
| 6 | 1.886419241678 | 1.850299940252 | 1.849719369141 |
| 10 | 1.296497455420 | 1.296496366387 | 1.296304915625 |
| 20 | 0.881641048491 | 0.881638385157 | 0.881543660339 |

Across the archived `n=2:0.1:20.3` grid, the largest relative change from the
historical output is **13.122649%**. The fixture is a historical comparison,
not an accuracy oracle. The new values are closer to the analytical roots in
the examples above; no claim is made that every parameter regime has been
fully validated. Conditioning, root filtering, and assembly all changed, so a
single cause is not assigned to every difference.

Unit ball at `n=40`, with analytical order-zero root `0.716709033375646`:

| N | Computed k | Absolute error |
| ---: | ---: | ---: |
| 8 | 0.717844989271 | 1.13596e-3 |
| 15 | 0.716982773619 | 2.73740e-4 |
| 20 | 0.716799377916 | 9.03445e-5 |

Normalized mass matrices agree with the identity within the suite's `1e-10`
tolerance; selected polynomial residuals are below `1e-10`. A small polynomial
residual does **not** imply a small finite-basis error, as the table illustrates.

### Angular branches matter

For a unit ball with `n=3`, the order-one determinant has a root near
`4.101812182726`, below the order-zero root `4.443358060871`. A script finding
the first order-zero root therefore must not label it the smallest eigenvalue
over every angular branch. The cleaned `N=8` Galerkin value is `4.233381027966`;
`N=20` gives `4.126465830835`. These are finite-basis approximations, not exact
roots. The code does not certify global completeness of the transmission spectrum.

## Analytical forward solvers and μ-diff

The disk comparison uses `k=2`, `n=3`, radius 1, and 17 incident/observation
angles. Relative Frobenius error against μ-diff is **8.956e-16**.
The checked μ-diff revision is
`6b7991424b570a3dbecfb9c3629d58b205f03e2e`.

The sphere analytical series satisfies boundary continuity and the lossless
partial-wave identity `Re(t_l)+|t_l|^2=0` to the test tolerances. Its truncation
was checked against a 50th-order expansion for the supplied `k=6`, `n=0.25`
example. These tests do not establish uniform accuracy at arbitrarily large
wavenumbers, coefficients, or radii.

LSM uses the exponential sampling RHS and the weighted matrix
`sqrt(W_obs) F sqrt(W_inc)`. Constant fundamental-solution prefactors are omitted
from the RHS; indicator magnitudes therefore should not be compared directly
to a publication using a different normalization or regularization convention.
A center-point spectral scan primarily probes the radial branch. Near an exact
singularity, finite regularization can produce a dip or split peak; the scripts
do not infer an eigenvalue automatically from the largest sampled norm.

## FEM-BEM convergence

The optional test uses a unit sphere, `k=1`, `n=2.25`, piecewise-linear FEM,
piecewise-constant boundary flux, and 18 far-field directions. Errors are
Euclidean relative errors against the analytical partial-wave solution.

| Mesh size bound | Vertices | Tetrahedra | Far-field error | GMRES relative residual |
| ---: | ---: | ---: | ---: | ---: |
| 0.45 | 60 | 175 | 14.4850% | 4.728e-9 |
| 0.30 | 176 | 651 | 6.1150% | 2.307e-9 |
| 0.20 | 504 | 2230 | 2.8496% | 4.763e-10 |

The finest mesh also passed a second incident-direction check. At four sample
points (two interior and two exterior), the relative total-field error was
**1.2416%**. The interior-field check caught and resolved a legacy DOLFIN array
layout issue: real and imaginary coefficient buffers must be contiguous.
Raw benchmark numbers are in [fem_validation.json](fem_validation.json).

This is a **low-frequency integration and convergence check**, not a production
accuracy guarantee. At larger `k`, refine the mesh and angular quadrature and
check solver conditioning. First-use just-in-time compilation can be slow.
Points exactly on the sphere boundary are masked in total-field plots; points
close to the boundary can be affected by the tetrahedral approximation and
potential quadrature. Mesh counts may vary across mesher/platform versions.

## Tested environments

| Component | Standard Python environment | Optional FEM environment |
| --- | --- | --- |
| Python | 3.12.14 | 3.11.16 |
| NumPy | 2.3.5 | 1.26.4 |
| SciPy | 1.17.0 | 1.13.1 |
| Matplotlib | 3.10.8 | 3.10.9 |
| DOLFIN | Not installed | 2019.1.0 |
| Bempp-cl | Not installed | 0.4.2 |
| Numba | Not required | 0.68.0 |
| meshio / pygalmesh | Not required | 5.3.5 / 0.10.7 |

The numerical `.m` checks used GNU Octave 10.3.0. The optional environment was
installed from conda-forge plus the Bempp-cl wheel. Harmless optional Gmsh,
legacy XML, and headless-graphics warnings occurred; no GUI or Gmsh dependency
is needed for the tested path.

## Not verified or not implemented

- MATLAB itself and its optional Symbolic Math Toolbox routines were not run.
  Octave coverage applies to the numerical core, not symbolic compatibility.
- GitHub Actions has been configured locally but has not yet executed on GitHub.
- Windows/macOS dependency installation, GPU/OpenCL execution, MPI parallel
  execution, and high-frequency FEM-BEM LSM sweeps were not exercised.
- No noisy measured data, limited-aperture experiment, heterogeneous material,
  non-spherical FEM geometry, uniqueness proof, or general-purpose inverse solver
  is supplied. Fixed regularization is a demonstrated choice, not an optimized
  parameter-selection rule.
- No claim is made that the corrected scripts exactly regenerate the original
  presentation figures. The presentation and historical images are retained as
  archival material with their original authorship.
