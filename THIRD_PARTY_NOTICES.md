# Licensing and attribution

The original project license is retained in [LICENSE](LICENSE), copyright
2026 hassanyazdanian. Original project code and new project utilities are MIT
licensed, except for the expressly identified derived files below. The root
license does not replace third-party terms.

## μ-diff-derived MATLAB examples

`Codes/2D/LSM/lsm_eigenvalues.m` and `Codes/2D/LSM/lsm_shape.m` reorganize the
project's archived μ-diff examples. Their original attribution is preserved:
copyright 2014–2020 X. Antoine and B. Thierry, University of Lorraine, CNRS,
France. These files are distributed under **GPL-3.0-or-later**; the complete
license is in [LICENSE-mu-diff.txt](Codes/2D/LSM/LICENSE-mu-diff.txt).
The upstream README says version 2 or later, while its actual `Doc/LICENSE.txt`
says version 3 or later. This package follows the included version 3-or-later
notice for these derived scripts and preserves that license file unchanged.

The μ-diff toolbox itself is an optional external dependency and is not bundled.
The analytical disk forward solver and generic regularization helper are new
project utilities. Their use by the GPL examples does not remove the examples'
GPL terms. When distributing a combined derivative, preserve the applicable
licenses and source notices.

Upstream: <https://github.com/mu-diff/mu-diff>.

## Bempp FEM-BEM coupling

`Codes/Python/far_field_pattern.py` reorganizes Hassan Yazdanian's archived
FEM-BEM workflow. The block coupling follows the Bempp Helmholtz/FEniCS example.
The Bempp Team's MIT notice is retained in
[Bempp-cl-LICENSE.txt](docs/licenses/Bempp-cl-LICENSE.txt).

Upstream: <https://github.com/bempp/bempp-cl>;
documentation: <https://bempp.com/handbook/>.

## External dependencies and research material

NumPy, SciPy, Matplotlib, legacy FEniCS/DOLFIN, Bempp-cl, meshio, pygalmesh, and
μ-diff retain their own licenses. Their installed packages are not included here.
Third-party papers and unrelated utilities from the supplied archive are not
redistributed. The project presentation and historical project figures retain
their existing authorship; no new third-party paper license is asserted.
