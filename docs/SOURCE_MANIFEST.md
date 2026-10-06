# Source provenance

Inputs: the public `hassanyazdanian/Inverse_Scattering` repository at commit
`1f1b6d713879131dc1bbe7cdf7a3f94c0c69e4f2`, the supplied `Archive.zip`, and
`Presentation_rv2.pdf`. The archive had 471 entries, including generated output,
duplicate experiments, notebooks, meshes, and third-party papers.

| Input | Cleaned destination / treatment |
| --- | --- |
| Repository `Codes/2D/main.m`, `Codes/3D/main_3D.m` | Retained entry points; shared numerical implementation in `Codes/common/` |
| Archive `2D/basis_functions.m` | Restored as `Codes/2D/basis_functions.m`; plotting/export block removed; respects `N` |
| Repository/archive `3D/basis_functions_3D_1.m` and archive `_2.m` | Callable `Codes/3D/basis_functions_3D_1.m`, with consistent Legendre factors |
| Archive `2D/A_function.m`, `B_function.m`, `C_function.m` and 3D counterparts | Full symbolic reference assembly via compatibility wrappers and `ite_symbolic_matrix.m` |
| Repository `A1/B1/C1` functions | Same public names, full reference assembly instead of hard-coded entries |
| Disk/ball determinant and Bessel utilities | Cleaned in their existing dimension folders; numeric origin limits and explicit root brackets |
| Archive `2D/EigenvalueFromFarfield/LSM/` | Consolidated into `Codes/2D/LSM/lsm_eigenvalues.m`, `lsm_shape.m`, and forward-solver helpers; original μ-diff notices retained |
| Archive `3D/linux/Inverse_Scattering/far_field_pattern.py` | Optional reusable `Codes/Python/far_field_pattern.py` |
| Archive `3D/linux/Inverse_Scattering/LSM_EigenValue.py` and `LSM_Penetrable.py` | Cleaned command-line workflows with shared quadrature and regularization |
| Archive `3D/linux/Inverse_Scattering/FEM_BEM_Helmholtz_sphere.py` | Cleaned total-field slice command; duplicate far-field variants consolidated |
| Archive `3D/linux/Inverse_Scattering/Spherical_domain_det_root.py` | Parameterized bracketed-root command |
| Archive `2D/Saved_results/k0_N_for2to20_3.mat` | Unchanged historical fixture `tests/data/archived_2d_sweep.mat` |
| Attached `Presentation_rv2.pdf` | Byte-identical to existing repository PDF; retained once at root |
| Five GitHub attachments in the original README | Downloaded unchanged to `docs/figures/historical_*.png` |

New work includes the normalized numerical Galerkin assembly, analytical sphere
benchmark, shared Python validation/helpers, tests, dependency specifications,
documentation, CI, and the two recomputed figures shown in the README.

## Deliberate omissions

Generated `.fig`/`.emf` files, local export directories, duplicate notebooks,
cached Python bytecode, generated mesh files, third-party research PDFs, and
`2D/EigenvalueFromFarfield/ChatGPT/` exploratory scripts are not included in the
release tree. No unsupported experimental script is presented as a validated solver.

The existing unused helpers `2D/d_sphbes.m`, `2D/sphbes.m`,
`2D/simplify_basis.m`, `3D/d_sp_besselj.m`, `3D/del_2_phi*.m`,
`3D/sp_laplacian.m`, `3D/simplify_basis_3D.m`, and `3D/polarplot3d.m`
were removed. The original versions remain available in git history and the
user's source archive. Their removal is reflected in the supplied git patch.

The complete source archive is not embedded in the cleaned package. Preserve
your original upload separately if you need the exploratory work or raw outputs.

## Presentation identity

`Presentation_rv2.pdf`: 41 pages, 3,314,283 bytes.

SHA-256:
`069f7dbbb44dee91a3653f1a20a7af2cc7b3008fabe6f2f266da492c5c185632`.
