# Apply the prepared package

The package contains a complete cleaned repository tree and a binary-safe
`github-update.patch`. Nothing has been pushed to GitHub.

## Recommended: apply the patch to a clean clone

Keep any existing local work safe. In a separate directory:

```bash
git clone https://github.com/hassanyazdanian/Inverse_Scattering.git
cd Inverse_Scattering
git switch -c cleanup/github-ready
git apply --check /path/to/github-update.patch
git apply /path/to/github-update.patch
git diff --check
git status --short
```

The patch was prepared against commit
`1f1b6d713879131dc1bbe7cdf7a3f94c0c69e4f2`. If the repository has changed,
`git apply --check` may reject it. Review the conflicts or compare with the
complete cleaned tree; do not force an overwrite of newer work.

Run the checks described in the README, review the changes and license notices,
then commit and push the branch:

```bash
git add -A
git commit -m "Clean inverse scattering workflows and add numerical validation"
git push -u origin cleanup/github-ready
```

Open a pull request from this branch to `main`. GitHub Actions should then run
the core Python and Octave suites. GitHub authentication uses your normal local
credentials; the package contains no credentials or deployment configuration.

## Using the complete source tree

The `Inverse_Scattering/` folder in the package can be opened and run directly.
It contains no `.git` directory. Prefer the patch when updating the existing
repository: it also removes obsolete files, whereas simply copying new files
on top of an old checkout can leave obsolete files behind.

`outputs/` is created on demand and excluded from git. Only selected documented
figures are committed under `docs/figures/`; local environments and cached data
are not part of the package.
