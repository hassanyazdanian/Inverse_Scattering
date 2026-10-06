# Historical regression fixture

`archived_2d_sweep.mat` is an unchanged copy of
`2D/Saved_results/k0_N_for2to20_3.mat` in the supplied `Archive.zip`.
It stores `k0_N`, a 1×184 vector for `n=2:0.1:20.3`.

The values are retained to quantify changes from the historical implementation.
They are **not** treated as exact eigenvalues. Analytical determinant roots and
convergence tests are the accuracy checks. The test suite reports the difference
from this fixture without forcing the cleaned solver to reproduce older errors.
