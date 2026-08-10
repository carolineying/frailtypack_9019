# Small diagnostic result scope

These are paired smoke/diagnostic results, not the requested full simulation.
They use three existing true-variance-2 datasets/subsets, the first five
families, and at most 10 outer iterations. The subsets contain only 67--97
subjects and have maximum family sizes 17--28. Two of the three subsets are
too small to estimate the variance reliably and most fits did not converge.

All frailtypack rows report an integration accuracy failure because the
diagnostic deliberately uses the package's historical `1e-100` tolerances to
force use of the available rule budget. This should be read together with the
separate cap indicator, not as a nonpositive-integral failure. There were no
mode failures or Hessian regularizations.

The nominal 7,500 legacy setting evaluated as many as 46,377 nodes because
the historical high-dimensional `HRMSYM` starts a complete indivisible rule
while still below the cap. That behavior was retained for backward
compatibility. Adaptive integration pre-counts the next rule and never
exceeded 7,500 (maximum 6,087 here).

On the only subset where both 7,500-method fits converged, the estimates were
5.31 (legacy) and 1.28 (adaptive) for truth 2. Across all three deliberately
tiny subsets, adaptive RMSE was lower than legacy at each nominal setting,
but the convergence rate was only 1/3. These results are encouraging as a
numerical diagnostic but do not establish that the adaptive method improves
the intended large-family estimator. The full degree/scenario simulation is
therefore intentionally not launched.

For replicate 1, a deterministic row permutation was applied to the data and
to both axes of `K`. The absolute variance-estimate changes were 8.24 and 3.62
for legacy budgets 750 and 7,500, versus 0.0105 and 0.00280 for adaptive
budgets 750 and 7,500. This is one paired dataset, so it is evidence of better
ordering stability here rather than a general simulation conclusion.
