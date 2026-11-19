# Notes on lmridge (Imdad Ullah, Aslam & Altaf, 2018)

Source: `dev/lmridge-RJ-2018-060.pdf`. Notes only (2026-11-19); nothing implemented.

Imdad Ullah, M., Aslam, M. & Altaf, S. (2018). lmridge: A Comprehensive R Package for Ridge
Regression. *The R Journal*, 10(2), 326--343. doi:10.32614/RJ-2018-060 (verified with Crossref).

## 1. What the paper/package is

A "complete suite" of ridge **statistics**: fitting (`lmridge()`, several scalings), `summary()`
with SEs, t-tests and an F-test, `kest()` (about 22 published estimators of k), `rstats1()`
(variance, bias^2, MSE, sigma^2, F, R^2, adj-R^2, condition number), `rstats2()` (C_k, DF, EP,
REDF, Eft, ISRM, m-scale, PRESS), `hatr()`, `infocr()` (AIC/BIC), `press()`, and `predict()`,
`fitted()`, `residuals()`. Plots are simple line plots of each statistic against k (ridge trace,
VIF trace, bias/variance/MSE, CV/GCV, AIC/BIC, ISRM/m-scale, df/RSS/PRESS).

Its emphasis is the opposite of genridge's: many numbers and criteria for *choosing k*, little
visualization of the *coefficients' precision* (no ellipses, no PCA/biplot view).

## 2. Should it be cited? Yes, briefly

- It is the most complete R toolkit for univariate ridge statistics, and its Table 1 compares
  ridge software, including genridge (0.6-6, 2017). That table credits genridge only with
  "ridge trace" and "VIF"; it now does much more, so a pointer the other way is fair.
- Where:
  - **"Getting started" vignette**: a short "Related packages" paragraph (MASS::lm.ridge, glmnet,
    lmridge, ridge) saying genridge focuses on visualizing bias and precision, and lmridge is the
    place to go for the long list of k estimators, tests and ridge statistics.
  - Optionally `@seealso` in `?ridge`: `lmridge::lmridge()` (would need lmridge in Suggests, or a
    plain-text mention without a link; plain text avoids the dependency).
  - README "related packages", when the README is revised.
- **Bib entry** (not yet in `statistics.bib`; only `R-lmridge` package entries exist in the
  `Rpackages-*.bib` files). Proposed, for statistics.bib when you OK it:

```bibtex
@Article{ImdadUllah-etal:2018,
  author  = {Imdad Ullah, Muhammad and Aslam, Muhammad and Altaf, Saima},
  title   = {{lmridge}: A Comprehensive {R} Package for Ridge Regression},
  journal = {The R Journal},
  year    = {2018},
  volume  = {10},
  number  = {2},
  pages   = {326--343},
  doi     = {10.32614/RJ-2018-060},
}
```

## 3. Ideas that could be added to genridge (ranked by fit with genridge's aims)

1. **Plug-in bias^2 / variance / MSE trace** (their Table 2, items 6--7; `bias.plot()`).
   From the SVD genridge already has: Var = sigma^2 sum d_j^2/(d_j^2+k)^2,
   bias^2 = k^2 sum alpha_j^2/(d_j^2+k)^2, with OLS alpha and sigma^2 plugged in. Cheap to add
   as components of the `ridge` object (or in `precision()`), plus a plot against df.
   - Pairs naturally with the simulation in the MLM article, which shows the *true* tradeoff.
     An instructive comparison: plug-in MSE vs true MSE (the plug-in uses the unreliable OLS
     estimates, so it can be badly off under collinearity -- worth showing, not hiding).
   - Multivariate version is immediate: sum over responses, or per response.
2. **More selection criteria in `best()` / `criteria`**: PRESS (closed form from the hat-matrix
   diagonal, no refitting), C_p (Mallows/Kennard), AIC/BIC with df(k). All are cheap given the
   SVD. A single "criteria trace" plot (each criterion rescaled to [0,1] vs df, with minima
   marked) would be a better display than lmridge's separate panels.
3. **Show many k estimators on the trace plot**, rather than as a table: a rug or tick marks on
   the λ axis of `traceplot()` for a chosen set of estimators (HKB, LW, GCV now; optionally
   Kibria AM/GM/MED, Thisted, Dwividi-Srivastava ...). The *spread* of these values -- often
   two orders of magnitude, as in their Hald example (0.003 to 15) -- is itself a useful message.
   Don't port all 22; a handful with clear references is enough.
4. **Basic model methods**: `predict()`, `fitted()`, `residuals()` for `ridge` (and `mlmridge`),
   needed for PRESS/CV anyway. genridge has none at present.
5. **Condition-number trace**, CN_k = sqrt((d_1^2+k)/(d_p^2+k)): fits the collinearity theme;
   could go in `precision()` alongside the existing det/trace/max-eigenvalue measures.
6. **`summary()` with SEs**: SEs are already available from `cov`. I would *not* add t/F tests:
   the ridge estimator is biased, so these tests are questionable (the paper's own Hald example
   finds "significant" coefficients whose sign depends on k). If added, say so plainly.

Not worth borrowing: ISRM / m-scale (Vinod 1976) and the "effectiveness index" -- obscure, and
they add little beyond df and the VIF trace; the separate line plots for each statistic.

## 4. Data

They use `Hald` (cement data; in lmridge and mctest). genridge's `Acetylene` and `longley`
already cover this ground; Hald could be a further small example (4 predictors, n = 13, VIF up
to 282), but it is not needed.
