## Test environments
* local Windows 11 x64 install, R version 4.6.1 (2026-06-24 ucrt)
* win-builder: R Under development (unstable) (2026-11-17 r90926 ucrt)

## R CMD check results
0 error(s) | 0 warning(s) | 0 note(s)

## Reverse dependencies checks

There are no reverse dependencies on CRAN.

## Comments

This release adds an experimental extension to multivariate ridge regression (a matrix of
responses), fixes several bugs, adds a vignette and a test suite, and now requires R >= 4.1.0.

### genridge 0.8.1

o **Experimental: multivariate ridge regression.** These functions are new in this version; their
  arguments and results may change in future releases.
    - `ridge()` now accepts a matrix response: `cbind(Y1, Y2) ~ ...` in a formula, a multivariate `lm()` fit,
      or a matrix `y`. All responses are fitted with common ridge constants, and the result has class `"mlmridge"`.
      The implementation was adapted from `lm_multiv_ridge()` in the **VARshrink** package.
    - The HKB, LW and GCV criteria are given for each response; GCV can also be pooled over responses
      (`pool = "scaled"`, `"raw"` or `"mahalanobis"`) to choose a common ridge constant.
    - `as.ridge()` extracts the results for one response as an ordinary `"ridge"` object.
    - `coef()`, `vcov()` and `print()` methods for `"mlmridge"` objects.
    - `best()` method for `"mlmridge"` objects: the optimal ridge constants for each response and criterion,
      plus GCV pooled over responses, as a data frame; `best(x, response = )` gives one response.
    - `traceplot()`, `plot()`, `pairs()`, `plot3d()`, `biplot()`, `pca()`, `precision()` and `vif()` methods
      apply the univariate methods to one response at a time (`response =`). `traceplot()` by default draws
      all responses, one panel each, marking the ridge constant chosen by pooled GCV (`GCV = TRUE`).
    - `plot(x, predictor = )` draws coefficient ellipses for one or more predictors in the space of two
      responses (`responses =`), traced over the ridge constants.
    - A one-column matrix response (e.g., `cbind(y) ~ ...`) gives an ordinary `"ridge"` result.

#### New
o New vignette, "Getting started with genridge", introducing the main functions with the `longley` and
  `Detroit` examples.
o `traceplot()` is now an S3 generic, with methods for `"ridge"` and `"mlmridge"` objects.
o `traceplot()` and the `plot()` method for `vif()` results now leave enough room for the variable labels, based on their width, so long names are no longer clipped at the side of the plot.

o New `ridge.lm()` method: `ridge()` can now be applied directly to a model fitted by `lm()`, e.g., `ridge(lm(y ~ x1 + x2, data), lambda = ...)`. The response, model matrix (with any factor contrasts and transformed terms) and offset are taken from the fitted model. `glm` and weighted `lm` fits give an error; a multivariate `lm()` fit gives a multivariate ridge regression (see above).

o Documentation converted to roxygen markdown (`Roxygen: list(markdown = TRUE)`): `\code{}`, `\link{}`, `\emph{}` and `\bold{}` replaced by markdown equivalents, and `\describe{}`, `\itemize{}`, `\enumerate{}` and `@return` `\item{}` lists by markdown lists. Also fixed a few small errors in the docs.

#### Bug fixes
o `ridge()` with a formula and `svd = FALSE` failed with "variable lengths differ".
o `ridge()` with a formula failed when the model had only one predictor; `vif()` failed for such a model.
o The HKB and LW criteria, which have the factor (p - 2), are now `NA` for fewer than three predictors
  (they were 0 for p = 2 and negative for p = 1), and the plot methods no longer mark them.
o `best()` was documented but its generic was not exported; it now is.
o `best()` reported the same (wrong) `df` for all three criteria; each now gets its own equivalent degrees of freedom.

o Replaced links to `datavis.ca`, which is offline, with DOIs (documentation) and `friendly.github.io/papers` (README).

#### Tests
o Added a testthat test suite: the multivariate fit against the univariate fits and the closed-form estimator, its
  covariance structure, the pooled GCV criteria, `best()`, the methods for `"mlmridge"` objects, and the univariate
  results against saved reference values from 0.8.1.

#### Other
o genridge now requires R >= 4.1.0 (for the native pipe, `|>`).
