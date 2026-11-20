## genridge 0.8.1

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

### New
o New vignette, "Getting started with genridge", introducing the main functions with the `longley` and
  `Detroit` examples.
o `traceplot()` is now an S3 generic, with methods for `"ridge"` and `"mlmridge"` objects.
o `traceplot()` and the `plot()` method for `vif()` results now leave enough room for the variable labels, based on their width, so long names are no longer clipped at the side of the plot.

o New `ridge.lm()` method: `ridge()` can now be applied directly to a model fitted by `lm()`, e.g., `ridge(lm(y ~ x1 + x2, data), lambda = ...)`. The response, model matrix (with any factor contrasts and transformed terms) and offset are taken from the fitted model. `glm` and weighted `lm` fits give an error; a multivariate `lm()` fit gives a multivariate ridge regression (see above).

o Documentation converted to roxygen markdown (`Roxygen: list(markdown = TRUE)`): `\code{}`, `\link{}`, `\emph{}` and `\bold{}` replaced by markdown equivalents, and `\describe{}`, `\itemize{}`, `\enumerate{}` and `@return` `\item{}` lists by markdown lists. Also fixed a few small errors in the docs.

### Bug fixes
o `ridge()` with a formula and `svd = FALSE` failed with "variable lengths differ".
o `ridge()` with a formula failed when the model had only one predictor; `vif()` failed for such a model.
o The HKB and LW criteria, which have the factor (p - 2), are now `NA` for fewer than three predictors
  (they were 0 for p = 2 and negative for p = 1), and the plot methods no longer mark them.
o `best()` was documented but its generic was not exported; it now is.
o `best()` reported the same (wrong) `df` for all three criteria; each now gets its own equivalent degrees of freedom.

o Replaced links to `datavis.ca`, which is offline, with DOIs (documentation) and `friendly.github.io/papers` (README).

### Tests
o Added a testthat test suite: the multivariate fit against the univariate fits and the closed-form estimator, its
  covariance structure, the pooled GCV criteria, `best()`, the methods for `"mlmridge"` objects, and the univariate
  results against saved reference values from 0.8.1.

### Other
o genridge now requires R >= 4.1.0 (for the native pipe, `|>`).

## genridge 0.8.0 (2024-11-30)

This is a major update to the package adding additional plotting methods

o Fleshed out and installed `plot.precision()`
o Illustrate `plot.precision()` in `README.Rmd`
o Added `diab` data, diabetes from CASI
o Fixed documentation for `plot.ridge()` so that `plot.pcaridge()` is documented explicitly.
o `vif.ridge()` now returns a `"vif.ridge"` object list to prepare for a plot method.
o Added `plot.vif.ridge()` to plot VIFs vs shrinkage
o Added `norm.diff` = sqrt((b_OLS - b_ridge)^2) as a measure of shrinkage to `precision()`
o Added `plot.precision()` for plots of shrinkage vs. precision
o Added `best()` to display the optimal shrinkage criteria

## genridge 0.7.1 (2024-11-07)

o Added links to gentalk.pdf
o Improved README and fixed some examples
o `precision()` result gains a class "precision" in preparation for a plot method
o Implemented `plot.precision()` for plots shrinkage vs. variance, using various criteria
o Fixed warning from `ridge()` related to contrasts
o `ridge()` now collects the optimal `criteria` in a named list, which can be added in `plot.precision()`

## genridge 0.7.0 (2023-07-31)

o Converted the package to roxygen2 documentation, correcting some infelicities with S3 methods
o Added an extended README example.
o fix link to genridge paper PDF

## genridge 0.6.8 (2023-07-11)
o a maintenance release, correcting a problem with an \alias{} in one .Rd file

## genridge 0.6.7 (2020-01-07)
o Remove references to ElemStatLearn, now defunct

## genridge 0.6-6 (2017-10-01)
o bump package version, preparing to move to Github

## genridge 0.6-5 (2014-11-24)
o Use rgl:: in calls to rgl package for CRAN
o Fixed long lines in .Rd files

## genridge 0.6-4 (2012-2-10)
o Added contourf()

## genridge 0.6-3 (2011-12-28)
o Minor tweaks to plot3d.ridge
o Added biplot.ridge method for symmetry, showing the PCA vectors in variable space
o Added examples to ridge.Rd using data(manpower, package="bestglm")
o Added Manpower data

## genridge 0.6-2 (2011-12-19)
o Made ridge() an S3 generic with formula and default methods
o Fixed labels problems in biplot.pcaridge and plot3d.ridge
o Added Detroit data with examples showing most methods

## genridge 0.6-1 (2011-12-06)
o Fixed order of classes in pca.ridge
o Made plot3d an S3 generic; added plot3d.pcaridge and plot.pcaridge
o Added biplot.pcaridge() to draw variable vectors on a plot.pcaridge() plot
o Added options to precision.* for greater flexibility
o Added GCV computation to ridge()
o Added which.lambda to select ellipses in plot.ridge(), plot3d.ridge() etc.
o Added labels argument to plot.ridge() and fixed pairs.ridge for this.

## genridge 0.6-0 (2011-12-04)
o Added plot3d.ridge() using rgl, completing the main stages in visualization methods for ridge-like problems

## genridge 0.5-3 (2011-11-30)
o Added precision.ridge() method for measures of precision and shrinkage
o Added vif.ridge() method, extending the car S3 generic for ridge objects
o Now Depends: car

## genridge 0.5-2 (2011-11-25)
o Added Acetylene data
o Added vcov.ridge() method

## genridge 0.5-1 (2011-08-22)
o Fixed Authors@R in DESCRIPTION

## genridge 0.5-0 (2011-08-22)

o Added NAMESPACE
o Initial package version, released to R-Forge
o Added pca.ridge
