# Package index

## Ridge functions

Functions for computing ridge regressions

- [`genridge-package`](https://friendly.github.io/genridge/reference/genridge-package.md)
  [`genridge`](https://friendly.github.io/genridge/reference/genridge-package.md)
  [`_PACKAGE`](https://friendly.github.io/genridge/reference/genridge-package.md)
  : Generalized ridge trace plots for ridge regression
- [`precision()`](https://friendly.github.io/genridge/reference/precision.md)
  : Measures of Precision and Shrinkage for Ridge Regression
- [`ridge()`](https://friendly.github.io/genridge/reference/ridge.md)
  [`coef(`*`<ridge>`*`)`](https://friendly.github.io/genridge/reference/ridge.md)
  [`print(`*`<ridge>`*`)`](https://friendly.github.io/genridge/reference/ridge.md)
  [`vcov(`*`<ridge>`*`)`](https://friendly.github.io/genridge/reference/ridge.md)
  [`best()`](https://friendly.github.io/genridge/reference/ridge.md) :
  Ridge Regression Estimates
- [`vif(`*`<ridge>`*`)`](https://friendly.github.io/genridge/reference/vif.ridge.md)
  [`print(`*`<vif.ridge>`*`)`](https://friendly.github.io/genridge/reference/vif.ridge.md)
  [`plot(`*`<vif.ridge>`*`)`](https://friendly.github.io/genridge/reference/vif.ridge.md)
  : Variance Inflation Factors for Ridge Regression

## Plot methods

Plotting ridge regression results

- [`contourf()`](https://friendly.github.io/genridge/reference/contourf.md)
  : Enhanced Contour Plots
- [`plot(`*`<ridge>`*`)`](https://friendly.github.io/genridge/reference/plot.ridge.md)
  [`plot(`*`<pcaridge>`*`)`](https://friendly.github.io/genridge/reference/plot.ridge.md)
  : Bivariate Ridge Trace Plots
- [`plot(`*`<precision>`*`)`](https://friendly.github.io/genridge/reference/plot.precision.md)
  : Plot Shrinkage vs. Variance for Ridge Precision
- [`plot3d()`](https://friendly.github.io/genridge/reference/plot3d.md)
  : 3D Ridge Trace Plots
- [`pairs(`*`<ridge>`*`)`](https://friendly.github.io/genridge/reference/pairs.ridge.md)
  : Scatterplot Matrix of Bivariate Ridge Trace Plots
- [`traceplot()`](https://friendly.github.io/genridge/reference/traceplot.md)
  : Univariate Ridge Trace Plots
- [`trans.colors()`](https://friendly.github.io/genridge/reference/trans.colors.md)
  : Make Colors Transparent

## PCA methods

Low-rank PCA methods

- [`pca()`](https://friendly.github.io/genridge/reference/pca.md) :
  Transform Ridge Estimates to PCA Space
- [`biplot(`*`<pcaridge>`*`)`](https://friendly.github.io/genridge/reference/biplot.pcaridge.md)
  : Biplot of Ridge Regression Trace Plot in SVD Space

## Multivariate ridge regression (experimental)

Ridge regression for several responses with common ridge constants:
[`ridge()`](https://friendly.github.io/genridge/reference/ridge.md) with
a matrix response returns an `"mlmridge"` object.

- [`as.ridge()`](https://friendly.github.io/genridge/reference/as.ridge.md)
  : Extract a Univariate Ridge Regression from a Multivariate One
- [`coef(`*`<mlmridge>`*`)`](https://friendly.github.io/genridge/reference/mlmridge-methods.md)
  [`vcov(`*`<mlmridge>`*`)`](https://friendly.github.io/genridge/reference/mlmridge-methods.md)
  [`print(`*`<mlmridge>`*`)`](https://friendly.github.io/genridge/reference/mlmridge-methods.md)
  : Methods for Multivariate Ridge Regression Objects
- [`traceplot(`*`<mlmridge>`*`)`](https://friendly.github.io/genridge/reference/mlmridge-plots.md)
  [`plot(`*`<mlmridge>`*`)`](https://friendly.github.io/genridge/reference/mlmridge-plots.md)
  [`pairs(`*`<mlmridge>`*`)`](https://friendly.github.io/genridge/reference/mlmridge-plots.md)
  [`plot3d(`*`<mlmridge>`*`)`](https://friendly.github.io/genridge/reference/mlmridge-plots.md)
  [`biplot(`*`<mlmridge>`*`)`](https://friendly.github.io/genridge/reference/mlmridge-plots.md)
  [`pca(`*`<mlmridge>`*`)`](https://friendly.github.io/genridge/reference/mlmridge-plots.md)
  [`precision(`*`<mlmridge>`*`)`](https://friendly.github.io/genridge/reference/mlmridge-plots.md)
  [`vif(`*`<mlmridge>`*`)`](https://friendly.github.io/genridge/reference/mlmridge-plots.md)
  : Plot and Summary Methods for Multivariate Ridge Regression

## Data sets

Data sets in the package

- [`Acetylene`](https://friendly.github.io/genridge/reference/Acetylene.md)
  : Acetylene Data
- [`Detroit`](https://friendly.github.io/genridge/reference/Detroit.md)
  : Detroit Homicide Data for 1961-1973
- [`Manpower`](https://friendly.github.io/genridge/reference/Manpower.md)
  : Hospital manpower data
- [`prostate`](https://friendly.github.io/genridge/reference/prostate.md)
  : Prostate Cancer Data
- [`diab`](https://friendly.github.io/genridge/reference/diab.md) :
  Diabetes Progression
