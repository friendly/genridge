[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/genridge)](http://cran.r-project.org/package=genridge)
[![](http://cranlogs.r-pkg.org/badges/grand-total/genridge)](https://cran.r-project.org/package=genridge)
[![Rdoc](http://www.rdocumentation.org/badges/version/genridge)](http://www.rdocumentation.org/packages/genridge)

# genridge

## Generalized Ridge Trace Plots for Ridge Regression

Version 0.6-6

### Overview

The genridge package introduces generalizations of the standard univariate
ridge trace plot used in ridge regression and related methods (Friendly, 2013)  These graphical methods
show both bias (actually, shrinkage) and precision, by plotting the covariance ellipsoids of the estimated
coefficients, rather than just the estimates themselves.  2D and 3D plotting methods are provided,
both in the space of the predictor variables and in the transformed space of the PCA/SVD of the
predictors.  

### Details

This package provides computational support for the graphical methods described in Friendly (2013). Ridge regression models may be fit using the function `ridge`, which incorporates features of `MASS::lm.ridge` and `ElemStatLearn::simple.ridge`. In particular, the shrinkage factors in ridge regression may be specified either in terms of the constant added to the diagonal of $X^T X$ matrix (lambda), or the equivalent number of degrees of freedom.

More importantly, the `ridge` function also calculates and returns the associated covariance matrices of each of the ridge estimates, allowing precision to be studied and displayed graphically.

This provides the support for the main plotting functions in the package:

* `plot.ridge`: Bivariate ridge trace plots

* `pairs.ridge`: All pairwise bivariate ridge trace plots

* `plot3d.ridge`: 3D ridge trace plots

* `traceplot`: Traditional univariate ridge trace plots

In addition, the function pca.ridge transforms the coefficients and covariance matrices of a ridge object from predictor space to the equivalent, but more interesting space of the PCA of $X^T X$ or the SVD of X. The main plotting functions also work for these objects, of class `c("ridge", "pcaridge")`.

Finally, the functions `precision` and `vif.ridge` provide other useful measures and plots.

## Installation

Get the released version from CRAN:

     install.packages("genridge")

The development version can be installed to your R library directly from this repo via:

     if (!require(devtools)) install.packages("devtools")
     library(devtools)
     install_github("friendly/genridge")


## References

Friendly, M. (2013).
The Generalized Ridge Trace Plot: Visualizing Bias *and* Precision.
*Journal of Computational and Graphical Statistics*, **22**(1), 50-68,
[doi](http://dx.doi.org/10.1080/10618600.2012.681237),
[genridge.pdf](http://euclid.psych.yorku.ca/datavis/papers/genridge.pdf)
