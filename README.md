[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/genridge)](http://cran.r-project.org/package=genridge)
[![](http://cranlogs.r-pkg.org/badges/grand-total/genridge)](https://cran.r-project.org/package=genridge)
[![Rdoc](http://www.rdocumentation.org/badges/version/genridge)](http://www.rdocumentation.org/packages/genridge)

# genridge

## Generalized Ridge Trace Plots for Ridge Regression

Version 0.6-6


The genridge package introduces generalizations of the standard univariate
ridge trace plot used in ridge regression and related methods (Friendly, 2013)  These graphical methods
show both bias (actually, shrinkage) and precision, by plotting the covariance ellipsoids of the estimated
coefficients, rather than just the estimates themselves.  2D and 3D plotting methods are provided,
both in the space of the predictor variables and in the transformed space of the PCA/SVD of the
predictors.  

## References

Friendly, M. (2013).
The Generalized Ridge Trace Plot: Visualizing Bias *and* Precision.
*Journal of Computational and Graphical Statistics*, **22**(1), 50-68,
[doi](http://dx.doi.org/10.1080/10618600.2012.681237),
[genridge.pdf](http://euclid.psych.yorku.ca/datavis/papers/genridge.pdf)
