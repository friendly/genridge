\name{genridge-package}
\alias{genridge-package}
\alias{genridge}
\docType{package}
\title{
Generalized ridge trace plots for ridge regression
}
\description{
The genridge package introduces generalizations of the standard univariate
	ridge trace plot used in ridge regression and related methods.  These graphical methods
	show both bias and precision, by plotting the covariance ellipsoids of the estimated
	coefficients, rather than just the estimates themselves.}
\details{
\tabular{ll}{
Package: \tab genridge\cr
Type: \tab Package\cr
Version: \tab 0.5-0\cr
Date: \tab 2011-08-22\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}
Ridge regression models may be fit using the function \code{\link{ridge}},
which incorporates features of  \code{\link[MASS]{lm.ridge}} and
 \code{\link[ElemStatLearn]{simple.ridge}}.  In particular, the shrinkage
factors in ridge regression may be specified either in terms of
the constant added to the diagonal of \eqn{X^T X} matrix (\code{lambda}),
or the equivalent number of degrees of freedom.  

More importantly, the  \code{\link{ridge}} function also calculates and
returns the associated covariance matrices of each of the ridge estimates.
This provides the support for the main plotting functions in the package:

\code{\link{plot.ridge}}: Bivariate ridge trace plots

\code{\link{pairs.ridge}}: All pairwise bivariate ridge trace plots

\code{\link{traceplot}}: Traditional univariate ridge trace plots

}

\author{
Michael Friendly

Maintainer: Michael Friendly <friendly@yorku.ca>
}
\references{
Friendly, M. (2011).
The Generalized Ridge Trace Plot: Visualizing Bias and Precision.
In preparation.

Arthur E. Hoerl and Robert W. Kennard (1970).
Ridge Regression: Biased Estimation for Nonorthogonal Problems,
\emph{Technometrics}, 12(1), pp. 55-67.

Arthur E. Hoerl and Robert W. Kennard (1970).
Ridge Regression: Applications to Nonorthogonal Problems
\emph{Technometrics}, 12(1), pp. 69-82.

}

\keyword{ package }
\seealso{
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
 \code{\link[MASS]{lm.ridge}},
 \code{\link[ElemStatLearn]{simple.ridge}}
 
}
\examples{
# see examples for ridge, etc.
}
