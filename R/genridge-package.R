
#' Generalized ridge trace plots for ridge regression
#' 
#' @description 
#' The genridge package introduces generalizations of the standard univariate
#' ridge trace plot used in ridge regression and related methods (Friendly,
#' 2012).  These graphical methods show both bias (actually, shrinkage) and
#' precision, by plotting the covariance ellipsoids of the estimated
#' coefficients, rather than just the estimates themselves.  2D and 3D plotting
#' methods are provided, both in the space of the predictor variables and in
#' the transformed space of the PCA/SVD of the predictors.
#'
#' @details  
#' 
#' This package provides computational support for the
#' graphical methods described in Friendly (2013). Ridge regression models may
#' be fit using the function \code{\link{ridge}}, which incorporates features
#' of \code{\link[MASS]{lm.ridge}}.  In particular, the shrinkage factors in
#' ridge regression may be specified either in terms of the constant added to
#' the diagonal of \eqn{X^T X} matrix (\code{lambda}), or the equivalent number
#' of degrees of freedom.
#' 
#' More importantly, the \code{\link{ridge}} function also calculates and
#' returns the associated covariance matrices of each of the ridge estimates,
#' allowing precision to be studied and displayed graphically.
#' 
#' This provides the support for the main plotting functions in the package:
#' 
#' \code{\link{plot.ridge}}: Bivariate ridge trace plots
#' 
#' \code{\link{pairs.ridge}}: All pairwise bivariate ridge trace plots
#' 
#' \code{\link{plot3d.ridge}}: 3D ridge trace plots
#' 
#' \code{\link{traceplot}}: Traditional univariate ridge trace plots
#' 
#' In addition, the function \code{\link{pca.ridge}} transforms the
#' coefficients and covariance matrices of a \code{ridge} object from predictor
#' space to the equivalent, but more interesting space of the PCA of \eqn{X^T
#' X} or the SVD of \bold{X}.  The main plotting functions also work for these
#' objects, of class \code{c("ridge", "pcaridge")}.
#' 
#' Finally, the functions \code{\link{precision}} and \code{\link{vif.ridge}}
#' provide other useful measures and plots.
#' 
#' @name genridge-package
#' @aliases genridge-package genridge
#' @docType package
#' @author Michael Friendly
#' 
#' Maintainer: Michael Friendly <friendly@@yorku.ca>
#' @seealso  \code{\link[MASS]{lm.ridge}}
#' 
#' @references 
#' Friendly, M. (2013). The Generalized Ridge Trace Plot:
#' Visualizing Bias \emph{and} Precision. \emph{Journal of Computational and
#' Graphical Statistics}, \bold{22}(1), 50-68,
#' doi:10.1080/10618600.2012.681237,
#' \url{https://www.datavis.ca/papers/genridge-jcgs.pdf}
#' 
#' Arthur E. Hoerl and Robert W. Kennard (1970). Ridge Regression: Biased
#' Estimation for Nonorthogonal Problems, \emph{Technometrics}, 12(1), pp.
#' 55-67.
#' 
#' Arthur E. Hoerl and Robert W. Kennard (1970). Ridge Regression: Applications
#' to Nonorthogonal Problems \emph{Technometrics}, 12(1), pp. 69-82.
#' 
#' @keywords package
#' @examples
#' 
#' # see examples for ridge, etc.
#' 
NULL
