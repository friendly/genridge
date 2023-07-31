
#' Transform Ridge Estimates to PCA Space
#'
#' @name pca
#' @aliases pca pca.ridge
#' 
#' @description  
#' The function \code{pca.ridge} transforms a \code{ridge} object from
#' parameter space, where the estimated coefficients are \eqn{\beta_k} with
#' covariance matrices \eqn{\Sigma_k}, to the principal component space defined
#' by the right singular vectors, \eqn{V}, of the singular value decomposition
#' of the scaled predictor matrix, \eqn{X}.
#' 
#' In this space, the transformed coefficients are \eqn{V \beta_k}, with
#' covariance matrices \deqn{V \Sigma_k V^T}.
#' 
#' This transformation provides alternative views of ridge estimates in
#' low-rank approximations. In particular, it allows one to see where the
#' effects of collinearity typically reside --- in the smallest PCA dimensions.
#' 
#' 
#' @param x A \code{ridge} object, as fit by \code{\link{ridge}}
#' @param \dots Other arguments passed down. Not presently used in this
#'        implementation.
#' @return An object of class \code{c("ridge", "pcaridge")}, with the same
#'        components as the original \code{ridge} object. 
#' 
#' @author Michael Friendly
#' @seealso \code{\link{ridge}}
#' @references 
#' Friendly, M. (2013). The Generalized Ridge Trace Plot:
#' Visualizing Bias \emph{and} Precision. \emph{Journal of Computational and
#' Graphical Statistics}, \bold{22}(1), 50-68,
#' doi:10.1080/10618600.2012.681237,
#' \url{https://www.datavis.ca/papers/genridge-jcgs.pdf}
#' @export 
#' @keywords dplot multivariate
#' @examples
#' 
#' longley.y <- longley[, "Employed"]
#' longley.X <- data.matrix(longley[, c(2:6,1)])
#' 
#' lambda <- c(0, 0.005, 0.01, 0.02, 0.04, 0.08)
#' lridge <- ridge(longley.y, longley.X, lambda=lambda)
#' 
#' plridge <- pca(lridge)
#' traceplot(plridge)
#' pairs(plridge)
#' # view in space of smallest singular values
#' plot(plridge, variables=5:6)
#' 
#' 
pca <- function(x, ...){
  UseMethod("pca")
}

#' @exportS3Method 
pca.ridge <- function(x, ...) {
	if (is.null(x$svd.V)) stop("ridge object must contain svd.V")
	
	# transform a variance covariance matrix S by L
	rot <- function(S,L) L %*% S %*% t(L)
	
	result <- x	
	result$coef <- result$coef %*% result$svd.V
	result$cov <- lapply(result$cov, rot,  t(result$svd.V))
	class(result) <- c("pcaridge", class(x))
	result
}
