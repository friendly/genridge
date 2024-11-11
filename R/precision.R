# Nov 10, 2024: result gains class "precision" for a plot method

#' Measures of Precision and Shrinkage for Ridge Regression
#' 
#' @name precision
#' @aliases precision precision.ridge precision.lm
#'
#' @description
#' 
#' Three measures of (inverse) precision based on the \dQuote{size} of the
#' covariance matrix of the parameters are calculated. Let \eqn{V_k} be the
#' covariance matrix for a given ridge constant, and let \eqn{\lambda_i , i= 1,
#' \dots p} be its eigenvalues 
#' \enumerate{ 
#'   \item \eqn{\log | V_k | = \log \prod \lambda} or \eqn{|V_k|^{1/p} =(\prod \lambda)^{1/p}} 
#'        measures the linearized
#'        volume of the covariance ellipsoid and corresponds conceptually to Wilks'
#'        Lambda criterion 
#'   \item \eqn{ \text{trace}( V_k ) = \sum \lambda} corresponds conceptually to Pillai's trace criterion 
#'   \item \eqn{ \lambda_1 = \max (\lambda)} corresponds to Roy's largest root criterion.  
#' }
#' 
#' @param object An object of class \code{ridge} or \code{lm}
#' @param det.fun Function to be applied to the determinants of the covariance
#'        matrices, one of \code{c("log","root")}.
#' @param normalize If \code{TRUE} the length of the coefficient vector is
#'        normalized to a maximum of 1.0.
#' @param \dots Other arguments (currently unused)
#' 
#' @return An object of class \code{c("precision", "data.frame")} with the following columns: 
#' \item{lambda}{The ridge constant} 
#' \item{df}{The equivalent effective degrees of freedom} 
#' \item{det}{The \code{det.fun} function of the determinant of the covariance matrix} 
#' \item{trace}{The trace of the covariance matrix}
#' \item{max.eig}{Maximum eigen value of the covariance matrix}
#' \item{norm.beta}{The root mean square of the estimated coefficients, possibly normalized} 
#' 
#' @note Models fit by \code{lm} and \code{ridge} use a different scaling for
#' the predictors, so the results of \code{precision} for an \code{lm} model
#' will not correspond to those for \code{ridge} with ridge constant = 0.
#' 
#' @author Michael Friendly
#' @export
#' @seealso \code{\link{ridge}}, \code{\link{plot.precision}}
#' @keywords regression models
#' @examples
#' 
#' longley.y <- longley[, "Employed"]
#' longley.X <- data.matrix(longley[, c(2:6,1)])
#' 
#' lambda <- c(0, 0.005, 0.01, 0.02, 0.04, 0.08)
#' lridge <- ridge(longley.y, longley.X, lambda=lambda)
#'
#' # same, using formula interface
#' lridge <- ridge(Employed ~ GNP + Unemployed + Armed.Forces + Population + Year + GNP.deflator, 
#' 		data=longley, lambda=lambda)
#' 
#' clr <- c("black", rainbow(length(lambda)-1, start=.6, end=.1))
#' coef(lridge)
#' 
#' (pdat <- precision(lridge))
#' # plot log |Var(b)| vs. length(beta)
#' with(pdat, {
#' 	plot(norm.beta, det, type="b", 
#' 	cex.lab=1.25, pch=16, cex=1.5, col=clr, lwd=2,
#' 	xlab='shrinkage: ||b|| / max(||b||)',
#' 	ylab='variance: log |Var(b)|')
#' 	text(norm.beta, det, lambda, cex=1.25, pos=c(rep(2,length(lambda)-1),4))
#' 	text(min(norm.beta), max(det), "Variance vs. Shrinkage", cex=1.5, pos=4)
#' 	})
#' 
#' # plot trace[Var(b)] vs. length(beta)
#' with(pdat, {
#' 	plot(norm.beta, trace, type="b",
#' 	cex.lab=1.25, pch=16, cex=1.5, col=clr, lwd=2,
#' 	xlab='shrinkage: ||b|| / max(||b||)',
#' 	ylab='variance: trace [Var(b)]')
#' 	text(norm.beta, trace, lambda, cex=1.25, pos=c(2, rep(4,length(lambda)-1)))
#' #	text(min(norm.beta), max(det), "Variance vs. Shrinkage", cex=1.5, pos=4)
#' 	})
#' 
#' 

precision <- function(object, det.fun, normalize, ...) {
	UseMethod("precision")
}

# DONE:  allow choice of log.det or det()^{1/p}

#' @exportS3Method 
precision.ridge <- function(object, 
                            det.fun=c("log","root"), 
                            normalize=TRUE, ...) {
	tr <- function(x) sum(diag(x))
	maxeig <- function(x) max(eigen(x)$values)
	
	V <- object$cov
	p <- ncol(coef(object))
	det.fun <- match.arg(det.fun)
	ldet <- unlist(lapply(V, det))
	ldet <- if(det.fun == "log") log(ldet) else ldet^(1/p)
	trace <- unlist(lapply(V, tr))
	meig <- unlist(lapply(V, maxeig))	
	norm <- sqrt(rowMeans(coef(object)^2))
	if (normalize) norm <- norm / max(norm)
	res <- data.frame(lambda=object$lambda, 
	           df=object$df, 
	           det=ldet, 
	           trace=trace, 
	           max.eig=meig, 
	           norm.beta=norm)
	class(res) <- c("precision", "data.frame")
	res
}

#' @exportS3Method 
precision.lm <- function(object, det.fun=c("log","root"), normalize=TRUE, ...) {
	V <- vcov(object)
	beta <- coefficients(object)
	if (names(beta[1]) == "(Intercept)") {
		V <- V[-1, -1]
		beta <- beta[-1]
	}
#	else warning("No intercept: precision may not be sensible.")
	p <- length(beta)
	det.fun <- match.arg(det.fun)
	ldet <- det(V)
	ldet <- if(det.fun == "log") log(ldet) else ldet^(1/p)
	trace <- sum(diag(V))
	meig <- max(eigen(V)$values)
	norm <- sqrt(mean(beta^2))
	if (normalize) norm <- norm / max(norm)
	res <- list(df=length(beta), 
	            det=ldet, 
	            trace=trace, 
	            max.eig=meig, 
	            norm.beta=norm)
	class(res) <- c("precision", "data.frame")
	unlist(res)
}

