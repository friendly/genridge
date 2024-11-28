#' Variance Inflation Factors for Ridge Regression
#' 
#' The function \code{vif.ridge} calculates variance inflation factors for the
#' predictors in a set of ridge regression models indexed by the
#' tuning/shrinkage factor, returning one row for each value of the $\\lambda$ parameter.
#' 
#' Variance inflation factors are calculated using the simplified formulation
#' in Fox & Monette (1992).
#' 
#' @param mod A \code{"ridge"} object
#' @param \dots Other arguments (unused) for compatibility with the generic from \code{\link[car]{vif}}
#' @return Returns a \code{"vif.ridge"} of variance inflation factors of the same size and
#'         shape as \code{coef{mod}}. The columns correspond to the predictors in the
#'         model and the rows correspond to the values of \code{lambda} in ridge
#'         estimation. 
#' @author Michael Friendly
#' @export
#' @importFrom car vif
#' @importFrom stats coef cov2cor vcov
#' @seealso \code{\link[car]{vif}}, \code{\link{precision}}
#' @references Fox, J. and Monette, G. (1992). Generalized collinearity
#' diagnostics. \emph{JASA}, \bold{87}, 178-183
#' @keywords models regression
#' @examples
#' 
#' data(longley)
#' lmod <- lm(Employed ~ GNP + Unemployed + Armed.Forces + Population + 
#'                       Year + GNP.deflator, data=longley)
#' vif(lmod)
#' 
#' lambda <- c(0, 0.005, 0.01, 0.02, 0.04, 0.08)
#' lridge <- ridge(Employed ~ GNP + Unemployed + Armed.Forces + 
#'                            Population + Year + GNP.deflator, 
#' 		             data=longley, lambda=lambda)
#'
#' coef(lridge)
#' 
#' # get VIFs for the shrunk estimates
#' vridge <- vif(lridge)
#' vridge
#' 
#' # plot VIFs
#' pch <- c(15:18, 7, 9)
#' clr <- c("black", rainbow(5, start=.6, end=.1))
#' 
#' matplot(rownames(vridge), vridge, type='b', 
#' 	xlab='Ridge constant (k)', ylab="Variance Inflation", 
#' 	xlim=c(0, 0.08), 
#' 	col=clr, pch=pch, cex=1.2)
#' text(0.0, vridge[1,], colnames(vridge), pos=4)
#' 
#' matplot(lridge$df, vridge, type='b', 
#' 	xlab='Degrees of freedom', ylab="Variance Inflation", 
#' 	col=clr, pch=pch, cex=1.2)
#' text(6, vridge[1,], colnames(vridge), pos=2)
#' 
#' # more useful to plot VIF on the sqrt scale
#' 
#' matplot(rownames(vridge), sqrt(vridge), type='b', 
#' 	xlab='Ridge constant (k)', ylab=expression(sqrt(VIF)), 
#' 	xlim=c(-0.01, 0.08), 
#' 	col=clr, pch=pch, cex=1.2, cex.lab=1.25)
#' text(-0.01, sqrt(vridge[1,]), colnames(vridge), pos=4, cex=1.2)
#' 
#' matplot(lridge$df, sqrt(vridge), type='b', 
#' 	xlab='Degrees of freedom', ylab=expression(sqrt(VIF)), 
#' 	col=clr, pch=pch, cex=1.2, cex.lab=1.25)
#' text(6, sqrt(vridge[1,]), colnames(vridge), pos=2, cex=1.2)
#' 
#' 
vif.ridge <- function(mod, ...) {

	Vif <- function(v) {
		R <- cov2cor(v)
		detR <- det(R)
		p <- nrow(v)
		res <- rep(0,p)
		for (i in 1:p) {
			res[i] <- R[i,i] * det(as.matrix(R[-i,-i])) / detR
		}
		res
	}

#	if(!require("car")) stop("Requires the car package for the vif generic")
	V <- vcov(mod)
	res <- t(sapply(V, Vif))
	colnames(res) <- colnames(coef(mod))
	rownames(res) <- rownames(coef(mod))
	res <- as.data.frame(res)
	res <- list(vif = res, lambda = mod$lambda, df = mod$df)
	class(res) <- "vif.ridge"
	res
}

#' @rdname vif.ridge
#' @exportS3Method print vif.ridge
print.vif.ridge <-
  function(x, digits = max(4, getOption("digits") - 5),...) {
    vif <- x$vif
    if (length(vif)) {
      cat("Variance inflaction factors:\n")
      print(format(vif, digits = digits), print.gap = 2, 
                    quote = FALSE)
    }
    invisible(x)
  }

