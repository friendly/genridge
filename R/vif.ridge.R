#' @title 
#' Variance Inflation Factors for Ridge Regression
#'
#' @description
#' The function \code{vif.ridge} calculates variance inflation factors for the
#' predictors in a set of ridge regression models indexed by the
#' tuning/shrinkage factor, returning one row for each value of the \eqn{\lambda} parameter.
#' 
#' Variance inflation factors are calculated using the simplified formulation
#' in Fox & Monette (1992).
#' 
#' @param mod A \code{"ridge"} object computed by \code{\link{ridge}}
#' @param \dots Other arguments passed to methods
#' @return Returns a \code{"vif.ridge"} of variance inflation factors of the same size and
#'         shape as \code{coef{mod}}. The columns correspond to the predictors in the
#'         model and the rows correspond to the values of \code{lambda} in ridge
#'         estimation. [Now returns a list!!! FIXME]
#' @author Michael Friendly
#' @export
#' @importFrom car vif
#' @importFrom stats coef cov2cor vcov
#' @seealso \code{\link[car]{vif}}, \code{\link{precision}}
#' @references Fox, J. and Monette, G. (1992). Generalized collinearity
#' diagnostics. \emph{JASA}, \bold{87}, 178-183, \doi{10.1080/01621459.1992.10475190}.
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
#' ### Transition examples, because the vif() method now returns a list structure 
#' ### rather than a data.frame 
#' vr <- vridge$vif
#' matplot(rownames(vr), vr, type='b', 
#' 	xlab='Ridge constant (k)', ylab="Variance Inflation", 
#' 	xlim=c(0, 0.08), 
#' 	col=clr, pch=pch, cex=1.2)
#' text(0.0, vr[1,], colnames(vr), pos=4)
#' 
#' # matplot(lridge$df, vridge, type='b', 
#' # 	xlab='Degrees of freedom', ylab="Variance Inflation", 
#' #	col=clr, pch=pch, cex=1.2)
#' # text(6, vridge[1,], colnames(vridge), pos=2)
#' 
#' # more useful to plot VIF on the sqrt scale
#' 
#' # matplot(rownames(vridge), sqrt(vridge), type='b', 
#' #	xlab='Ridge constant (k)', ylab=expression(sqrt(VIF)), 
#' #	xlim=c(-0.01, 0.08), 
#' #	col=clr, pch=pch, cex=1.2, cex.lab=1.25)
#' # text(-0.01, sqrt(vridge[1,]), colnames(vridge), pos=4, cex=1.2)
#' 
#' # matplot(lridge$df, sqrt(vridge), type='b', 
#' #	xlab='Degrees of freedom', ylab=expression(sqrt(VIF)), 
#' #	col=clr, pch=pch, cex=1.2, cex.lab=1.25)
#' # text(6, sqrt(vridge[1,]), colnames(vridge), pos=2, cex=1.2)
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
	res <- list(vif = res, lambda = mod$lambda, df = mod$df, criteria = mod$criteria)
	class(res) <- "vif.ridge"
	res
}

#' @param digits Number of digits to display in the \code{print} method
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


#' @description
#' 
#' The \code{plot.vif.ridge} method plots variance inflation factors for a \code{"vif.ridge"} object
#' in a similar style to what is provided by \code{\link{traceplot}}. That is, it plots the VIF for each
#' coefficient in the model against either the ridge \eqn{\lambda} tuning constant or it's equivalent
#' effective degrees of freedom.
#' 
#' @inheritParams traceplot
#' @rdname vif.ridge
#' @exportS3Method plot vif.ridge
plot.vif.ridge <-
  function(x, 
           X=c("lambda","df"), 
           col = c("black", "red", "darkgreen", "blue","darkcyan","magenta", "brown","darkgray"), 
           pch = c(15:18, 7, 9, 12, 13),
           xlab, ylab="Variance Inflation", 
           xlim, ylim, ... ) {
    
    type <- X <- match.arg(X)
    if (type=="lambda") {
      X <- x$lambda
      if (missing(xlab)) xlab <- "Ridge constant"
      labels <- "left"
    }
    else {
      X <- x$df
      if (missing(xlab)) xlab <- "Degrees of freedom"
      labels <- "right"
    }
    vif <- x$vif
    K <- nrow(vif)
    if (missing(xlim)) xlim <- range(X)
    if (missing(ylim)) ylim <- range(vif)
    
    #	labels <- match.arg(labels)
    if (labels == "left") {
      xlim[1] <- xlim[1] - .1 * diff(xlim)
      labx <- X[1]
      laby <- vif[1,]
    }
    else {
      xlim[2] <- xlim[2] + .1 * diff(xlim)
      labx <- X[1]
      laby <- vif[1,]
    }
    
    matplot(X, vif, 	type="b", xlim=xlim, ylim=ylim, ylab=ylab, xlab=xlab, col=col, pch=pch, ...)
    abline(h=0, lty=3)
    
    if (type=="lambda") {
      criteria <- x$criteria
      abline(v=criteria, col="gray", lty=2)
      text(criteria, ylim[1], names(criteria), pos=3)
    }
    vnames <- colnames(vif)
    text(labx, laby, colnames(vif), pos=c(2,4)[1+(labels=="right")], xpd=TRUE)
  }
