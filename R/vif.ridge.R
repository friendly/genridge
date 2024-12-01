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
#' @return \code{vif} returns a \code{"vif.ridge"} object, which is a list of four components
#'   \item{vif}{a data frame of the same size and
#'         shape as \code{coef{mod}}. The columns correspond to the predictors in the
#'         model and the rows correspond to the values of \code{lambda} in ridge
#'         estimation.}
#'   \item{lambda}{the vector of ridge constants from the original call to \code{\link{ridge}} }
#'   \item{df}{the vector of effective degrees of freedom corresponding to \code{lambda}}
#'   \item{criteria}{the optimal values of \code{lambda}}
#' 
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
#' names(vridge)
#' 
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
#' @param Y What to plot as the vertical coordinate, one of \code{c("vif", "sqrt")}, where the latter plots \eqn{\sqrt{VIF}}.
#' @rdname vif.ridge
#' @exportS3Method plot vif.ridge
#' @examples
#' # plot VIFs
#' pch <- c(15:18, 7, 9)
#' clr <- c("black", rainbow(5, start=.6, end=.1))
#' 
#' plot(vridge, 
#'      col=clr, pch=pch, cex = 1.2,
#'      xlim = c(-0.02, 0.08))
#'
#' plot(vridge, X = "df",
#'      col=clr, pch=pch, cex = 1.2,
#'      xlim = c(4, 6.5))
#'
#' # Better to plot sqrt(VIF). Plot against degrees of freedom
#' plot(vridge, X = "df", Y="sqrt",
#'      col=clr, pch=pch, cex = 1.2,
#'      xlim = c(4, 6.5))
#'
#' 
#' 
plot.vif.ridge <-
  function(x, 
           X = c("lambda","df"),
           Y = c("vif", "sqrt"),
           col = c("black", "red", "darkgreen", "blue","darkcyan","magenta", "brown","darkgray"), 
           pch = c(15:18, 7, 9, 12, 13),
           xlab, 
           ylab, 
           xlim, ylim, ... ) {
    
    type <- X <- match.arg(X)
    if (type=="lambda") {
      X <- x$lambda
      if (missing(xlab)) xlab <- "Ridge constant (k)"
      labels <- "left"
    }
    else {
      X <- x$df
      if (missing(xlab)) xlab <- "Degrees of freedom"
      labels <- "right" 
    }

    Y <- match.arg(Y)
    vif <- if(Y == "vif") x$vif else sqrt(x$vif)
    if (missing(ylab)) ylab <- if(Y == "vif") "Variance Inflation" else expression(sqrt(VIF))

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
