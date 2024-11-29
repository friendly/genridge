# TODO: allow for a longer list of colors

#' Univariate Ridge Trace Plots
#' 
#' The \code{traceplot} function extends and simplifies the univariate ridge
#' trace plots for ridge regression provided in the \code{plot} method for
#' \code{\link[MASS]{lm.ridge}}
#' 
#' For ease of interpretation, the variables are labeled at the side of the
#' plot (left, right) where the coefficient estimates are expected to be most
#' widely spread.  If \code{xlim} is not specified, the range of the \code{X}
#' variable is extended slightly to accommodate the variable names.
#' 
#' @param x A \code{ridge} object, as fit by \code{\link{ridge}}
#' @param X What to plot as the horizontal coordinate, one of \code{c("lambda", "df")}
#' @param col A numeric or character vector giving the colors used to plot the
#'        ridge trace curves.  Recycled as necessary.
#' @param pch Vector of plotting characters used to plot the ridge trace
#'        curves.  Recycled as necessary.
#' @param xlab Label for horizontal axis
#' @param ylab Label for vertical axis
#' @param xlim,ylim x, y limits for the plot. You may need to adjust these to allow for the variable labels.
#' @param \dots Other arguments passed to \code{\link[graphics]{matplot}}
#' @return None. Used for its side effect of plotting.
#' @author Michael Friendly
#' @importFrom graphics matplot
#' @export
#' @seealso 
#' \code{\link{ridge}} for details on ridge regression as implemented here
#' 
#' \code{\link{plot.ridge}}, \code{\link{pairs.ridge}} for other plotting
#' methods
#' 
#' @references Friendly, M. (2013). The Generalized Ridge Trace Plot:
#' Visualizing Bias \emph{and} Precision. \emph{Journal of Computational and
#' Graphical Statistics}, \bold{22}(1), 50-68,
#' doi:10.1080/10618600.2012.681237,
#' \url{https://www.datavis.ca/papers/genridge-jcgs.pdf}
#' 
#' Hoerl, A. E.  and Kennard R. W. (1970). "Ridge Regression: Applications to
#' Nonorthogonal Problems", \emph{Technometrics}, 12(1), 69-82.
#' @keywords hplot
#' @examples
#' 
#' longley.y <- longley[, "Employed"]
#' longley.X <- data.matrix(longley[, c(2:6,1)])
#' 
#' lambda <- c(0, 0.005, 0.01, 0.02, 0.04, 0.08)
#' lridge <- ridge(longley.y, longley.X, lambda=lambda)
#' 
#' traceplot(lridge)
#' #abline(v=lridge$kLW, lty=3)
#' #abline(v=lridge$kHKB, lty=3)
#' #text(lridge$kLW, -3, "LW")
#' #text(lridge$kHKB, -3, "HKB")
#' 
#' traceplot(lridge, X="df")
#' 
#' 
traceplot <-
function(x, 
	X=c("lambda","df"), 
#	labels=c("left", "right"),
	col = c("black", "red", "darkgreen", "blue","darkcyan","magenta", "brown","darkgray"), 
	pch = c(15:18, 7, 9, 12, 13),
	xlab, ylab="Coefficient", 
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
	coef <- coef(x)
	K <- nrow(coef)
	if (missing(xlim)) xlim <- range(X)
	if (missing(ylim)) ylim <- range(coef)

#	labels <- match.arg(labels)
	if (labels == "left") {
		xlim[1] <- xlim[1] - .1 * diff(xlim)
		labx <- X[1]
		laby <- coef[1,]
	}
	else {
		xlim[2] <- xlim[2] + .1 * diff(xlim)
		labx <- X[1]
		laby <- coef[1,]
	}

	matplot(X, coef, 	type="b", xlim=xlim, ylim=ylim, ylab=ylab, xlab=xlab, col=col, pch=pch, ...)
	abline(h=0, lty=3)
	if (type=="lambda") {
		abline(v=x$kHKB, col="gray", lty=2)
		text(x$kHKB, ylim[1], "HKB", pos=3)
		abline(v=x$kLW, col="gray", lty=2)
		text(x$kLW, ylim[1], "LW", pos=3)
	}
	vnames <- colnames(coef)
	text(labx, laby, colnames(coef), pos=c(2,4)[1+(labels=="right")])
}

