# TODO: allow for a longer list of colors

#' Univariate Ridge Trace Plots
#' 
#' The `traceplot` function extends and simplifies the univariate ridge
#' trace plots for ridge regression provided in the `plot` method for
#' [MASS::lm.ridge()]
#' 
#' For ease of interpretation, the variables are labeled at the side of the
#' plot (left, right) where the coefficient estimates are expected to be most
#' widely spread.  If `xlim` is not specified, the range of the `X`
#' variable is extended slightly to accommodate the variable names.
#' 
#' @param x A `ridge` object, as fit by [ridge()]
#' @param X What to plot as the horizontal coordinate, one of `c("lambda", "df")`
#' @param col A numeric or character vector giving the colors used to plot the
#'        ridge trace curves.  Recycled as necessary.
#' @param pch Vector of plotting characters used to plot the ridge trace
#'        curves.  Recycled as necessary.
#' @param xlab Label for horizontal axis
#' @param ylab Label for vertical axis
#' @param xlim,ylim x, y limits for the plot. You may need to adjust these to allow for the variable labels.
#' @param ... Other arguments passed to [graphics::matplot()]
#' @return None. Used for its side effect of plotting.
#' @author Michael Friendly
#' @importFrom graphics matplot strwidth par
#' @export
#' @seealso 
#' [ridge()] for details on ridge regression as implemented here
#' 
#' [plot.ridge()], [pairs.ridge()] for other plotting
#' methods
#' 
#' @references Friendly, M. (2013). The Generalized Ridge Trace Plot:
#' Visualizing Bias *and* Precision. *Journal of Computational and
#' Graphical Statistics*, **22**(1), 50-68,
#' \doi{10.1080/10618600.2012.681237}.
#' 
#' Hoerl, A. E.  and Kennard R. W. (1970). "Ridge Regression: Applications to
#' Nonorthogonal Problems", *Technometrics*, 12(1), 69-82.
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
traceplot <- function(x, ...) UseMethod("traceplot")

#' @rdname traceplot
#' @exportS3Method
traceplot.ridge <-
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

	# make room for the variable labels at the side where they are drawn
	xlim <- label_xlim(xlim, colnames(coef), side = labels)
	labx <- X[1]
	laby <- coef[1,]

	matplot(X, coef, 	type="b", xlim=xlim, ylim=ylim, ylab=ylab, xlab=xlab, col=col, pch=pch, ...)
	abline(h=0, lty=3)
	if (type=="lambda") {
		abline(v=x$kHKB, col="gray", lty=2)
		text(x$kHKB, ylim[1], "HKB", pos=3)
		abline(v=x$kLW, col="gray", lty=2)
		text(x$kLW, ylim[1], "LW", pos=3)
	}
	vnames <- colnames(coef)
	text(labx, laby, colnames(coef), pos=c(2,4)[1+(labels=="right")], xpd=TRUE)
}

# Extend xlim on one side so that labels drawn at the first point (with pos = 2 or 4) fit inside
# the plot region. The space needed is the width of the longest label (plus the text offset) as
# a fraction of the plot width, allowing for the default 4% axis extension (xaxs = "r") at each
# end. Never extends by less than the 10% used previously.
label_xlim <- function(xlim, labels, side = c("left", "right"), cex = par("cex")) {
  side <- match.arg(side)
  r <- diff(xlim)
  w <- max(strwidth(labels, units = "inches", cex = cex)) + strwidth("MM", units = "inches", cex = cex)
  f <- min(w / par("pin")[1], 0.45)                  # fraction of the plot width needed
  e <- r * (1.08 * f - 0.04) / (1.04 - 1.08 * f)     # extension giving that fraction
  e <- max(e, 0.1 * r)
  if (side == "left") xlim[1] <- xlim[1] - e else xlim[2] <- xlim[2] + e
  xlim
}

