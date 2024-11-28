
#' Bivariate Ridge Trace Plots
#' 
#' @description 
#' The bivariate ridge trace plot displays 2D projections of the covariance
#' ellipsoids for a set of ridge regression estimates indexed by a ridge tuning
#' constant.
#' 
#' The centers of these ellipses show the bias induced for each parameter, and
#' also how the change in the ridge estimate for one parameter is related to
#' changes for other parameters.
#' 
#' The size and shapes of the covariance ellipses show directly the effect on
#' precision of the estimates as a function of the ridge tuning constant.
#' 
#' 
#' @aliases plot.ridge plot.pcaridge
#' @param x A \code{ridge} object, as fit by \code{\link{ridge}}
#' @param variables Predictors in the model to be displayed in the plot: an
#'        integer or character vector of length 2, giving the indices or names of the
#'        variables. Defaults to the first two predictors for \code{ridge} objects or
#'        the \emph{last} two dimensions for \code{pcaridge} objects.
#' @param radius Radius of the ellipse-generating circle for the covariance
#'        ellipsoids. The default, \code{radius=1} gives a standard \dQuote{unit}
#'        ellipsoid. Typically, values \code{radius<1} gives less cluttered displays.
#' @param which.lambda A vector of indices used to select the values of
#'        \code{lambda} for which ellipses are plotted. The default is to plot
#'        ellipses for all values of \code{lambda} in the \code{ridge} object.
#' @param labels A vector of character strings or expressions used as labels
#'        for the ellipses. Use \code{labels=NULL} to suppress these.
#' @param pos,cex Scalars or vectors of positions (relative to the ellipse
#'        centers) and character size used to label the ellipses
#' @param lwd,lty Line width and line type for the covariance ellipsoids.
#'        Recycled as necessary.
#' @param xlim,ylim X, Y limits for the plot, each a vector of length 2.  If
#'        missing, the range of the covariance ellipsoids is used.
#' @param col A numeric or character vector giving the colors used to plot the
#'        covariance ellipsoids.  Recycled as necessary.
#' @param center.pch Plotting character used to show the bivariate ridge
#'        estimates. Recycled as necessary.
#' @param center.cex Size of the plotting character for the bivariate ridge
#'        estimates
#' @param fill Logical vector: Should the covariance ellipsoids be filled?
#'        Recycled as necessary.
#' @param fill.alpha Numeric vector: alpha transparency value(s) in the range
#'        (0, 1) for filled ellipsoids. Recycled as necessary.
#' @param ref Logical: whether to draw horizontal and vertical reference lines
#'        at 0.
#' @param ref.col Color of reference lines.
#' @param \dots Other arguments passed down to
#'        \code{\link[graphics]{plot.default}}, e.g., \code{xlab}, \code{ylab}, and
#'        other graphic parameters.
#' @return None. Used for its side effect of plotting.
#' @author Michael Friendly
#' @importFrom graphics abline polygon text
#' @importFrom grDevices gray
#' @seealso \code{\link{ridge}} for details on ridge regression as implemented
#' here; \code{\link{pairs.ridge}}, \code{\link{traceplot}}, for basic plots.
#' 
#' \code{\link{pca.ridge}} for transformation of ridge regression estimates to PCA space.
#' \code{\link{biplot.pcaridge}} and \code{\link{plot3d.ridge}} for other
#' plotting methods
#' 
#' @references Friendly, M. (2013). The Generalized Ridge Trace Plot:
#' Visualizing Bias \emph{and} Precision. \emph{Journal of Computational and
#' Graphical Statistics}, \bold{22}(1), 50-68,
#' doi:10.1080/10618600.2012.681237,
#' \url{https://www.datavis.ca/papers/genridge-jcgs.pdf}
#' 
#' @keywords hplot
#' @examples
#' 
#' longley.y <- longley[, "Employed"]
#' longley.X <- data.matrix(longley[, c(2:6,1)])
#' 
#' lambda <- c(0, 0.005, 0.01, 0.02, 0.04, 0.08)
#' lambdaf <- c("", ".005", ".01", ".02", ".04", ".08")
#' lridge <- ridge(longley.y, longley.X, lambda=lambda)
#' 
#' op <- par(mfrow=c(2,2), mar=c(4, 4, 1, 1)+ 0.1)
#' for (i in 2:5) {
#' 	plot(lridge, variables=c(1,i), radius=0.5, cex.lab=1.5)
#' 	text(lridge$coef[1,1], lridge$coef[1,i], expression(~widehat(beta)^OLS), 
#' 	     cex=1.5, pos=4, offset=.1)
#' 	if (i==2) text(lridge$coef[-1,1:2], lambdaf[-1], pos=3, cex=1.25)
#' }
#' par(op)
#' 
#' data(prostate)
#' py <- prostate[, "lpsa"]
#' pX <- data.matrix(prostate[, 1:8])
#' pridge <- ridge(py, pX, df=8:1)
#' 
#' plot(pridge)
#' plot(pridge, fill=c(TRUE, rep(FALSE,7)))
#' 
#' 
#' @exportS3Method 
plot.ridge <-
function(x, variables=1:2, radius=1, which.lambda=1:length(x$lambda), 
		labels=lambda, pos=3, cex=1.2,
		lwd=2, lty=1, xlim, ylim,
		col = c("black", "red", "darkgreen", "blue","darkcyan","magenta", "brown","darkgray"), 
		center.pch = 16, center.cex=1.5,
		fill=FALSE, fill.alpha=0.3, ref=TRUE, ref.col=gray(.70), ...) {

	ell <- function(center, shape, radius, segments=60) {
		angles <- (0:segments)*2*pi/segments
		circle <- radius * cbind( cos(angles), sin(angles))
		warn <- options(warn=-1)
		on.exit(options(warn))
		Q <- chol(shape, pivot=TRUE)
		order <- order(attr(Q, "pivot"))
		t( c(center) + t( circle %*% Q[,order]))
	}

	vnames <- dimnames(x$coef)[[2]]
	if (!is.numeric(variables)) {
		vars <- variables
		variables <- match(vars, vnames)
		check <- is.na(variables)
		if (any(check)) stop(paste(vars[check], collapse=", "), 
					" not among the predictor variables.") 
	}
	else {
		if (any (variables > length(vnames))) stop("There are only ", 
					length(vnames), " predictor variables.")
		vars <- vnames[variables]
	}
	if(length(variables)>2) {
		warning("Only two variables will be plotted. Perhaps you want plot3d.ridge()")
		variables <- variables[1:2]
	}
	lambda <- x$lambda[which.lambda]
	coef <- x$coef[which.lambda,variables]
	cov <- x$cov[which.lambda]
	n.ell <- length(lambda)
	lambda <- signif(lambda, 3)   # avoid many decimals when used as labels

	ells <- as.list(rep(0, n.ell))
# generate the ellipses for each lambda, to get xlim & ylim
	for (i in 1:n.ell) {
		ells[[i]] <- ell(center=coef[i,], shape=(cov[[i]])[variables,variables], radius=radius)
	}
	max <- apply(sapply(ells, function(X) apply(X, 2, max)), 1, max)
	min <- apply(sapply(ells, function(X) apply(X, 2, min)), 1, min)
	xlim <- if(missing(xlim)) c(min[1], max[1]) else xlim
	ylim <- if(missing(ylim)) c(min[2], max[2]) else ylim

	col <- rep(col, n.ell)		
	lwd <- rep(lwd, n.ell)		
	lty <- rep(lty, n.ell)		
	# handle filled ellipses
	fill <- rep(fill, n.ell)
	fill.alpha <- rep(fill.alpha, n.ell)
	fill.col <- trans.colors(col, fill.alpha)
	fill.col <- ifelse(fill, fill.col, NA)

	plot(coef, type='b', pch=center.pch, cex=center.cex, col=col, xlim=xlim, ylim=ylim, ...)
	if (ref) abline(v=0, h=0, col=ref.col)
	for (i in 1:n.ell) {
#		lines(ells[[i]], col=col[i], lwd=lwd[i], lty=lty[i])
		polygon(ells[[i]], col=fill.col[i], border=col[i],  lty=lty[i], lwd=lwd[i])
	}
	if(!is.null(labels)) text(coef, labels=labels, pos=pos, cex=cex)
}

# for pcaridge objects, default to last 2 variables
#'
#' @description
#' \code{plot.pcaridge} does these bivariate ridge trace plots for \code{"pcaridge"} objects, defaulting to plotting the
#' two smallest components.
#' @rdname plot.ridge
#' @exportS3Method 
plot.pcaridge <-
  function(x, variables=(p-1):p, labels=NULL, ...) {
    p <- dim(coef(x))[2]
    plot.ridge(x, variables, labels=labels, ...)
  }


