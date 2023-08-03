# Plot methods to add variable vectors showing the 
# original variables in PCA/SVD space.

# Thx: Uwe Ligges for the code for calculating scale...

# for ridge objects, default to first 2 variables
# and show PCA vectors in variable space....


#' Biplot of Ridge Regression Trace Plot in SVD Space
#' 
#' \code{biplot.pcaridge} supplements the standard display of the covariance
#' ellipsoids for a ridge regression problem in PCA/SVD space with labeled
#' arrows showing the contributions of the original variables to the dimensions
#' plotted.
#' 
#' The biplot view showing the dimensions corresponding to the two
#' \emph{smallest} singular values is particularly useful for understanding how
#' the predictors contribute to shrinkage in ridge regression.
#' 
#' This is only a biplot in the loose sense that results are shown in two
#' spaces simultaneously -- the transformed PCA/SVD space of the original
#' predictors, and vectors representing the predictors projected into this
#' space.
#' 
#' \code{biplot.ridge} is a similar extension of \code{\link{plot.ridge}},
#' adding vectors showing the relation of the PCA/SVD dimensions to the plotted
#' variables.
#' 
#' \code{class("ridge")} objects use the transpose of the right singular
#' vectors, \code{t(x$svd.V)} for the dimension weights plotted as vectors. 
#' 
#' 
#' @aliases biplot.pcaridge biplot.ridge
#' @param x A \code{pcaridge} object computed by \code{\link{pca.ridge}} or a \code{ridge} object.
#' @param variables The dimensions or variables to be shown in the the plot.
#'        By default, the \emph{last} two dimensions, corresponding to the smallest
#'        singular values, are plotted for \code{class("pcaridge")} objects or the
#'         \emph{first} two variables for \code{class("ridge")} objects.
#' @param labels A vector of character strings or expressions used as labels
#'        for the ellipses. Use \code{labels=NULL} to suppress these.
#' @param asp Aspect ratio for the plot. The default value, \code{asp=1} helps
#'       ensure that lengths and angles are preserved in these plots. Use
#'       \code{asp=NA} to override this.
#' @param origin The origin for the variable vectors in this plot, a vector of
#'        length 2. If not specified, the function calculates an origin to make the
#'        variable vectors approximately centered in the plot window.
#' @param scale The scale factor for variable vectors in this plot. If not
#'        specified, the function calculates a scale factor to make the variable
#'        vectors approximately fill the plot window.
#' @param var.lab Labels for variable vectors. The default is the names of the
#'        predictor variables.
#' @param var.lwd,var.col,var.cex Line width, color and character size used to
#'        draw and label the arrows representing the variables in this plot.
#' @param xlab,ylab Labels for the plot dimensions.  If not specified,
#'        \code{prefix} and \code{suffix} are used to construct informative dimension
#'        labels.
#' @param prefix Prefix for labels of the plot dimensions.
#' @param suffix Suffix for labels of the plot dimensions. If
#'        \code{suffix=TRUE} the percent of variance accounted for by each dimension
#'        is added to the axis label.
#' @param \dots Other arguments, passed to \code{\link{plot.pcaridge}}
#' @return None 
#' 
#' @author Michael Friendly, with contributions by Uwe Ligges
#' @seealso \code{\link{plot.ridge}}, \code{\link{pca.ridge}}
#' @references 
#' Friendly, M. (2013). The Generalized Ridge Trace Plot:
#' Visualizing Bias \emph{and} Precision. 
#' \emph{Journal of Computational and Graphical Statistics}, \bold{22}(1), 50-68,
#' \doi{10.1080/10618600.2012.681237},
#' \url{https://datavis.ca/papers/genridge-jcgs.pdf}
#'
#' @importFrom graphics arrows
#' @importFrom stats coef
#' @keywords hplot
#' @examples
#' 
#' longley.y <- longley[, "Employed"]
#' longley.X <- data.matrix(longley[, c(2:6,1)])
#' 
#' lambda <- c(0, 0.005, 0.01, 0.02, 0.04, 0.08)
#' lridge <- ridge(longley.y, longley.X, lambda=lambda)
#' 
#' plridge <- pca(lridge)
#' 
#' plot(plridge, radius=0.5)
#' 
#' # same, with variable vectors
#' biplot(plridge, radius=0.5)
#' # add some other options
#' biplot(plridge, radius=0.5, var.col="brown", var.lwd=2, var.cex=1.2, prefix="Dimension ")
#' 
#' # biplots for ridge objects, showing PCA vectors
#' plot(lridge, radius=0.5)
#' biplot(lridge, radius=0.5)
#' biplot(lridge, radius=0.5, asp=NA)
#' 
#' 
#' @exportS3Method  
biplot.pcaridge <- function(x, variables=(p-1):p, labels=NULL, asp=1, 
		origin, scale, 
		var.lab=rownames(V), 
		var.lwd=1, var.col="black", var.cex=1,
		xlab, ylab,      # override prefix/suffix?
		prefix = "Dim ", # prefix for labels of PCA dimensions
		suffix = TRUE,   # add label suffix with PCA % ?
		...) {
	
	# more convenient versions of arrows() and text()
  Arrows <- function(xy, lenxy, length, angle, col, lwd=1) {
    arrows(xy[1], xy[2], xy[1] + lenxy[, 1], xy[2] + lenxy[, 2], 
           length = length, angle = angle, lwd = lwd, col = col)
  }
  Text <- function(xy, lenxy, text, col="black", cex=1) {
    text(xy[1] + lenxy[, 1], xy[2] + lenxy[, 2], text, col = col, cex = cex)
  }	

  coef <- coef(x)
  p <- ncol(coef)
  if (is.null(x$svd.V)) stop("x must have an svd.V component")
  V <- x$svd.V[, variables]	

	# add
	pct <- 100*x$svd.D^2 /(sum(x$svd.D^2))
	if (is.logical(suffix) & suffix)
		suffix <- paste( " (", round(pct[variables],3), "%)", sep="" ) else suffix <- NULL
	dimlab <- paste(prefix, variables, suffix, sep="")
	if (missing(xlab)) xlab=dimlab[1]
	if (missing(ylab)) ylab=dimlab[2]
	
	plot(x, variables=variables, labels=labels, asp=asp, xlab=xlab, ylab=ylab, ...)
	
	bbox <- matrix(par("usr"), 2, 2, dimnames=list(c("min", "max"),c("x", "y")))
	if(missing(origin)) origin <- colMeans(bbox)
	
	# plot variable vectors
	if(missing(scale)) {
		scale <- c(sapply(bbox[,"x"] - origin[1], function(dist) dist/V[,1]),
				sapply(bbox[,"y"] - origin[2], function(dist) dist/V[,2]))
		scale <- 0.95* min(scale[scale > 0])
		cat("Vector scale factor set to ", scale, "\n")
	}
	
	Arrows(origin, scale*V, angle=8, length=.1, col=var.col, lwd=var.lwd)
	Text(origin, 1.01*scale*V, var.lab, col=var.col, cex=var.cex)
}

#' @exportS3Method 
biplot.ridge <-
  function(x, variables=1:2, xlab, ylab, ...) {
    x$svd.V <- t(x$svd.V)
    vnames <- colnames(coef(x))[variables]
    if(missing(xlab)) xlab <- vnames[1]
    if(missing(ylab)) ylab <- vnames[2]
    
    biplot.pcaridge(x, variables, xlab=xlab, ylab=ylab, ...)
  }

