# tweaked labels
# axis3s(... nticks ...)


#' 3D Ridge Trace Plots
#' 
#' @name plot3d
#' @aliases plot3d plot3d.ridge plot3d.pcaridge
#' 
#' @description
#' The 3D ridge trace plot displays 3D projections of the covariance ellipsoids
#' for a set of ridge regression estimates indexed by a ridge tuning constant.
#' 
#' The centers of these ellipses show the bias induced for each parameter, and
#' also how the change in the ridge estimate for one parameter is related to
#' changes for other parameters.
#' 
#' The size and shapes of the covariance ellipsoids show directly the effect on
#' precision of the estimates as a function of the ridge tuning constant.
#' 
#' \code{plot3d.ridge} and \code{plot3d.pcaridge} differ only in the defaults
#' for the variables plotted.
#' 
#' @param x A \code{ridge} object, as fit by \code{\link{ridge}} or a
#'        \code{pcaridge} object as transformed by \code{\link{pca.ridge}}
#' @param variables Predictors in the model to be displayed in the plot: an
#'        integer or character vector of length 3, giving the indices or names of the
#'        variables. Defaults to the first three predictors for \code{ridge} objects
#'        or the \emph{last} three dimensions for \code{pcaridge} objects.
#' @param radius Radius of the ellipse-generating circle for the covariance
#'        ellipsoids.  The default, \code{radius=1} gives a standard \dQuote{unit}
#'        ellipsoid. Typically, \code{radius<1} gives less cluttered displays.
#' @param which.lambda A vector of indices used to select the values of
#'        \code{lambda} for which ellipsoids are plotted. The default is to plot
#'        ellipsoids for all values of \code{lambda} in the \code{ridge} object.
#' @param lwd,lty Line width and line type for the covariance ellipsoids.
#'        Recycled as necessary.
#' @param xlim,ylim,zlim X, Y, Z limits for the plot, each a vector of length
#'        2.  If missing, the range of the covariance ellipsoids is used.
#' @param xlab,ylab,zlab Labels for the X, Y, Z variables in the plot. If
#'        missing, the names of the predictors given in \code{variables} is used.
#' @param col A numeric or character vector giving the colors used to plot the
#'        covariance ellipsoids.  Recycled as necessary.
#' @param labels A numeric or character vector giving the labels to be drawn at
#'        the centers of the covariance ellipsoids.
#' @param ref Logical: whether to draw horizontal and vertical reference lines
#'        at 0. This is not yet implemented.
#' @param ref.col Color of reference lines.
#' @param segments Number of line segments used in drawing each dimension of a
#'        covariance ellipsoid.
#' @param shade a logical scalar or vector, indicating whether the ellipsoids
#'        should be rendered with \code{\link[rgl]{shade3d}}. Recycled as necessary.
#' @param shade.alpha a numeric value in the range [0,1], or a vector of such
#'        values, giving the alpha transparency for ellipsoids rendered with
#'        \code{shade=TRUE}.
#' @param wire a logical scalar or vector, indicating whether the ellipsoids
#'        should be rendered with \code{\link[rgl]{wire3d}}. Recycled as necessary.
#' @param aspect a scalar or vector of length 3, or the character string "iso",
#'        indicating the ratios of the x, y, and z axes of the bounding box.  The
#'        default, \code{aspect=1} makes the bounding box display as a cube
#'        approximately filling the display. See \code{\link[rgl]{aspect3d}} for
#'       details.
#' @param add if \code{TRUE}, add to the current \code{rgl} plot; the default
#'        is \code{FALSE}.
#' @param \dots Other arguments passed down
#' @return None. Used for its side-effect of plotting
#' 
#' @note This is an initial implementation.  The details and arguments are
#' subject to change.
#' @author Michael Friendly
#' @seealso \code{\link{plot.ridge}}, \code{\link{pairs.ridge}},
#' \code{\link{pca.ridge}}
#' @references Friendly, M. (2013). The Generalized Ridge Trace Plot:
#' Visualizing Bias \emph{and} Precision. \emph{Journal of Computational and
#' Graphical Statistics}, \bold{22}(1), 50-68,
#' doi:10.1080/10618600.2012.681237,
#' \url{https://www.datavis.ca/papers/genridge-jcgs.pdf}
#' @keywords hplot
#' @examples
#' 
#' lmod <- lm(Employed ~ GNP + Unemployed + Armed.Forces + Population + 
#'                       Year + GNP.deflator, data=longley)
#' longley.y <- longley[, "Employed"]
#' longley.X <- model.matrix(lmod)[,-1]
#' 
#' lambda <- c(0, 0.005, 0.01, 0.02, 0.04, 0.08)
#' lambdaf <- c("0", ".005", ".01", ".02", ".04", ".08")
#' lridge <- ridge(longley.y, longley.X, lambda=lambda)
#' 
#' \donttest{
#' plot3d(lridge, var=c(1,4,5), radius=0.5)
#' 
#' # view in SVD/PCA space
#' plridge <- pca(lridge)
#' plot3d(plridge, radius=0.5)
#' }
#' 
#' 
#'
#' @export
plot3d <-
  function (x, ...) {
    UseMethod("plot3d")
  }

# for pcaridge objects, default to last 3 variables
#' @rdname plot3d
#' @exportS3Method 
plot3d.pcaridge <-
  function(x, variables=(p-2):p, ...) {
    p <- dim(coef(x))[2]
    plot3d.ridge(x, variables, ...)
  }


#' @rdname plot3d
#' @import rgl
#' @exportS3Method 
plot3d.ridge <-
function(x, variables=1:3, radius=1, which.lambda=1:length(x$lambda),
		lwd=1, lty=1, 
		xlim, ylim, zlim,
		xlab, ylab, zlab,
		col = c("black", "red", "darkgreen", "blue","darkcyan","magenta", "brown","darkgray"),
#		c("black", rainbow(n.ell, start=.6, end=.1)),    # or, use rainbow colors by default?
		labels=lambda, 
#		center.pch = 16, center.cex=1.5, 
		ref=TRUE, ref.col=gray(.70), 
		segments=40,          # line segments in each ellipse
		shade=TRUE,           # use shade3d to render ellipsoids?
		shade.alpha=0.1,      # alpha transparency for shaded3d
		wire=FALSE,   				# use wire3d to render ellipsoids?
		aspect=1,							# ratios of the x, y, and z axes of the bounding box
#				bg.col=c("white", "black"),  # background colour
#				fogtype=c("none", "exp2", "linear", "exp"), # fog -- for depth cueing
#				fov=30,   # field of view (for perspective)
		add=FALSE,      # add to existing plot?
		...) {

	ellipsoid <- function(center, shape, radius=1, label="", col, lwd=1, shade=TRUE, alpha=0.1, wire=TRUE){
		# adapted from the shapes3d demo in the rgl package and from the heplots package
		degvec <- seq(0, 2*pi, length=segments)
		ecoord2 <- function(p) c(cos(p[1])*sin(p[2]), sin(p[1])*sin(p[2]), cos(p[2]))
		v <- t(apply(expand.grid(degvec,degvec), 1, ecoord2))

		warn <- options(warn=-1)
		on.exit(options(warn))
		Q <- chol(shape, pivot=TRUE)
		order <- order(attr(Q, "pivot"))
		v <- center + radius * t(v %*% Q[, order])
		v <- rbind(v, rep(1,ncol(v))) 
		e <- expand.grid(1:(segments-1), 1:segments)
		i1 <- apply(e, 1, function(z) z[1] + segments*(z[2] - 1))
		i2 <- i1 + 1
		i3 <- (i1 + segments - 1) %% segments^2 + 1
		i4 <- (i2 + segments - 1) %% segments^2 + 1
		i <- rbind(i1, i2, i4, i3)
		x <- rgl::asEuclidean(t(v))
		ellips <- rgl::qmesh3d(v, i)

		if(shade) rgl::shade3d(ellips, col=col, alpha=alpha, lit=TRUE)
		if(wire) rgl::wire3d(ellips, col=col, size=lwd, lit=FALSE)

		if (!is.null(label) && label !="")
			if(is.numeric(label)) label <- signif(label, 3)
			rgl::texts3d(center, adj=0.5, texts=label, color=col, lit=FALSE)
		bbox <- matrix(rgl::par3d("bbox"), nrow=2)
		rownames(bbox) <- c("min", "max")
		return(bbox)
	}

#	if (!require(rgl)) stop("rgl package is required.")
	if (!requireNamespace("rgl", quietly = TRUE)) stop("rgl package is required.")

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
	if(length(variables)>3) {
		warning("Only three variables will be plotted")
		variables <- variables[1:3]
	}
	vnames <- vnames[variables]
	lambda <- x$lambda[which.lambda]
	df <- x$df[which.lambda]                       # allow use of labels=df
	coef <- x$coef[which.lambda,variables]         # extract coefficients
	cov <- x$cov[which.lambda]                     # extract covariance matrices
	n.ell <- length(lambda)                        # number of ellipsoids to plot
	if (missing(xlab)) xlab <- vars[1]
	if (missing(ylab)) ylab <- vars[2]
	if (missing(zlab)) zlab <- vars[3]

#	fogtype <- match.arg(fogtype)
#	bg.col <- match.arg(bg.col)    
	if (!add){ 
		rgl::open3d()   
#		rgl.clear()
#		rgl.viewpoint(fov=fov)
#		rgl.bg(col=bg.col, fogtype=fogtype)    
	} 

	col <- rep(col, n.ell)		
	lwd <- rep(lwd, n.ell)		
	lty <- rep(lty, n.ell)		
	shade <- rep(shade, n.ell)
	wire <- rep(wire, n.ell)
	shade.alpha <- rep(shade.alpha, n.ell)

	ells <- as.list(rep(0, n.ell))
# generate the ellipses for each lambda, to get xlim & ylim
	for (i in 1:n.ell) {
		V <- (cov[[i]])[variables,variables]
		ells[[i]] <- ellipsoid(center=coef[i,], shape=V, radius=radius, label=labels[i],
			col=col[i], lwd=lwd, shade=shade[i], alpha=shade.alpha[i], wire=wire[i])
	}

	rgl::lines3d(coef, color="black")
	rgl::points3d(coef, size=5)
	rgl::aspect3d(aspect)		
	
	max <- apply(sapply(ells, function(X) apply(X, 2, max)), 1, max)
	min <- apply(sapply(ells, function(X) apply(X, 2, min)), 1, min)
	xlim <- if(missing(xlim)) c(min[1], max[1]) else xlim
	ylim <- if(missing(ylim)) c(min[2], max[2]) else ylim
	zlim <- if(missing(zlim)) c(min[3], max[3]) else zlim

	# handle xlim, ylim, zlim
	## enforce that the specified limits are at least as large as the bbox
	if (!missing(xlim) | !missing(ylim) | !missing(zlim)) {
		bbox <- matrix(rgl::par3d("bbox"),3,2,byrow=TRUE)
		xlim <- if(missing(xlim)) bbox[1,] else c(min(xlim[1],bbox[1,1]), max(xlim[2],bbox[1,2]))
		ylim <- if(missing(ylim)) bbox[2,] else c(min(ylim[1],bbox[2,1]), max(ylim[2],bbox[2,2]))
		zlim <- if(missing(zlim)) bbox[3,] else c(min(zlim[1],bbox[3,1]), max(zlim[2],bbox[3,2]))
		rgl::decorate3d(xlim=xlim, ylim=ylim, zlim=zlim, box=FALSE, axes=FALSE, xlab=NULL, ylab=NULL, zlab=NULL, top=FALSE)
	}
#	decorate3d(xlab=xlab, ylab=ylab, zlab=zlab, box=FALSE, axes=FALSE)
	frame <- rgl::axis3d("x-", col="black", tick=FALSE, nticks=0)
	frame <- c(frame, rgl::mtext3d(xlab, "x-", col="black", line=1.5))
	frame <- c(frame, rgl::axis3d("y-", col="black"), tick=FALSE, nticks=0)
	frame <- c(frame, rgl::mtext3d(ylab, "y-", col="black", line=1.5))
	frame <- c(frame, rgl::axis3d("z-", col="black"), tick=FALSE, nticks=0)
	frame <- c(frame, rgl::mtext3d(zlab, "z-", col="black", line=1.5))
	frame <- c(frame, rgl::box3d(col="black"))

}
