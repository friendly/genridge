# Simple version of a plot method to add variable vectors showing the 
# original variables in PCA/SVD space.

# Thx: Uwe Ligges for the code for calculating scale...

biplot.pcaridge <- function(x, variables=(p-1):p, asp=1, origin, scale, ...) {

	# more convenient versions of arrows() and text()
	Arrows <- function(xy, lenxy, length, angle, col, lwd=1) {
		arrows(xy[1], xy[2], xy[1]+lenxy[,1], xy[2]+lenxy[,2], length=length, angle=angle)
	}
	Text <- function(xy, lenxy, text, col="black") {
		text(xy[1]+lenxy[,1], xy[2]+lenxy[,2], text, col=col)
	}
	
	coef <- coef(x)
	p <- ncol(coef)
	if(is.null(x$svd.V)) stop("x must have an svd.V component")
	V <- x$svd.V[,variables]
	
	plot(x, variables=variables, asp=asp, ...)
	
	bbox <- matrix(par("usr"), 2, 2, dimnames=list(c("min", "max"),c("x", "y")))
	if(missing(origin)) origin <- colMeans(bbox)

	# plot variable vectors
	if(missing(scale)) {
		scale <- c(sapply(bbox[,"x"] - origin[1], function(dist) dist/V[,1]),
		           sapply(bbox[,"y"] - origin[2], function(dist) dist/V[,2]))
		scale <- 0.95* min(scale[scale > 0])
		cat("Vector scale factor set to ", scale, "\n")
	}

	Arrows(origin, scale*V, angle=8, length=.1)
	Text(origin, scale*V, rownames(V))
	
}

