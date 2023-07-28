#' Scatterplot Matrix of Bivariate Ridge Trace Plots
#' 
#' Displays all possible pairs of bivariate ridge trace plots for a given set
#' of predictors.
#' 
#' 
#' @param x A \code{ridge} object, as fit by \code{\link{ridge}}
#' @param variables Predictors in the model to be displayed in the plot: an
#'        integer or character vector, giving the indices or names of the variables.
#' @param radius Radius of the ellipse-generating circle for the covariance
#'        ellipsoids.
#' @param lwd,lty Line width and line type for the covariance ellipsoids.
#'        Recycled as necessary.
#' @param col A numeric or character vector giving the colors used to plot the
#'        covariance ellipsoids.  Recycled as necessary.
#' @param center.pch Plotting character used to show the bivariate ridge
#'        estimates. Recycled as necessary.
#' @param center.cex Size of the plotting character for the bivariate ridge
#'        estimates
#' @param fill Logical vector: Should the covariance ellipsoids be filled?
#'        Recycled as necessary.
#' @param fill.alpha Numeric vector: alpha transparency value(s) for filled
#'        ellipsoids. Recycled as necessary.
#' @param digits Number of digits to be displayed as the (min, max) values in
#'        the diagonal panels
#' @param diag.cex Character size for predictor labels in diagonal panels
#' @param diag.panel Function to draw diagonal panels.  Not yet implemented:
#'        just uses internal \code{panel.label} to write the variable name and ranges.
#' @param \dots Other arguments passed down
#' @return None. Used for its side effect of plotting.
#' @author Michael Friendly
#' @importFrom graphics par text box barplot
#' @export
#' @seealso \code{\link{ridge}} for details on ridge regression as implemented here
#' 
#' \code{\link{plot.ridge}}, \code{\link{traceplot}} for other plotting methods
#' 
#' @references Friendly, M. (2013). The Generalized Ridge Trace Plot:
#' Visualizing Bias \emph{and} Precision. \emph{Journal of Computational and
#' Graphical Statistics}, \bold{22}(1), 50-68,
#' doi:10.1080/10618600.2012.681237,
#' \url{http://euclid.psych.yorku.ca/datavis/papers/genridge.pdf}
#' 
#' @keywords hplot
#' @examples
#' 
#' longley.y <- longley[, "Employed"]
#' longley.X <- data.matrix(longley[, c(2:6,1)])
#' 
#' lambda <- c(0, 0.005, 0.01, 0.02, 0.04, 0.08)
#' lridge <- ridge(longley.y, longley.X, lambda=lambda)
#' 
#' pairs(lridge, radius=0.5, diag.cex=1.75)
#' 
#' data(prostate)
#' py <- prostate[, "lpsa"]
#' pX <- data.matrix(prostate[, 1:8])
#' pridge <- ridge(py, pX, df=8:1)
#' 
#' pairs(pridge)
#' 
pairs.ridge <-
function(x, variables, radius=1, lwd=1, lty=1,
		col = c("black", "red", "darkgreen", "blue","darkcyan","magenta", "brown","darkgray"), 
		center.pch = 16, center.cex=1.25, digits=getOption("digits") - 3,
		diag.cex= 2, diag.panel = panel.label,
		fill=FALSE, fill.alpha=0.3, ...) {

	panel.label <- function(x, ...) {
		op <- par(xpd=TRUE)
		on.exit(par(op))
		plot(c(min, max),c(min, max), type="n", axes=FALSE)
		text(0.5, 0.5, vars[i], cex=diag.cex)
		text(1, 0, signif(range[1, i], digits=digits), adj=c(1, 0))
		text(0, 1, signif(range[2, i], digits=digits), adj=c(0, 1)) 
		box()
	}	

#	why doesn't this work??
	panel.barplot <- function(x, ...) {
		barplot(x$coef[,i], axes=FALSE, col=col, ...)
		box()
	}

	vars <- dimnames(x$coef)[[2]]
  if (!missing(variables)){
      if (is.numeric(variables)) {
          vars <- vars[variables]
          if (any(is.na(vars))) stop("Bad variable selection.")
          }
      else {
          check <- !(variables %in% vars)
          if (any(check)) stop(paste("The following", 
              if (sum(check) > 1) "variables are" else "variable is",
              "not in the model:", paste(variables[check], collapse=", ")))
          vars <- variables
          }
      }
  else variables <- vars
	nvar <- length(vars)

  range <- apply(x$coef, 2, range)
	min=0
	max=1	
  old.par <- par(mfrow=c(nvar, nvar), mar=rep(0,4))
  on.exit(par(old.par))

#	browser()
	for (i in 1:nvar){
    for (j in 1:nvar){
      if (i == j)
				diag.panel(x)
      else {
        plot.ridge(x, variables=c(vars[j], vars[i]), radius=radius,
			labels=NULL,
        	col=col, lwd=lwd, lty=lty, center.cex=center.cex,
        	axes=FALSE,
          fill=fill, fill.alpha=fill.alpha, ...)
        box()
            }
        }
    }

}

