\name{plot3d.ridge}
\alias{plot3d}
\alias{plot3d.ridge}
\alias{plot3d.pcaridge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
3D Ridge Trace Plots
}
\description{
The 3D ridge trace plot displays 3D projections of the covariance ellipsoids
for a set of ridge regression estimates indexed by a ridge tuning constant.

The centers of these ellipses show the bias induced for each parameter, and
also how the change in the ridge estimate for one parameter is related to
changes for other parameters.

The size and shapes of the covariance ellipsoids show directly the effect on 
precision of the estimates as a function of the ridge tuning constant.
}
\usage{
plot3d(x, ...)

\method{plot3d}{ridge}(x, variables = 1:3, radius = 1,  which.lambda=1:length(x$lambda),
	lwd = 1, lty = 1, 
	xlim, ylim, zlim, 
	xlab, ylab, zlab, 
	col = c("black", "red", "darkgreen", "blue", 
	        "darkcyan", "magenta", "brown", "darkgray"), 
	labels = lambda, 
	ref = TRUE, ref.col = gray(0.7), 
	segments = 40, shade = TRUE, shade.alpha = 0.1, 
	wire = FALSE, aspect=1, add = FALSE, ...)

\method{plot3d}{pcaridge}(x, variables = (p-2):p, ...)

}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{ridge} object, as fit by \code{\link{ridge}} or a \code{pcaridge} object
as transformed by \code{\link{pca.ridge}}
}
  \item{variables}{
Predictors in the model to be displayed in the plot: an integer or character
vector of length 3, giving the indices or names of the variables.
Defaults to the first three predictors for \code{ridge} objects or the
\emph{last} three dimensions for \code{pcaridge} objects.
}
  \item{radius}{
Radius of the ellipse-generating circle for the covariance ellipsoids. 
The default, \code{radius=1} gives a standard \dQuote{unit} ellipsoid.
Typically, \code{radius<1} gives less cluttered displays.
}
  \item{which.lambda}{
A vector of indices used to select the values of \code{lambda} for which ellipsoids are plotted.
The default is to plot ellipsoids for all values of \code{lambda} in the \code{ridge} object.
}
  \item{lwd, lty}{
Line width and line type for the covariance ellipsoids. Recycled as necessary.
}
  \item{xlim, ylim, zlim}{
X, Y, Z limits for the plot, each a vector of length 2.  
If missing, the range of the covariance ellipsoids is used.
}
  \item{xlab, ylab, zlab}{
Labels for the X, Y, Z variables in the plot. If missing, the names of the
predictors given in \code{variables} is used.
}
  \item{col}{
A numeric or character vector giving the colors used to plot the covariance ellipsoids.
 Recycled as necessary.
}
  \item{labels}{
A numeric or character vector giving the labels to be drawn at the centers of the covariance ellipsoids.
}
  \item{ref}{
Logical: whether to draw horizontal and vertical reference lines at 0.
This is not yet implemented.
}
  \item{ref.col}{Color of reference lines.}

  \item{segments}{Number of line segments used in drawing each dimension of a covariance ellipsoid.}
  
  \item{shade}{a logical scalar or vector, indicating whether the ellipsoids should be
    rendered with \code{\link[rgl]{shade3d}}. Recycled as necessary.}

  \item{shade.alpha}{a numeric value in the range [0,1], or a vector of such values, giving
    the alpha transparency for ellipsoids rendered with \code{shade=TRUE}.}

  \item{wire}{a logical scalar or vector, indicating whether the ellipsoids should be
    rendered with \code{\link[rgl]{wire3d}}. Recycled as necessary.}

  \item{aspect}{a scalar or vector of length 3, or the character string "iso", indicating
  	the ratios of the x, y, and z axes of the bounding box.  The default, \code{aspect=1}
  	makes the bounding box display as a cube approximately filling the display. See
  	\code{\link[rgl]{aspect3d}} for details.}

  \item{add}{if \code{TRUE}, add to the current  \code{rgl} plot; the default is \code{FALSE}.
}
  \item{\dots}{
Other arguments passed down
}
}
\details{
\code{plot3d.ridge} and \code{plot3d.pcaridge} differ only in the defaults for the variables
plotted.
}
\value{
None
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Friendly, M. (2012).
The Generalized Ridge Trace Plot: Visualizing Bias and Precision.
In press, \emph{Journal of Computational and Graphical Statistics}, 21.
}
\author{
Michael Friendly
}
\note{
This is an initial implementation.  The details and arguments are subject to change.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.ridge}}, \code{\link{pairs.ridge}}, \code{\link{pca.ridge}}
}
\examples{
lmod <- lm(Employed ~ GNP + Unemployed + Armed.Forces + Population + 
                      Year + GNP.deflator, data=longley)
longley.y <- longley[, "Employed"]
longley.X <- model.matrix(lmod)[,-1]

lambda <- c(0, 0.005, 0.01, 0.02, 0.04, 0.08)
lambdaf <- c("0", ".005", ".01", ".02", ".04", ".08")
lridge <- ridge(longley.y, longley.X, lambda=lambda)

plot3d(lridge, var=c(1,4,5), radius=0.5)

# view in SVD/PCA space
plridge <- pca.ridge(lridge)
plot3d(plridge, radius=0.5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
