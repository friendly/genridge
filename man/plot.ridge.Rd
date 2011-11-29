\name{plot.ridge}
\alias{plot.ridge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bivariate Ridge Trace Plots
}
\description{
The bivariate ridge trace plot displays 2D projections of the covariance ellipsoids
for a set of ridge regression estimates indexed by a ridge tuning constant.

The centers of these ellipses show the bias induced for each parameter, and
also how the change in the ridge estimate for one parameter is related to
changes for other parameters.

The size and shapes of the covariance ellipses show directly the effect on 
precision of the estimates as a function of the ridge tuning constant.

}
\usage{
\method{plot}{ridge}(x, variables = 1:2, radius = 1, lwd = 2, lty = 1, xlim, ylim, 
   col = c("black", "red", "darkgreen", "blue", "darkcyan", "magenta", "brown", "darkgray"), 
   center.pch = 16, center.cex = 1.5, fill = FALSE, fill.alpha = 0.3, ref=TRUE, ref.col=gray(.70), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{ridge} object, as fit by \code{\link{ridge}}
}
  \item{variables}{
Predictors in the model to be displayed in the plot: an integer or character
vector of length 2, giving the indices or names of the variables.
}
  \item{radius}{
Radius of the ellipse-generating circle for the covariance ellipsoids.  
}
  \item{lwd, lty}{
Line width and line type for the covariance ellipsoids. Recycled as necessary.
}
  \item{xlim, ylim}{
X, Y limits for the plot.  If missing, the range of the covariance ellipsoids is used.
}
  \item{col}{
A numeric or character vector giving the colors used to plot the covariance ellipsoids.
 Recycled as necessary.
}
  \item{center.pch}{
Plotting character used to show the bivariate ridge estimates. Recycled as necessary.
}
  \item{center.cex}{
Size of the plotting character for the bivariate ridge estimates
}
  \item{fill}{
Logical vector: Should the covariance ellipsoids be filled?  Recycled as necessary.
}
  \item{fill.alpha}{
Numeric vector: alpha transparency value(s) for filled ellipsoids. Recycled as necessary.
}
  \item{ref}{
Logical: whether to draw horizontal and vertical reference lines at 0.
}
  \item{ref.col}{
Color of reference lines.
}
  \item{\dots}{
Other arguments passed down
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{
	None. Used for its side effect of plotting.
}
\references{
Friendly, M. (2011).
The Generalized Ridge Trace Plot: Visualizing Bias and Precision.
Under review.
}
\author{
Michael Friendly
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ridge}} for details on ridge regression as implemented here

\code{\link{pairs.ridge}}, \code{\link{traceplot}} for other plotting methods
}
\examples{
longley.y <- longley[, "Employed"]
longley.X <- data.matrix(longley[, c(2:6,1)])

lambda <- c(0, 0.005, 0.01, 0.02, 0.04, 0.08)
lambdaf <- c("", ".005", ".01", ".02", ".04", ".08")
lridge <- ridge(longley.y, longley.X, lambda=lambda)

op <- par(mfrow=c(2,2), mar=c(4, 4, 1, 1)+ 0.1)

for (i in 2:5) {
	plot.ridge(lridge, variables=c(1,i), radius=0.5, cex.lab=1.5)
	text(lridge$coef[1,1], lridge$coef[1,i], expression(~widehat(beta)^OLS), cex=1.5, pos=4, offset=.1)
	if (i==2) text(lridge$coef[-1,1:2], lambdaf[-1], pos=3, cex=1.25)
}
par(op)

if (require("ElemStatLearn")) {
	py <- prostate[, "lpsa"]
	pX <- data.matrix(prostate[, 1:8])
	pridge <- ridge(py, pX, df=8:1)
	
	plot(pridge)
	plot(pridge, fill=c(TRUE, rep(FALSE,7)))
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
