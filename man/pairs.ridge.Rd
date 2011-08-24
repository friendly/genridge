\name{pairs.ridge}
\alias{pairs.ridge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scatterplot Matrix of Bivariate Ridge Trace Plots
}
\description{
Displays all possible pairs of bivariate ridge trace plots for a given set of predictors.
}
\usage{
\method{pairs}{ridge}(x, variables, radius = 1, lwd = 1, lty = 1, 
   col = c("black", "red", "darkgreen", "blue", "darkcyan", "magenta", "brown", "darkgray"), 
   center.pch = 16, center.cex = 1.25, digits = getOption("digits") - 3, 
   diag.cex = 2, diag.panel = panel.label, fill = FALSE, fill.alpha = 0.3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{ridge} object, as fit by \code{\link{ridge}}
}
  \item{variables}{
Predictors in the model to be displayed in the plot: an integer or character
vector, giving the indices or names of the variables.
}
  \item{radius}{
Radius of the ellipse-generating circle for the covariance ellipsoids.  
}
  \item{lwd, lty}{
Line width and line type for the covariance ellipsoids. Recycled as necessary.
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
  \item{digits}{
Number of digits to be displayed as the (min, max) values in the diagonal panels
}
  \item{diag.cex}{
Character size for predictor labels in diagonal panels
}
  \item{diag.panel}{
Function to draw diagonal panels.  Not yet implemented: just uses internal \code{panel.label}
to write the variable name and ranges.
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
In preparation.
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

\code{\link{plot.ridge}}, \code{\link{traceplot}} for other plotting methods
}
\examples{
longley.y <- longley[, "Employed"]
longley.X <- data.matrix(longley[, c(2:6,1)])

lambda <- c(0, 0.005, 0.01, 0.02, 0.04, 0.08)
lridge <- ridge(longley.y, longley.X, lambda=lambda)

pairs(lridge, radius=0.5, diag.cex=1.75)

if (require("ElemStatLearn")) {
	py <- prostate[, "lpsa"]
	pX <- data.matrix(prostate[, 1:8])
	pridge <- ridge(py, pX, df=8:1)
	
	pairs(pridge)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
