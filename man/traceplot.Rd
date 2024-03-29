% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/traceplot.R
\name{traceplot}
\alias{traceplot}
\title{Univariate Ridge Trace Plots}
\usage{
traceplot(
  x,
  X = c("lambda", "df"),
  col = c("black", "red", "darkgreen", "blue", "darkcyan", "magenta", "brown",
    "darkgray"),
  pch = c(15:18, 7, 9, 12, 13),
  xlab,
  ylab = "Coefficient",
  xlim,
  ylim,
  ...
)
}
\arguments{
\item{x}{A \code{ridge} object, as fit by \code{\link{ridge}}}

\item{X}{What to plot as the horizontal coordinate, one of \code{c("lambda", "df")}}

\item{col}{A numeric or character vector giving the colors used to plot the
ridge trace curves.  Recycled as necessary.}

\item{pch}{Vector of plotting characters used to plot the ridge trace
curves.  Recycled as necessary.}

\item{xlab}{Label for horizontal axis}

\item{ylab}{Label for vertical axis}

\item{xlim, ylim}{x, y limits for the plot}

\item{\dots}{Other arguments passed to \code{\link[graphics]{matplot}}}
}
\value{
None. Used for its side effect of plotting.
}
\description{
The \code{traceplot} function extends and simplifies the univariate ridge
trace plots for ridge regression provided in the \code{plot} method for
\code{\link[MASS]{lm.ridge}}
}
\details{
For ease of interpretation, the variables are labeled at the side of the
plot (left, right) where the coefficient estimates are expected to be most
widely spread.  If \code{xlim} is not specified, the range of the \code{X}
variable is extended slightly to accommodate the variable names.
}
\examples{

longley.y <- longley[, "Employed"]
longley.X <- data.matrix(longley[, c(2:6,1)])

lambda <- c(0, 0.005, 0.01, 0.02, 0.04, 0.08)
lridge <- ridge(longley.y, longley.X, lambda=lambda)

traceplot(lridge)
#abline(v=lridge$kLW, lty=3)
#abline(v=lridge$kHKB, lty=3)
#text(lridge$kLW, -3, "LW")
#text(lridge$kHKB, -3, "HKB")

traceplot(lridge, X="df")


}
\references{
Friendly, M. (2013). The Generalized Ridge Trace Plot:
Visualizing Bias \emph{and} Precision. \emph{Journal of Computational and
Graphical Statistics}, \bold{22}(1), 50-68,
doi:10.1080/10618600.2012.681237,
\url{https://www.datavis.ca/papers/genridge-jcgs.pdf}

Hoerl, A. E.  and Kennard R. W. (1970). "Ridge Regression: Applications to
Nonorthogonal Problems", \emph{Technometrics}, 12(1), 69-82.
}
\seealso{
\code{\link{ridge}} for details on ridge regression as implemented here

\code{\link{plot.ridge}}, \code{\link{pairs.ridge}} for other plotting
methods
}
\author{
Michael Friendly
}
\keyword{hplot}
