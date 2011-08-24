\name{traceplot}
\alias{traceplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Univariate ridge trace plots
}
\description{
The \code{traceplot} function extends and simplifies the univariate ridge trace plots
for ridge regresssion provided in the \code{plot} method for
\code{\link[MASS]{lm.ridge}} 

}
\usage{
traceplot(x, X = c("lambda", "df"), 
   col = c("black", "red", "darkgreen", "blue", "darkcyan", "magenta", "brown", "darkgray"), 
   pch = c(15:18, 7, 9, 12, 13), xlab, ylab = "Coefficient", xlim, ylim, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{ridge} object, as fit by \code{\link{ridge}}
}
  \item{X}{
What to plot as the horizontal coordinate, one of \code{c("lambda", "df")}
}
  \item{col}{
A numeric or character vector giving the colors used to plot the ridge trace curves.
 Recycled as necessary.
}
  \item{pch}{
Vector of plotting characters used to plot the ridge trace curves.
 Recycled as necessary.
}
  \item{xlab}{
Label for horizontal axis
}
  \item{ylab}{
Label for vertical axis
}
  \item{xlim, ylim}{
x, y limits for the plot
}

  \item{\dots}{
Other arguments passed to \code{\link[graphics]{matplot}}
}
}
\details{
For ease of interpretation, the variables are labeled at the side of the plot (left, right)
where the coefficient estimates are expected to be most widely spread.  If \code{xlim} is
not specified, the range of the \code{X} variable is extended slightly to accommodate the
variable names.
}

\value{
	None. Used for its side effect of plotting.
}
\references{
Hoerl, A. E.  and  Kennard R. W. (1970).
"Ridge Regression: Applications to Nonorthogonal Problems",
\emph{Technometrics},
12(1), 69-82.
}
\author{
Michael Friendly
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{ridge}} for details on ridge regression as implemented here

\code{\link{plot.ridge}}, \code{\link{pairs.ridge}} for other plotting methods
}
\examples{
longley.y <- longley[, "Employed"]
longley.X <- data.matrix(longley[, c(2:6,1)])

lambda <- c(0, 0.005, 0.01, 0.02, 0.04, 0.08)
lridge <- ridge(longley.y, longley.X, lambda=lambda)

traceplot(lridge)
traceplot(lridge, X="df")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
