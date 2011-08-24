\name{ridge}
\alias{ridge}
\alias{coef.ridge}
\alias{print.ridge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ridge Regression Estimates
}
\description{
The function \code{ridge} fits linear models by ridge regression,
returning an object of class \code{ridge} designed to be used with
the plotting methods in this package.
}
\usage{

ridge(y, X, lambda = 0, df, svd = TRUE)

\method{print}{ridge}(x, digits = max(5, getOption("digits") - 5), ...)

\method{coef}{ridge}(object, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{y}{
A numeric vector containing the response variable. NAs not allowed.
}
  \item{X}{
A matrix of predictor variables. NA's not allowed. Should not include a column of 1's
for the intercept
}
  \item{lambda}{
A scalar or vector of ridge constants
}
  \item{df}{
A scalar or vector of effective degrees of freedom corresponding to \code{lambda}
}
  \item{svd}{
If \code{TRUE} the SVD of the centered and scaled \code{X} matrix is returned in
the \code{ridge} object.
}

  \item{x, object}{
An object of class \code{ridge}
}
  \item{\dots}{
Other arguments, passed down to methods
}

  \item{digits}{
For the \code{print} method, the number of digits to print.


}
}
\details{
Ridge regression shrinkage can be parameterized in several ways.
If a vector of \code{lambda} values is supplied, these are used directly in 
the ridge regression computations. Otherwise, if a vector \code{df} is supplied
the equivalent values of \code{lambda}.  In either case, both \code{lambda} and
\code{df} are returned in the \code{ridge} object.
}
\value{
A list with the following components:
%%  If it is a LIST, use
\item{lambda}{The vector of ridge constants}
\item{df}{The vector of effective degrees of freedom corresponding to \code{lambda}}
\item{coef}{The matrix of estimated ridge regression coefficients}
\item{scales}{scalings used on the X matrix}
\item{kHKB}{HKB estimate of the ridge constant}
\item{kLW}{L-W estimate of the ridge constant}

If \code{svd==TRUE}, the following are also included:

\item{svd.D}{Singular values of the \code{svd} of the scaled X matrix}
\item{svd.U}{Left singular vectors of the \code{svd} of the scaled X matrix}
\item{svd.V}{Right singular vectors of the \code{svd} of the scaled X matrix}
%% ...
}


\references{
Hoerl, A. E., Kennard, R. W., and Baldwin, K. F.
(1975), "Ridge Regression: Some Simulations," \emph{Communications
in Statistics}, 4, 105-123.

Lawless, J.F., and Wang, P. (1976), "A Simulation Study of Ridge
and Other Regression Estimators," \emph{Communications in Statistics},
5, 307-323.
}
\author{
Michael Friendly
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link[MASS]{lm.ridge}},
 \code{\link[ElemStatLearn]{simple.ridge}}
 for other implementations of ridge regression

\code{\link{plot.ridge}}, \code{\link{pairs.ridge}}, \code{\link{traceplot}} for plotting methods

}
\examples{
longley.y <- longley[, "Employed"]
longley.X <- data.matrix(longley[, c(2:6,1)])

lambda <- c(0, 0.005, 0.01, 0.02, 0.04, 0.08)
lridge <- ridge(longley.y, longley.X, lambda=lambda)
coef(lridge)

if (require("ElemStatLearn")) {
	py <- prostate[, "lpsa"]
	pX <- data.matrix(prostate[, 1:8])
	pridge <- ridge(py, pX, df=8:1)
	pridge
	
	plot(pridge)
	pairs(pridge)
	traceplot(pridge)
	traceplot(pridge, X="df")
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}% __ONLY ONE__ keyword per line
