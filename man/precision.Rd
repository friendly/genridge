\name{precision}
\alias{precision}
\alias{precision.ridge}
\alias{precision.lm}
\title{
Measures of Precision and Shrinkage for Ridge Regression
}
\description{
Calculates measures of precision based on the size of the 
estimated covariance matrices of the parameters and
shrinkage of the parameters in a ridge regression model.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{

precision(object, ...)

\method{precision}{ridge}(object, det.fun=c("log","root"), normalize=TRUE, ...)
\method{precision}{lm}(object, det.fun=c("log","root"), normalize=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{object}{
An object of class \code{ridge} or \code{lm}
}
  \item{det.fun}{
Function to be applied to the determinants of the covariance matrices, one of
\code{c("log","root")}.
}
  \item{normalize}{
If \code{TRUE} the length of the coefficient vector is normalized to a maximum of 1.0.
}
  \item{\dots}{
Other arguments (currently unused)
}

}
\details{
Three measures of (inverse) precision based on the \dQuote{size} of the
covariance matrix of the parameters are calculated. Let \eqn{V_k} be
the covariance matrix for a given ridge constant, and let
\eqn{\lambda_i , i= 1, \dots p} be its eigenvalues
\enumerate{
         \item \eqn{\log | V_k | = \log \prod \lambda} or \eqn{|V_k|^{1/p} =(\prod \lambda)^{1/p}}
         measures the linearized volume of the covariance ellipsoid
         and corresponds conceptually to Wilks' Lambda criterion
         \item \eqn{ trace( V_k ) =  \sum \lambda} corresponds conceptually to Pillai's trace criterion
         \item \eqn{ \lambda_1 =  max (\lambda)} corresponds to Roy's largest root criterion.
       }

}
\value{
A data.frame with the following columns
%%  If it is a LIST, use
  \item{lambda}{The ridge constant}
  \item{df}{The equivalent effective degrees of freedom}
  \item{det}{The \code{det.fun} function of the determinant of the covariance matrix}
  \item{trace}{The trace of the covariance matrix}
  \item{max.eig}{Maximum eigen value of the covariance matrix}
  \item{norm.beta}{The root mean square of the estimated coefficients, possibly normalized}
%% ...
}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
Michael Friendly
}
\note{
Models fit by \code{lm} and \code{ridge} use a different scaling for the
predictors, so the results of \code{precision} for an \code{lm} model
will not correspond to those for \code{ridge} with ridge constant = 0.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ridge}}, 
}
\examples{
longley.y <- longley[, "Employed"]
longley.X <- data.matrix(longley[, c(2:6,1)])

lambda <- c(0, 0.005, 0.01, 0.02, 0.04, 0.08)
lridge <- ridge(longley.y, longley.X, lambda=lambda)
coef(lridge)

(pdat <- precision(lridge))
with(pdat, {
	plot(norm.beta, det, type="b", 
	cex.lab=1.25, pch=16, cex=1.5, col=clr, lwd=2,
	xlab='shrinkage: ||b|| / max(||b||)',
	ylab='variance: log |Var(b)|')
	text(norm.beta, det, lambda, cex=1.25, pos=c(rep(2,length(lambda)-1),4))
	text(min(norm.beta), max(det), "Variance vs. Shrinkage", cex=1.5, pos=4)
	})


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
\keyword{models}% __ONLY ONE__ keyword per line
