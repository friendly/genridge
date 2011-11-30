\name{pca.ridge}
\alias{pca.ridge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform Ridge Estimates to PCA Space
}
\description{
The function \code{pca.ridge} transforms a \code{ridge} object
from parameter space, where the estimated coefficients are
\eqn{\beta_k} with covariance matrices \eqn{\Sigma_k}, to the
principal component space defined by the right singular vectors, \eqn{V},
of the singular value decomposition of the scaled predictor matrix, X.

In this space, the transformed coefficients are \eqn{V \beta_k},
with covariance matrices
\deqn{V \Sigma_k V^T}

This transformation provides alternative views of ridge estimates
in low-rank approximations.

}
\usage{
pca.ridge(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{ridge} object, as fit by \code{\link{ridge}}
}
  \item{\dots}{
Other arguments passed down
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{
An object of class \code{c("ridge", "pcaridge")}, with the same
components as the original \code{ridge} object.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
\code{\link{ridge}}
}
\examples{
longley.y <- longley[, "Employed"]
longley.X <- data.matrix(longley[, c(2:6,1)])

lambda <- c(0, 0.005, 0.01, 0.02, 0.04, 0.08)
lridge <- ridge(longley.y, longley.X, lambda=lambda)

plridge <- pca.ridge(lridge)
pairs(plridge)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{dplot}
\keyword{multivariate}% __ONLY ONE__ keyword per line
