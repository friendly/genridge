\name{biplot.pcaridge}
\alias{biplot.pcaridge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Biplot of Ridge Regression Trace Plot in SVD Space
}
\description{
\code{biplot.pcaridge} supplements the standard display of the 
covariance ellipsoids for a ridge regression problem in PCA/SVD
space with labeled arrows showing the contributions of the
original variables to the dimensions plotted.

This is only a biplot in the loose sense that results are shown
in two spaces simultaneously -- the transformed PCA/SVD space
of the original predictors, and vectors representing the predictors
projected into this space.
}
\usage{
\method{biplot}{pcaridge}(x, variables = (p - 1):p, asp = 1, origin, scale, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{pcaridge} object
}
  \item{variables}{
The dimensions to be shown in the the plot.  By default, the
\emph{last} two dimensions, corresponding to the smallest
singular values, are plotted
}
  \item{asp}{
Aspect ratio for the plot. The default value, \code{asp=1}
helps ensure that lengths and angles are preserved in these plots.
}
  \item{origin}{
The origin for the variable vectors in this plot.
}
  \item{scale}{
The scale factor for variable vectors in this plot.
}
  \item{\dots}{
Other arguments, passed to \code{\link{plot.pcaridge}}
}
}
\details{
None yet.
%%  ~~ If necessary, more details than the description above ~~
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
Friendly, M. (2011).
The Generalized Ridge Trace Plot: Visualizing Bias and Precision.
Submitted, JCGS.
}
\author{
Michael Friendly, with contributions by Uwe Ligges
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.ridge}}, \code{\link{pca.ridge}}
}
\examples{
longley.y <- longley[, "Employed"]
longley.X <- data.matrix(longley[, c(2:6,1)])

lambda <- c(0, 0.005, 0.01, 0.02, 0.04, 0.08)
lridge <- ridge(longley.y, longley.X, lambda=lambda)

plridge <- pca.ridge(lridge)

plot(plridge, radius=0.5)
# same, with variable vectors
biplot(plridge, radius=0.5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
