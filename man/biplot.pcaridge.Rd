\name{biplot.pcaridge}
\alias{biplot.pcaridge}
\alias{biplot.ridge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Biplot of Ridge Regression Trace Plot in SVD Space
}
\description{
\code{biplot.pcaridge} supplements the standard display of the 
covariance ellipsoids for a ridge regression problem in PCA/SVD
space with labeled arrows showing the contributions of the
original variables to the dimensions plotted.

The biplot view showing the dimensions corresponding to the two
\emph{smallest} singular values is particularly useful for understanding
how the predictors contribute to shrinkage in ridge regression.

This is only a biplot in the loose sense that results are shown
in two spaces simultaneously -- the transformed PCA/SVD space
of the original predictors, and vectors representing the predictors
projected into this space.

\code{biplot.ridge} is a similar extension of \code{\link{plot.ridge}},
adding vectors showing the relation of the PCA/SVD dimensions to the
plotted variables.
}
\usage{
\method{biplot}{pcaridge}(x, variables = (p - 1):p, labels=NULL, asp = 1, 
	origin, scale, 
	var.lab = rownames(V), var.lwd = 1, var.col = "black", var.cex = 1, 
	xlab, ylab, prefix = "Dim ", suffix = TRUE, ...)

\method{biplot}{ridge}(x, variables = 1:2, xlab, ylab, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{pcaridge} object computed by \code{\link{pca.ridge}}
or a \code{ridge} object.
}
  \item{variables}{
The dimensions or variables to be shown in the the plot.  By default, the
\emph{last} two dimensions, corresponding to the smallest
singular values, are plotted for \code{class("pcaridge")} objects
or the \emph{first} two variables for \code{class("ridge")} objects.
}
  \item{labels}{
A vector of character strings or expressions used as labels for the ellipses.
Use \code{labels=NULL} to suppress these.
}
  \item{asp}{
Aspect ratio for the plot. The default value, \code{asp=1}
helps ensure that lengths and angles are preserved in these plots.
Use \code{asp=NA} to override this.
}
  \item{origin}{
The origin for the variable vectors in this plot, a vector of length 2. If
  not specified, the function calculates an origin to make the variable vectors 
  approximately centered in  the plot window. 
}
  \item{scale}{
The scale factor for variable vectors in this plot. If
  not specified, the function calculates a scale factor to make the variable vectors 
  approximately fill the plot window.  
}
  \item{var.lab}{
Labels for variable vectors. The default is the names of the predictor variables.
}
  \item{var.lwd, var.col, var.cex}{
Line width, color and character size used to draw and label the arrows representing the variables
in this plot.
}
   \item{xlab, ylab}{
Labels for the plot dimensions.  If not specified, \code{prefix} and \code{suffix} are used to
construct informative dimension labels.
}
  \item{prefix}{
 Prefix for labels of the plot dimensions.
}
  \item{suffix}{
Suffix for labels of the plot dimensions. If \code{suffix=TRUE}
  the percent of variance accounted for by each dimension is added to the axis label.
}
  \item{\dots}{
Other arguments, passed to \code{\link{plot.pcaridge}}
}
}
\details{
\code{class("ridge")} objects use the transpose of the right singular vectors, \code{t(x$svd.V)} for the dimension weights
plotted as vectors. 
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
# add some other options
biplot(plridge, radius=0.5, var.col="brown", var.lwd=2, var.cex=1.2, prefix="Dimension ")

# biplots for ridge objects, showing PCA vectors
plot(lridge, radius=0.5)
biplot(lridge, radius=0.5)
biplot(lridge, radius=0.5, asp=NA)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
