\name{vif.ridge}
\alias{vif.ridge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variance Inflation Factors for Ridge Regression
}
\description{
The function \code{vif.ridge} calculates variance inflation factors for the
predictors in a set of ridge regression models indexed by the tuning/shrinkage factor.
}
\usage{
\method{vif}{ridge}(mod, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{
A \code{ridge} object
}
  \item{\dots}{
Other arguments (unused)
}
}
\details{
Variance inflation factors are calculated using the simplified formulation
in Fox & Monette (1992).
}
\value{
Returns a matrix of variance inflation factors
of the same size and shape as \code{coef{mod}}.
The columns correspond to the predictors in the model and the
rows correspond to the values of \code{lambda} in ridge estimation.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Fox, J. and Monette, G. (1992). Generalized collinearity diagnostics. \emph{JASA}, \bold{87}, 178-183
}
\author{
Michael Friendly
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[car]{vif}}, \code{\link{precision}}
}
\examples{
data(longley)
lmod <- lm(Employed ~ GNP + Unemployed + Armed.Forces + Population + Year + GNP.deflator, data=longley)
vif(lmod)

longley.y <- longley[, "Employed"]
longley.X <- data.matrix(longley[, c(2:6,1)])

lambda <- c(0, 0.005, 0.01, 0.02, 0.04, 0.08)
lridge <- ridge(longley.y, longley.X, lambda=lambda)
coef(lridge)


vridge <- vif(lridge)
vridge

# plot VIFs
pch <- c(15:18, 7, 9)
clr <- c("black", rainbow(5, start=.6, end=.1))

matplot(rownames(vridge), vridge, type='b', 
	xlab='Ridge constant (k)', ylab="Variance Inflation", 
	xlim=c(0, 0.08), 
	col=clr, pch=pch, cex=1.2)
text(0.0, vridge[1,], colnames(vridge), pos=4)

matplot(lridge$df, vridge, type='b', 
	xlab='Degrees of freedom', ylab="Variance Inflation", 
	col=clr, pch=pch, cex=1.2)
text(6, vridge[1,], colnames(vridge), pos=2)

# more useful to plot VIF on the sqrt scale

matplot(rownames(vridge), sqrt(vridge), type='b', 
	xlab='Ridge constant (k)', ylab=expression(sqrt(VIF)), 
	xlim=c(-0.01, 0.08), 
	col=clr, pch=pch, cex=1.2, cex.lab=1.25)
text(-0.01, sqrt(vridge[1,]), colnames(vridge), pos=4, cex=1.2)

matplot(lridge$df, sqrt(vridge), type='b', 
	xlab='Degrees of freedom', ylab=expression(sqrt(VIF)), 
	col=clr, pch=pch, cex=1.2, cex.lab=1.25)
text(6, sqrt(vridge[1,]), colnames(vridge), pos=2, cex=1.2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}% __ONLY ONE__ keyword per line
