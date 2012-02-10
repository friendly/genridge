\name{contourf}
\alias{contourf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Enhanced Contour Plots
}
\description{
This is an enhancement to \code{\link[graphics]{contour}}, written as a wrapper for that function. 
It creates a contour plot, or adds contour lines to an existing plot, allowing the contours to
be filled and returning the list of contour lines.
}
\usage{
contourf(x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, length.out = ncol(z)), z, 
	nlevels = 10, levels = pretty(zlim, nlevels), 
	zlim = range(z, finite = TRUE), 
  col = par("fg"),
  color.palette = colorRampPalette(c("white", col)),
  fill.col = color.palette(nlevels+1),
	fill.alpha = 0.5, 
	add = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{locations of grid lines at which the values in \code{z} are measured. 
  	These must be in ascending order. 
  	By default, equally spaced values from 0 to 1 are used. 
  	If \code{x} is a list, its components \code{x$x} and \code{x$y} are used for x and y, respectively. 
  	If the list has component \code{x$z} this is used for \code{z}.
}
  \item{z}{
a matrix containing the values to be plotted (NAs are allowed). 
Note that \code{x} can be used instead of \code{z} for convenience.
}
  \item{nlevels}{number of contour levels desired \bold{iff} levels is not supplied}
  \item{levels}{numeric vector of levels at which to draw contour lines}

  \item{zlim}{
z-limits for the plot. x-limits and y-limits can be passed through \dots
}
  \item{col}{
color for the lines drawn
}
  \item{color.palette}{
a color palette function to be used to assign fill colors in the plot
}

  \item{fill.col}{
a call to the \code{color.palette} function or an	
an explicit set of colors to be used in the plot.
Use \code{fill.col=NULL} to suppress the filled polygons. 
a vector of fill colors corresponding to levels.  
By default, a set of possibly transparent colors is calculated
ranging from white to \code{col}, using transparency given by \code{fill.alpha}
}
  \item{fill.alpha}{
transparency value for \code{fill.col}, either a hex character string, or
a numeric value between 0 and 1.  Use \code{fill.alpha=NA} to suppress
transparency.
}
  \item{add}{logical. If \code{TRUE}, add to a current plot.}

  \item{\dots}{
additional arguments passed to \code{\link[graphics]{contour}}, including
all arguments of \code{\link[graphics]{contour.default}}  
not mentioned above, as well as additional graphical parameters passed
by \code{\link[graphics]{contour.default}} to more basic functions.
}
}
%\details{
%
%}
\value{
Returns invisibly the list of contours lines, with components \code{levels}, \code{x}, \code{y}.
See \code{\link[grDevices]{contourLines}}.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
Michael Friendly
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[graphics]{contour}}, \code{\link[grDevices]{contourLines}}

\code{\link[lattice]{contourplot}} from package lattice.
}
\examples{
x <- 10*1:nrow(volcano)
y <- 10*1:ncol(volcano)
contourf(x,y,volcano, col="blue")
contourf(x,y,volcano, col="blue", nlevels=6)

# return value, unfilled, other graphic parameters
res <- contourf(x,y,volcano, col="blue", fill.col=NULL, lwd=2)
# levels used in the plot
sapply(res, function(x) x[[1]])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
