# Add facility to fill contours in contour() and to return the contourLines()

#  All arguments of contour.default commented out can be passed to contour() as ...
#  fill.col is a vector of fill colors corresponding to levels,
#    or a call to the color.palette() function.
#    if not supplied, calculate a set of transparent colors
#    based on col and fill.alpha
#  DONE:  add logic for fill.col=NULL to nothing extra, other than return contourLines
#  DONE:  make fill.alpha accept a more standard 0-1 numeric



#' Enhanced Contour Plots
#' 
#' This is an enhancement to \code{\link[graphics]{contour}}, written as a
#' wrapper for that function.  It creates a contour plot, or adds contour lines
#' to an existing plot, allowing the contours to be filled and returning the
#' list of contour lines.
#' 
#' 
#' @param x,y locations of grid lines at which the values in \code{z} are
#'         measured.  These must be in ascending order.  By default, equally spaced
#'         values from 0 to 1 are used.  If \code{x} is a list, its components
#'         \code{x$x} and \code{x$y} are used for x and y, respectively.  If the list
#'         has component \code{x$z} this is used for \code{z}.
#' @param z a matrix containing the values to be plotted (NAs are allowed).
#'         Note that \code{x} can be used instead of \code{z} for convenience.
#' @param nlevels number of contour levels desired \bold{iff} levels is not
#'        supplied
#' @param levels numeric vector of levels at which to draw contour lines
#' @param zlim z-limits for the plot. x-limits and y-limits can be passed
#'       through \dots{}
#' @param col color for the lines drawn
#' @param color.palette a color palette function to be used to assign fill
#'        colors in the plot
#' @param fill.col a call to the \code{color.palette} function or an an
#'        explicit set of colors to be used in the plot. Use \code{fill.col=NULL} to
#'        suppress the filled polygons.  a vector of fill colors corresponding to
#'        levels.  By default, a set of possibly transparent colors is calculated
#'        ranging from white to \code{col}, using transparency given by
#'        \code{fill.alpha}
#' @param fill.alpha transparency value for \code{fill.col}, either a hex
#'        character string, or a numeric value between 0 and 1.  Use
#'        \code{fill.alpha=NA} to suppress transparency.
#' @param add logical. If \code{TRUE}, add to a current plot.
#' @param \dots additional arguments passed to \code{\link[graphics]{contour}},
#'        including all arguments of \code{\link[graphics]{contour.default}} not
#'        mentioned above, as well as additional graphical parameters passed by
#'        \code{\link[graphics]{contour.default}} to more basic functions.
#'        
#' @return Returns invisibly the list of contours lines, with components
#'        \code{levels}, \code{x}, \code{y}. See
#'        \code{\link[grDevices]{contourLines}}. 
#' 
#' @author Michael Friendly
#' @seealso \code{\link[graphics]{contour}},
#' \code{\link[grDevices]{contourLines}}
#' 
#' \code{\link[lattice]{contourplot}} from package lattice.
#' @keywords hplot
#' @examples
#' 
#' x <- 10*1:nrow(volcano)
#' y <- 10*1:ncol(volcano)
#' contourf(x,y,volcano, col="blue")
#' contourf(x,y,volcano, col="blue", nlevels=6)
#' 
#' # return value, unfilled, other graphic parameters
#' res <- contourf(x,y,volcano, col="blue", fill.col=NULL, lwd=2)
#' # levels used in the plot
#' sapply(res, function(x) x[[1]])
#' 
#' 
#' @export
contourf <- function(
                     x = seq(0, 1, length.out = nrow(z)),
                     y = seq(0, 1, length.out = ncol(z)),
                     z,
                     nlevels = 10,
                     levels = pretty(zlim, nlevels),
                     #    labels = NULL,
                     #    xlim = range(x, finite = TRUE),
                     #    ylim = range(y, finite = TRUE),
                     zlim = range(z, finite = TRUE),
                     #    labcex = 0.6, drawlabels = TRUE, method = "flattest",
                     #    vfont, axes = TRUE, frame.plot = axes,
                     col = par("fg"),
                     color.palette = colorRampPalette(c("white", col)),
                     fill.col = color.palette(nlevels + 1),
                     fill.alpha = 0.5, # alpha transparency
                     #    lty = par("lty"), lwd = par("lwd"),
                     add = FALSE, ...) {
  contour(x, y, z, nlevels = nlevels, levels = levels, zlim = zlim, col = col, add = add, ...)
  line.list <- contourLines(x, y, z, nlevels = nlevels, levels = levels)
  # contourLines returns a list of lists, each with components
  # 'level', 'x', 'y'

  if (!is.null(fill.col)) {
    if (!is.na(fill.alpha)) {
      if (is.numeric(fill.alpha) && fill.alpha >= 0 && fill.alpha <= 1)
      fill.alpha  <- as.hexmode(round(255 * fill.alpha))
      fill.col <- paste(fill.col, fill.alpha, sep = "")
    }

    levs <- sapply(line.list, function(x) x[[1]])
    for (i in seq_along(line.list)) {
      clev <- which(levs[i] == unique(levs))
      polygon(line.list[[i]][2:3], col = fill.col[clev], border = NA)
    }
  }
  invisible(line.list)
}
