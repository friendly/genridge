#' Plot Bias vs Variance for Ridge Precision
#'
#' @param x A data frame of class \code{"precision"} resulting from \code{link{precision}} called
#'        on a \code{"ridge"} object. Named \code{x} only to conform with the \code{\link{plot}} generic.
#' @param xvar    The character name of the column to be used for the horizontal axis.
#' @param yvar    The character name of the column to be used for the vertical axis.
#' @param labels The character name of the column to be used for point labels.
#' @param criteria The vector of optimal shrinkage criteria from the \code{link{ridge}} call to be added
#'        as labels in the plot. Not yet implemented.
#' @param pch  Plotting character for points
#' @param cex  Character size for points
#' @param col  Point colors
#' @param main Plot title
#' @param xlab Label for horizontal axis
#' @param ylab Label for vertical axis
#' @param ...  Other arguments passed to \code{link{plot}}.
#'
#' @return     Returns nothing. Used for the side effect of plotting.
#' @importFrom colorspace qualitative_hcl
#' @exportS3Method 
#'
#' @examples
#' # none yet
plot.precision <- function(
    x, 
    xvar = "norm.beta", 
    yvar = c("det", "trace", "max.eig"),
    labels = c("lambda", "df"),
    criteria = NULL,
    pch = 16,
    cex = 1.5,
    col,
    main = NULL,
    xlab, ylab,
    ...) {
  yvar <- match.arg(yvar)
  labels <- match.arg(labels)
  if (missing(xlab)) xlab <- paste("shrinkage:", xvar)
  if (missing(ylab)) ylab <- paste("variance:", yvar)
  
  xvar <- x[, xvar]
  yvar <- x[, yvar]
  labs <- x[, labels]
  if (labels=="df") labs <- round(labs, 2)
  labs[1] <- "OLS"     # expression(~widehat(beta)^OLS)  # or, maybe just "OLS"
  nl <- length(labs)
  
  if (missing(col)) col <- c("black", 
                             colorspace::qualitative_hcl(nl-1L, palette = "Dark 3"))
  
  plot(xvar, yvar, type = "b", 
       pch = pch,
       col = col,
       cex = cex,
       lwd = 2,
       xlab = xlab, ylab = ylab,
       main = main,
       ...)
  text(xvar, yvar, labs, pos=c(4, rep(2, nl)), cex = 1.25, xpd = TRUE)
}

if (FALSE) {
  lambda <- c(0, 0.001, 0.005, 0.01, 0.02, 0.04, 0.08)
  lambdaf <- c(expression(~widehat(beta)^OLS), lambda[-1])
  lridge <- ridge(Employed ~ GNP + Unemployed + Armed.Forces + 
                    Population + Year + GNP.deflator, 
                  data=longley, lambda=lambda)
  pridge <- precision(lridge)
  
  plot(pridge)
  plot(pridge, labels = "df")
  
  plot(pridge, yvar="trace")
  
}