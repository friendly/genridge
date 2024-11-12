# TODO: Add label.prefix code

#' Plot Bias vs Variance for Ridge Precision
#' 
#' This function uses the results of \code{\link{precision}} to
#' plot a measure of shrinkage of the coefficients in ridge regression against a selected measure
#' of their estimated sampling variance, so as to provide a direct visualization of the tradeoff
#' between bias and precision.
#'
#' @param x A data frame of class \code{"precision"} resulting from \code{\link{precision}} called
#'        on a \code{"ridge"} object. Named \code{x} only to conform with the \code{\link{plot}} generic.
#' @param xvar    The character name of the column to be used for the horizontal axis. Typically, this is the normalized sum 
#'        of squares of the coefficients (\code{"norm.beta"}) used as a measure of shrinkage / bias.
#' @param yvar    The character name of the column to be used for the vertical axis. One of 
#'        \code{c("det", "trace", "max.eig")}. See \code{\link{precision}} for definitions of these measures.
#' @param labels  The character name of the column to be used for point labels. One of \code{c("lambda", "df")}.
#' @param label.cex Character size for point labels.
#' @param label.prefix Character or expression prefix for the point labels. Not yet implemented.
#' @param criteria The vector of optimal shrinkage criteria from the \code{\link{ridge}} call to be added
#'        as points in the plot. 
#' @param pch  Plotting character for points
#' @param cex  Character size for points
#' @param col  Point colors
#' @param main Plot title
#' @param xlab Label for horizontal axis
#' @param ylab Label for vertical axis
#' @param ...  Other arguments passed to \code{\link{plot}}.
#'
#' @return     Returns nothing. Used for the side effect of plotting.
#' @author Michael Friendly
#' @seealso
#'   \code{\link{ridge}} for details on ridge regression as implemented here.
#'   \code{\link{precision}} for definitions of the measures
#'
#' @importFrom colorspace qualitative_hcl
#' @importFrom splines bs
#' @exportS3Method 
#'
#' @examples
#' lambda <- c(0, 0.001, 0.005, 0.01, 0.02, 0.04, 0.08)
#' lridge <- ridge(Employed ~ GNP + Unemployed + Armed.Forces + 
#'                   Population + Year + GNP.deflator, 
#'                 data=longley, lambda=lambda)
#' 
#' criteria <- lridge$criteria |> print()
#' 
#' pridge <- precision(lridge) |> print()
#' 
#' plot(pridge)
#' # also show optimal criteria
#' plot(pridge, criteria = criteria)
#'
#' # use degrees of freedom as point labels 
#' plot(pridge, labels = "df")
#' # show the trace measure
#' plot(pridge, yvar="trace")


plot.precision <- function(
    x, 
    xvar = "norm.beta", 
    yvar = c("det", "trace", "max.eig"),
    labels = c("lambda", "df"),
    label.cex = 1.25,
    label.prefix,
    criteria = NULL,
    pch = 16,
    cex = 1.5,
    col,
    main = NULL,
    xlab, ylab,
    ...) {
  
  if (!inherits(x, "precision")) stop('Object must be of class "precision", not ', class(x))
  yvar <- match.arg(yvar)
  labels <- match.arg(labels)
  if (missing(xlab)) xlab <- paste("shrinkage:", xvar)
  if (missing(ylab)) ylab <- paste("variance:", yvar)
  
  xvar <- x[, xvar]
  yvar <- x[, yvar]
  labs <- x[, labels]
  lambda <- x[, "lambda"]
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
  text(xvar, yvar, labs, 
       pos=c(4, rep(2, nl)), 
       cex = label.cex, 
       xpd = TRUE)
  
  if (!is.null(criteria)) {
    mod <- lm(cbind(yvar, xvar) ~ splines::bs(lambda, df=5), 
              data=pridge)
    pts  <- data.frame(lambda=criteria) 
    fit <- predict(mod, pts) 
    points(fit[,2:1], pch=15, col=gray(.50), cex=1.6)
    text(fit[,2:1], rownames(fit), pos=c(3, 1, 3), 
         cex=1.25, col=gray(.50))
  }
}

if (FALSE) {
  lambda <- c(0, 0.001, 0.005, 0.01, 0.02, 0.04, 0.08)
  lambdaf <- c(expression(~widehat(beta)^OLS), lambda[-1])
  lridge <- ridge(Employed ~ GNP + Unemployed + Armed.Forces + 
                    Population + Year + GNP.deflator, 
                  data=longley, lambda=lambda)
  
  criteria <- lridge$criteria
  pridge <- precision(lridge)
  
  plot(pridge)
  plot(pridge, criteria = criteria)
  
  mod <- lm(cbind(det, norm.beta) ~ splines::bs(lambda, df=5), 
            data=pridge)
  pts  <- data.frame(lambda=c(lridge$kHKB, 
                           lridge$kLW))
  fit <- predict(mod, pts)
  points(fit[,2:1], pch=15, col=gray(.50), cex=1.6)
  text(fit[,2:1], c("HKB", "LW"), pos=3, cex=1.5, col=gray(.50))
  
  plot(pridge, labels = "df")
  
  plot(pridge, yvar="trace")
  
}