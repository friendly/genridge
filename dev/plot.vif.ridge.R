plot.vif.ridge <- function(
    x,
    X=c("lambda","df"), 
    #	labels=c("left", "right"),
    col = c("black", "red", "darkgreen", "blue","darkcyan","magenta", "brown","darkgray"), 
    pch = c(15:18, 7, 9, 12, 13),
    xlab, ylab="Variance Inflation", 
    xlim, ylim, ... ) {

  type <- X <- match.arg(X)
  if (type=="lambda") {
    X <- x$lambda
    if (missing(xlab)) xlab <- "Ridge constant"
    labels <- "left"
  }
  else {
    X <- x$df
    if (missing(xlab)) xlab <- "Degrees of freedom"
    labels <- "right"
  }
  if (missing(xlim)) xlim <- range(X)
  if (missing(ylim)) ylim <- range(x)
  
  #	labels <- match.arg(labels)
  if (labels == "left") {
    xlim[1] <- xlim[1] - .1 * diff(xlim)
    labx <- X[1]
    laby <- x[1,]
  }
  else {
    xlim[2] <- xlim[2] + .1 * diff(xlim)
    labx <- X[1]
    laby <- x[1,]
  }
  
}