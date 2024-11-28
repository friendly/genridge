
#' @description
#' 
#' The \code{plot.vif.ridge} method plots variance inflation factors for a \code{"vif.ridge"} object
#' in a similar style to what is provided by \code{\link{traceplot}}. That is, it plots the VIF for each
#' coefficient in the model against either the ridge \eqn{\lambda} tuning constant or it's equivalent
#' effective degrees of freedom.
#' 
#' @inheritParams traceplot
#' @rdname vif.ridge
#' @exportS3Method plot vif.ridge
plot.vif.ridge <-
  function(x, 
           X=c("lambda","df"), 
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
    vif <- x$vif
    K <- nrow(vif)
    if (missing(xlim)) xlim <- range(X)
    if (missing(ylim)) ylim <- range(vif)
    
    #	labels <- match.arg(labels)
    if (labels == "left") {
      xlim[1] <- xlim[1] - .1 * diff(xlim)
      labx <- X[1]
      laby <- vif[1,]
    }
    else {
      xlim[2] <- xlim[2] + .1 * diff(xlim)
      labx <- X[1]
      laby <- vif[1,]
    }
    
    matplot(X, vif, 	type="b", xlim=xlim, ylim=ylim, ylab=ylab, xlab=xlab, col=col, pch=pch, ...)
    abline(h=0, lty=3)

    if (type=="lambda") {
      criteria <- x$criteria
      abline(v=criteria, col="gray", lty=2)
      text(criteria, ylim[1], names(criteria), pos=3)
    }
    vnames <- colnames(vif)
    text(labx, laby, colnames(vif), pos=c(2,4)[1+(labels=="right")], xpd=TRUE)
  }
