# lm method for ridge regression

#' @rdname ridge
#' @exportS3Method 
ridge.lm <-
  function(formula, data, lambda=0, df, svd=TRUE, ...){
    
    #code from MASS:::lm.ridge
    m <- match.call(expand.dots = FALSE)
    m$model <- m$x <- m$y <- m$contrasts <- m$... <- m$lambda <-m$df <- NULL
    m[[1L]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    Y <- model.response(m)
    X <- model.matrix(Terms, m, contrasts)
    n <- nrow(X)
    p <- ncol(X)
    offset <- model.offset(m)
    if (!is.null(offset)) 
      Y <- Y - offset
    if (Inter <- attr(Terms, "intercept")) {
      Xm <- colMeans(X[, -Inter])
      Ym <- mean(Y)
      p <- p - 1
      X <- X[, -Inter] - rep(Xm, rep(n, p))
      Y <- Y - Ym
    }
    ridge.default(Y, X, lambda=lambda, df=df, svd=svd)
  }


if(FALSE) {
longley.lm <- lm(Employed ~ GNP + Unemployed + Armed.Forces + 
                   Population + Year + GNP.deflator,		data=longley)

ridge(longley.lm)

mf <- model.frame(longley.lm)
names(mf)

y <- mf[, 1]
X <- mf[, -1]

frame <- model.frame(longley.lm)
y <- model.response(longley.lm)
X <- model.matrix(longley.lm)

}