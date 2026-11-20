# Fitting a multivariate ridge regression: ridge() with a matrix response -> class "mlmridge"

test_that("a matrix response gives an 'mlmridge' object with the documented structure", {
  m <- mt_fit()
  p <- length(mt_xnames); q <- length(mt_ynames); nl <- length(mt_lambda)
  expect_identical(class(m), "mlmridge")
  expect_type(m$coef, "list")
  expect_length(m$coef, nl)
  expect_identical(dim(m$coef[[1]]), c(p, q))
  expect_identical(dimnames(m$coef[[1]]), list(mt_xnames, mt_ynames))
  expect_length(m$cov, nl)
  expect_identical(dim(m$cov[[1]]), c(p, p))
  expect_identical(dim(m$Sigma), c(q, q))
  expect_identical(dim(m$mse), c(nl, q))
  expect_identical(dim(m$GCV), c(nl, q))
  expect_identical(dim(m$criteria), c(q, 3L))
  expect_length(m$GCV.pooled, nl)
  expect_identical(m$pool, "scaled")
})

test_that("a one-column matrix response still gives a univariate 'ridge' object", {
  # regression test: several lambdas with a one-column matrix failed ("non-conformable arrays")
  r <- ridge(mt_Y[, 1, drop = FALSE], mt_X, lambda = mt_lambda)
  expect_identical(class(r), "ridge")
  expect_equal(r, ridge(mt_Y[, 1], mt_X, lambda = mt_lambda))
  expect_equal(coef(ridge(cbind(mpg) ~ cyl + disp + hp + wt + carb, data = mtcars, lambda = mt_lambda)),
               coef(ridge(mpg ~ cyl + disp + hp + wt + carb, data = mtcars, lambda = mt_lambda)))
})

test_that("each response equals the univariate ridge() for that response", {
  m <- mt_fit()
  for (k in mt_ynames) {
    u <- ridge(mt_Y[, k], mt_X, lambda = mt_lambda)
    a <- as.ridge(m, k)
    expect_equal(a, u)                               # all components: coef, cov, mse, GCV, criteria, svd
    expect_identical(names(a), names(u))
  }
  # by number as well as by name
  expect_equal(as.ridge(m, 2), as.ridge(m, "qsec"))
})

test_that("the df parameterization is shared by all responses", {
  df <- c(5, 4.5, 4, 3.5, 3)
  m <- ridge(mt_Y, mt_X, df = df)
  expect_equal(m$df, df)
  expect_equal(as.ridge(m, "qsec"), ridge(mt_Y[, "qsec"], mt_X, df = df))
})

test_that("coefficients equal the closed-form estimator (X'X + lambda I)^-1 X'Y", {
  # This is VARshrink::lm_multiv_ridge()'s estimator with its lambda = genridge's lambda / n,
  # for X scaled as genridge does (Breiman & Friedman, 1997, eq. 6.3).
  m <- mt_fit()
  Xs <- std_X(mt_X)
  Yc <- scale(mt_Y, scale = FALSE)
  n <- nrow(Xs)
  for (i in seq_along(mt_lambda)) {
    B <- solve(crossprod(Xs) + mt_lambda[i] * diag(ncol(Xs)), crossprod(Xs, Yc))
    expect_equal(unname(m$coef[[i]]), unname(B))
    lambda_v <- mt_lambda[i] / n                     # VARshrink's scale
    Bv <- solve(crossprod(Xs) + n * lambda_v * diag(ncol(Xs)), crossprod(Xs, Yc))
    expect_equal(unname(m$coef[[i]]), unname(Bv))
  }
})

test_that("Sigma is the residual covariance matrix of the OLS multivariate fit", {
  m <- mt_fit()
  fit <- lm(mt_Y ~ mt_X)
  expect_equal(unname(m$Sigma), unname(crossprod(residuals(fit)) / fit$df.residual))
  expect_identical(dimnames(m$Sigma), list(mt_ynames, mt_ynames))
})

test_that("cov holds V_lambda, and Var(vec B) = Sigma %x% V_lambda", {
  m <- mt_fit()
  Xs <- std_X(mt_X)
  p <- ncol(Xs)
  i <- 4
  XPXI <- solve(crossprod(Xs) + mt_lambda[i] * diag(p))
  expect_equal(unname(m$cov[[i]]), unname(XPXI %*% crossprod(Xs) %*% XPXI))

  V <- vcov(m)[[i]]
  expect_identical(dim(V), c(p * 3L, p * 3L))
  expect_equal(unname(V), unname(kronecker(m$Sigma, m$cov[[i]])))
  # diagonal blocks are the univariate covariance matrices ...
  for (k in 1:3) {
    idx <- (k - 1) * p + seq_len(p)
    expect_equal(unname(V[idx, idx]),
                 unname(ridge(mt_Y[, k], mt_X, lambda = mt_lambda)$cov[[i]]))
  }
  # ... and an off-diagonal block is Sigma[k, l] * V_lambda
  expect_equal(unname(V[1:p, p + 1:p]), unname(m$Sigma[1, 2] * m$cov[[i]]))
})

test_that("the formula, mlm and matrix interfaces agree", {
  m  <- mt_fit()
  mf <- ridge(mt_formula, data = mtcars, lambda = mt_lambda)
  mm <- ridge(lm(mt_formula, data = mtcars), lambda = mt_lambda)
  strip <- function(z) { z$svd.U <- NULL; z }      # svd.U rownames: data rownames vs NULL
  expect_equal(strip(mf), strip(m))
  expect_equal(mm, mf)
  # pool = is passed through the formula and lm methods
  expect_identical(ridge(mt_formula, data = mtcars, lambda = mt_lambda, pool = "raw")$pool, "raw")
  expect_identical(ridge(lm(mt_formula, data = mtcars), lambda = mt_lambda, pool = "raw")$pool, "raw")
})

test_that("a model with a single predictor can be fitted (formula, matrix and lm interfaces)", {
  # regression test: ridge.formula() used X[, -Inter] without drop = FALSE
  u1 <- ridge(mpg ~ wt, data = mtcars, lambda = mt_lambda)
  u2 <- ridge(mtcars$mpg, as.matrix(mtcars["wt"]), lambda = mt_lambda)
  u3 <- ridge(lm(mpg ~ wt, data = mtcars), lambda = mt_lambda)
  expect_equal(coef(u1), coef(u2))
  expect_equal(coef(u1), coef(u3))
  m1 <- ridge(cbind(mpg, qsec) ~ wt, data = mtcars, lambda = mt_lambda)
  expect_s3_class(m1, "mlmridge")
  expect_equal(coef(m1, response = "mpg"), coef(u1))
  expect_identical(dim(vif(u1)$vif), c(length(mt_lambda), 1L))
  expect_equal(unlist(vif(u1)$vif), rep(1, length(mt_lambda)), ignore_attr = TRUE)
})
