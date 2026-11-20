# Criteria for choosing lambda in a multivariate ridge regression: per response, pooled GCV, best()

test_that("per-response criteria equal the univariate ones", {
  m <- mt_fit()
  for (k in seq_along(mt_ynames)) {
    u <- ridge(mt_Y[, k], mt_X, lambda = mt_lambda)
    expect_equal(m$kHKB[[k]], u$kHKB)
    expect_equal(m$kLW[[k]], u$kLW)
    expect_equal(unname(m$GCV[, k]), unname(u$GCV))
    expect_equal(m$kGCV[[k]], u$kGCV)
    expect_equal(m$criteria[k, ], u$criteria)
  }
})

test_that("pooled GCV: 'scaled', 'raw' and 'mahalanobis' follow their definitions", {
  m  <- mt_fit()                                     # pool = "scaled", the default
  mr <- mt_fit(pool = "raw")
  mh <- mt_fit(pool = "mahalanobis")

  # scaled: sum over responses of GCV_k / sigma^2_kk
  expect_equal(unname(m$GCV.pooled), unname(rowSums(sweep(m$GCV, 2, diag(m$Sigma), "/"))))
  # raw: sum over responses of GCV_k
  expect_equal(unname(mr$GCV.pooled), unname(rowSums(m$GCV)))
  # mahalanobis: tr(Sigma^-1 R'R) / (n - df)^2
  n <- nrow(mt_X)
  Xs <- std_X(mt_X); Yc <- scale(mt_Y, scale = FALSE)
  Sinv <- solve(m$Sigma)
  gcv.m <- sapply(seq_along(mt_lambda), function(i) {
    R <- Yc - Xs %*% m$coef[[i]]
    sum(diag(Sinv %*% crossprod(R))) / (n - m$df[i])^2
  })
  expect_equal(unname(mh$GCV.pooled), gcv.m)

  # kGCV.pooled is the lambda minimizing GCV.pooled; the estimates don't depend on pool
  for (f in list(m, mr, mh))
    expect_equal(f$kGCV.pooled, mt_lambda[which.min(f$GCV.pooled)])
  expect_equal(mr$coef, m$coef)
  expect_error(mt_fit(pool = "nonsense"))
})

test_that("pool = 'scaled' does not depend on the units of the responses", {
  Y2 <- mt_Y %*% diag(c(1, 100, 0.01))
  colnames(Y2) <- mt_ynames
  m  <- mt_fit()
  m2 <- ridge(Y2, mt_X, lambda = mt_lambda)
  expect_equal(m2$GCV.pooled, m$GCV.pooled)
  expect_equal(m2$kGCV.pooled, m$kGCV.pooled)
})

test_that("best() for 'mlmridge': one row per response and criterion, plus pooled GCV", {
  m <- mt_fit()
  b <- best(m)
  q <- length(mt_ynames)
  expect_s3_class(b, "data.frame")
  expect_identical(names(b), c("response", "criterion", "k", "df"))
  expect_identical(nrow(b), 3L * q + 1L)
  expect_identical(levels(b$response), c(mt_ynames, "(pooled)"))
  expect_identical(levels(b$criterion), c("HKB", "LW", "GCV"))

  # each response's rows equal best() for the univariate fit
  for (k in mt_ynames) {
    u <- best(ridge(mt_Y[, k], mt_X, lambda = mt_lambda))
    bk <- b[b$response == k, ]
    expect_identical(as.character(bk$criterion), rownames(u))
    expect_equal(bk$k, u$k)
    expect_equal(bk$df, u$df)
  }
  # the pooled row
  last <- b[nrow(b), ]
  expect_identical(as.character(last$response), "(pooled)")
  expect_identical(as.character(last$criterion), "GCV")
  expect_equal(last$k, m$kGCV.pooled)
  expect_equal(last$df, m$df[match(m$kGCV.pooled, mt_lambda)])

  # response = gives the univariate table
  expect_identical(best(m, response = "qsec"), best(as.ridge(m, "qsec")))
  expect_identical(best(m, response = 2), best(m, response = "qsec"))
  expect_error(best(m, response = "nonsense"), "response must be one of")
})

test_that("best() for 'mlmridge' has no df column when svd = FALSE", {
  b <- best(mt_fit(svd = FALSE))
  expect_identical(names(b), c("response", "criterion", "k"))
})
