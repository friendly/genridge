# Methods for "mlmridge" objects: extractors, as.ridge(), and the plot/summary methods

test_that("coef(): list, 3-d array, and one response", {
  m <- mt_fit()
  expect_identical(coef(m), m$coef)
  arr <- coef(m, type = "array")
  expect_identical(dim(arr), c(length(mt_lambda), length(mt_xnames), length(mt_ynames)))
  expect_identical(names(dimnames(arr)), c("lambda", "predictor", "response"))
  for (i in seq_along(mt_lambda))
    expect_equal(arr[i, , ], m$coef[[i]], ignore_attr = TRUE)
  expect_equal(coef(m, response = "qsec"), coef(ridge(mt_Y[, "qsec"], mt_X, lambda = mt_lambda)))
})

test_that("vcov(): pq x pq matrices with response:predictor names, or one response", {
  m <- mt_fit()
  V <- vcov(m)
  expect_length(V, length(mt_lambda))
  expect_identical(rownames(V[[1]])[1:2], c("mpg:cyl", "mpg:disp"))
  expect_equal(vcov(m, response = "drat"), vcov(ridge(mt_Y[, "drat"], mt_X, lambda = mt_lambda)))
})

test_that("as.ridge() rejects an unknown response, and is the identity on a 'ridge' object", {
  m <- mt_fit()
  expect_error(as.ridge(m, "nonsense"), "response must be one of")
  expect_error(as.ridge(m, 4), "response must be one of")
  expect_error(as.ridge(m, c(1, 2)), "response must be one of")
  u <- ridge(mt_Y[, 1], mt_X, lambda = mt_lambda)
  expect_identical(as.ridge(u), u)
})

test_that("print() shows each response and the pooled GCV, and returns its argument", {
  m <- mt_fit()
  out <- capture.output(res <- print(m))
  expect_identical(res, m)
  for (k in mt_ynames) expect_true(any(grepl(paste0("Ridge Coefficients for ", k), out)))
  expect_true(any(grepl("GCV pooled over responses \\(scaled\\)", out)))
})

test_that("summary methods delegate to the univariate methods for one response", {
  m <- mt_fit()
  a2 <- as.ridge(m, 2)
  expect_equal(precision(m, response = 2), precision(a2))
  expect_equal(pca(m, response = 2), pca(a2))
  # VIFs depend only on X: the same for every response, with predictor and lambda names
  v1 <- vif(m)
  expect_equal(v1$vif, vif(m, response = 3)$vif)
  expect_equal(v1$vif, vif(ridge(mt_Y[, 1], mt_X, lambda = mt_lambda))$vif)
  expect_identical(colnames(v1$vif), mt_xnames)
})

test_that("methods that use one response say which, when it isn't given", {
  withr::local_pdf(NULL)
  m <- mt_fit()
  expect_message(precision(m), "Using response 'mpg'")
  expect_message(plot(m), "Using response 'mpg'")
  expect_silent(precision(m, response = "mpg"))
})

test_that("traceplot() draws all responses and restores the graphics layout", {
  withr::local_pdf(NULL)
  m <- mt_fit()
  op <- par("mfrow")
  expect_no_error(traceplot(m))
  expect_identical(par("mfrow"), op)
  expect_no_error(traceplot(m, X = "df"))
  expect_no_error(traceplot(m, GCV = FALSE))
  expect_no_error(traceplot(m, response = c("qsec", "mpg")))
  expect_no_error(traceplot(m, response = 2))
  expect_error(traceplot(m, response = "nonsense"), "response must be in")
})

test_that("plot(predictor = ): ellipse centres are the coefficients for two responses", {
  withr::local_pdf(NULL)
  m <- mt_fit()
  res <- plot(m, predictor = "wt", responses = c("mpg", "qsec"))
  ctr <- res$centers$wt
  expect_identical(dim(ctr), c(length(mt_lambda), 2L))
  for (i in seq_along(mt_lambda))
    expect_equal(ctr[i, ], m$coef[[i]]["wt", c("mpg", "qsec")], ignore_attr = TRUE)

  # every ellipse has the shape of Sigma (Kronecker structure): the ellipse for lambda i,
  # centred and divided by sqrt(V_lambda[j, j]), is the same for all lambda
  shape <- function(i) {
    e <- res$ellipses$wt[[i]]
    sweep(e, 2, ctr[i, ]) / sqrt(m$cov[[i]]["wt", "wt"])
  }
  expect_equal(shape(1), shape(length(mt_lambda)))

  # several predictors, chosen by number; 'response' is ignored, with a warning
  res2 <- plot(m, predictor = 1:3)
  expect_identical(names(res2$centers), mt_xnames[1:3])
  expect_warning(plot(m, predictor = "wt", response = 1), "'response' is ignored")
  expect_error(plot(m, predictor = "nonsense"), "'predictor' must be among")
  expect_error(plot(m, predictor = "wt", responses = "mpg"), "'responses' must be two of")
})

test_that("other plot methods run for one response", {
  withr::local_pdf(NULL)
  m <- mt_fit()
  expect_no_error(plot(m, response = "qsec", variables = c("disp", "wt")))
  expect_no_error(pairs(m, response = "qsec"))
  expect_no_error(capture.output(biplot(m, response = "qsec")))
})
