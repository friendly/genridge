# The univariate ridge() and its methods must be unchanged by the multivariate extension.
# fixtures/ridge-ref-0.8.1.rds was saved by dev/make-ridge-ref.R (genridge 0.8.1, before
# ridge.default() was generalized to a matrix response); the fits below repeat that script.

test_that("univariate fits and methods match the saved 0.8.1 reference results", {
  ref <- readRDS(test_path("fixtures", "ridge-ref-0.8.1.rds"))
  prostate <- genridge::prostate; Manpower <- genridge::Manpower; Acetylene <- genridge::Acetylene

  lambda <- c(0, 0.005, 0.01, 0.02, 0.04, 0.08)
  lform <- Employed ~ GNP + Unemployed + Armed.Forces + Population + Year + GNP.deflator
  fits <- list(
    longley_formula = ridge(lform, data = longley, lambda = lambda),
    longley_default = ridge(longley[, "Employed"], data.matrix(longley[, c(2:6, 1)]), lambda = lambda),
    longley_lm      = ridge(lm(lform, data = longley), lambda = lambda),
    longley_nosvd   = ridge(longley[, "Employed"], data.matrix(longley[, c(2:6, 1)]), lambda = lambda, svd = FALSE),
    longley_formula_nosvd = ridge(lform, data = longley, lambda = lambda, svd = FALSE),
    prostate_df     = ridge(prostate[, "lpsa"], data.matrix(prostate[, 1:8]), df = 8:1),
    manpower_df     = ridge(Hours ~ ., data = Manpower, df = seq(5, 3.75, -.25)),
    acetylene       = ridge(yield ~ temp + ratio + time + I(time^2) + temp:time, data = Acetylene,
                            lambda = c(0, 0.0005, 0.001, 0.002, 0.005, 0.01))
  )
  expect_identical(names(fits), names(ref$fits))

  for (nm in names(fits)) {
    f <- fits[[nm]]
    expect_equal(f, ref$fits[[nm]], label = nm)
    r <- ref$methods[[nm]]
    expect_equal(coef(f), r$coef, label = paste(nm, "coef"))
    expect_equal(vcov(f), r$vcov, label = paste(nm, "vcov"))
    expect_equal(best(f), r$best, label = paste(nm, "best"))
    expect_equal(vif(f), r$vif, label = paste(nm, "vif"))
    expect_equal(precision(f), r$precision, label = paste(nm, "precision"))
    if (!is.null(f$svd.V)) expect_equal(pca(f), r$pca, label = paste(nm, "pca"))
  }
})
