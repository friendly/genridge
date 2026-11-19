# Checks for the multivariate branch of ridge.default() (class "mlmridge").
# Run from the package root:  Rscript dev/check-mlmridge.R
# (The univariate reference results are checked separately by dev/check-ridge-ref.R.)

suppressMessages(devtools::load_all(".", quiet = TRUE))
data(schooldata, package = "heplots")

ok <- TRUE
check <- function(label, cond) {
  cat(sprintf("%-62s %s\n", label, if (isTRUE(cond)) "OK" else "FAIL"))
  if (!isTRUE(cond)) { ok <<- FALSE; if (!is.logical(cond)) print(cond) }
}

ynames <- c("reading", "mathematics", "selfesteem")
xnames <- c("education", "occupation", "visit", "counseling", "teacher")
Y <- as.matrix(schooldata[, ynames])
X <- as.matrix(schooldata[, xnames])
lambda <- c(0, 1, 2, 5, 10, 20)

m <- ridge(Y, X, lambda = lambda)
check("class is 'mlmridge'", identical(class(m), "mlmridge"))
check("coef: list of p x q matrices, one per lambda",
      is.list(m$coef) && length(m$coef) == length(lambda) && all(dim(m$coef[[1]]) == c(5, 3)))

# 1. each response equals the univariate ridge() on that response
for (k in ynames) {
  u <- ridge(Y[, k], X, lambda = lambda)
  a <- as.ridge(m, k)
  check(sprintf("as.ridge(m, '%s') all.equal univariate ridge()", k), all.equal(a, u))
  check(sprintf("  ... same components in same order (%s)", k), identical(names(a), names(u)))
}

# 2. df parameterization
mdf <- ridge(Y, X, df = c(5, 4.5, 4, 3.5))
check("df =: each response equals univariate",
      all.equal(as.ridge(mdf, 2), ridge(Y[, 2], X, df = c(5, 4.5, 4, 3.5))))

# 3. the three interfaces agree
f <- cbind(reading, mathematics, selfesteem) ~ education + occupation + visit + counseling + teacher
mf <- ridge(f, data = schooldata, lambda = lambda)
mm <- ridge(lm(f, data = schooldata), lambda = lambda)
strip <- function(z) { z$svd.U <- NULL; z }     # svd.U rownames differ (data rownames vs NULL)
check("formula interface == matrix interface", all.equal(strip(mf), strip(m)))
check("mlm interface == formula interface", all.equal(mm, mf))

# 4. pooled GCV
check("pool = 'scaled' == sum_k GCV_k / MSE_k",
      all.equal(unname(m$GCV.pooled), unname(rowSums(sweep(m$GCV, 2, diag(m$Sigma), "/")))))
mr <- ridge(Y, X, lambda = lambda, pool = "raw")
check("pool = 'raw' == sum_k GCV_k", all.equal(unname(mr$GCV.pooled), unname(rowSums(m$GCV))))
mh <- ridge(Y, X, lambda = lambda, pool = "mahalanobis")
check("pool = 'mahalanobis' runs, finite", all(is.finite(mh$GCV.pooled)))
check("pool passed through formula interface",
      identical(ridge(f, data = schooldata, lambda = lambda, pool = "raw")$pool, "raw"))
cat("\nkGCV per response:", format(m$kGCV), "; pooled:", m$kGCV.pooled,
    "(scaled)", mr$kGCV.pooled, "(raw)", mh$kGCV.pooled, "(mahalanobis)\n\n")

# 5. covariance structure: Sigma[k,k] * cov == univariate cov
check("Sigma[k,k] * cov[[i]] == univariate cov[[i]]",
      all.equal(m$Sigma[3, 3] * m$cov[[4]], ridge(Y[, 3], X, lambda = lambda)$cov[[4]]))

# 6. agreement with VARshrink's estimator (lambda_v = lambda_g / n, X scaled as genridge does)
source("dev/ridge_mlm.R")
n <- nrow(X)
Xc <- scale(X, scale = FALSE); Xg <- Xc / rep(sqrt(colMeans(Xc^2)), each = n)
v <- ridge_mlm(scale(Y, scale = FALSE), Xg, lambda = lambda / n)
check("coef equal VARshrink ridge_mlm (lambda_v = lambda_g / n)",
      all.equal(lapply(m$coef, unname), lapply(unname(v$B), unname), check.attributes = FALSE))

# 7. a 1-column matrix y still gives a univariate "ridge" object
check("1-column matrix y -> class 'ridge'", identical(class(ridge(Y[, 1, drop = FALSE], X)), "ridge"))

# 8. errors
check("as.ridge() with bad response name -> error",
      inherits(tryCatch(as.ridge(m, "nope"), error = function(e) e), "error"))

cat(if (ok) "\nALL CHECKS PASSED\n" else "\nSOME CHECKS FAILED\n")
