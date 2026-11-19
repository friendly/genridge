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

# 8. coef(), vcov(), print() methods
check("coef(m) is the list", identical(coef(m), m$coef))
check("coef(m, response = k) == univariate coef",
      all.equal(coef(m, response = 2), coef(ridge(Y[, 2], X, lambda = lambda))))
arr <- coef(m, type = "array")
check("coef(type = 'array'): dim lambda x p x q", identical(dim(arr), c(6L, 5L, 3L)))
check("  array[i, , k] == coef list[[i]][, k]",
      all(sapply(seq_along(lambda), function(i) isTRUE(all.equal(arr[i, , ], m$coef[[i]], check.attributes = FALSE)))))
check("  array[, , k] == coef(response = k)",
      all.equal(unname(arr[, , "selfesteem"]), unname(coef(m, response = "selfesteem"))))
vc <- vcov(m)
check("vcov(m): list of pq x pq", is.list(vc) && identical(dim(vc[[1]]), c(15L, 15L)))
ur <- ridge(Y[, "mathematics"], X, lambda = lambda)
check("  diagonal block k of vcov == univariate cov (math, lambda[5])",
      all.equal(unname(vc[[5]][6:10, 6:10]), unname(vcov(ur)[[5]])))
check("  off-diagonal block == Sigma[j,k] * cov",
      all.equal(unname(vc[[3]][1:5, 11:15]), unname(m$Sigma[1, 3] * m$cov[[3]])))
check("vcov(m, response = k) == univariate vcov", all.equal(vcov(m, response = 2), vcov(ur)))
pr <- capture.output(print(m))
check("print() runs and shows pooled GCV", any(grepl("pooled", pr)))

# 9. errors
check("as.ridge() with bad response name -> error",
      inherits(tryCatch(as.ridge(m, "nope"), error = function(e) e), "error"))

# 10. traceplot(): all responses by default
pdf(NULL)
op_before <- par("mfrow")
msg <- tryCatch({ withCallingHandlers(traceplot(m), message = function(e) stop("message: ", conditionMessage(e))); "ok" },
                error = function(e) conditionMessage(e))
check("traceplot(m): all responses, no message, no error", identical(msg, "ok"))
check("  par(mfrow) restored afterwards", identical(par("mfrow"), op_before))
check("  draws q = 3 panels (one plot.new per response)", {
  n <- 0; setHook("plot.new", function() n <<- n + 1)
  traceplot(m); setHook("plot.new", NULL, "replace"); n == 3 })
par(mfrow = c(2, 2))
traceplot(m, response = c("reading", "selfesteem"), layout = FALSE)
check("  layout = FALSE keeps the user's layout", identical(par("mfrow"), c(2L, 2L)))
par(mfrow = c(1, 1))
check("traceplot(m, response = 2:3) runs", is.null(traceplot(m, response = 2:3)))
check("traceplot(m, response = 'nope') -> error",
      inherits(tryCatch(traceplot(m, response = "nope"), error = function(e) e), "error"))
dev.off()

# 11. plot(predictor = ): cross-response coefficient ellipses
pdf(NULL)
cr <- plot(m, predictor = "visit", responses = c("reading", "selfesteem"))
check("plot(predictor=): centres == coef rows across lambda",
      all.equal(unname(cr$centers$visit), unname(t(sapply(m$coef, function(B) B["visit", c(1, 3)])))))
# ellipse at lambda[4]: points satisfy (z - c)' S^-1 (z - c) = radius^2 with S = V_jj * Sigma[r, r]
S4 <- m$cov[[4]]["visit", "visit"] * m$Sigma[c(1, 3), c(1, 3)]
z <- sweep(cr$ellipses$visit[[4]], 2, cr$centers$visit[4, ])
check("  ellipse shape == [V_lambda]_jj * Sigma[r, r] (radius 1)",
      all.equal(rowSums((z %*% solve(S4)) * z), rep(1, nrow(z))))
cr2 <- plot(m, predictor = 1:5, which.lambda = c(1, 6), responses = 1:2)
check("  several predictors by number, which.lambda subset",
      identical(names(cr2$centers), rownames(m$coef[[1]])) && nrow(cr2$centers[[1]]) == 2)
check("  fill = TRUE, vector fill.alpha runs",
      is.list(plot(m, predictor = "visit", fill = TRUE, fill.alpha = c(0.1, 0.3))))
check("  bad predictor -> error",
      inherits(tryCatch(plot(m, predictor = "nope"), error = function(e) e), "error"))
check("  bad responses -> error",
      inherits(tryCatch(plot(m, predictor = 1, responses = c("reading", "x")), error = function(e) e), "error"))
w <- tryCatch(plot(m, predictor = 1, response = 2), warning = function(w) conditionMessage(w))
check("  response = with predictor = -> warning", is.character(w) && grepl("ignored", w))
check("plot(m, response = 2) still the within-response plot",
      is.null(suppressMessages(plot(m, response = 2))))
dev.off()

cat(if (ok) "\nALL CHECKS PASSED\n" else "\nSOME CHECKS FAILED\n")
