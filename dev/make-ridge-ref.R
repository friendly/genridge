# Reference results from the current (univariate) ridge() and its methods, saved before
# generalizing ridge.default() to a matrix response (see dev/ridge_mlm-plan.md, §9).
# After changing the code, run dev/check-ridge-ref.R: every result must be identical.

library(genridge)   # use devtools::load_all() instead when checking work in progress
stopifnot(packageVersion("genridge") >= "0.8.1")

data(prostate, Manpower, Acetylene, package = "genridge")

lambda <- c(0, 0.005, 0.01, 0.02, 0.04, 0.08)
lform <- Employed ~ GNP + Unemployed + Armed.Forces + Population + Year + GNP.deflator

fits <- list(
  longley_formula = ridge(lform, data = longley, lambda = lambda),
  longley_default = ridge(longley[, "Employed"], data.matrix(longley[, c(2:6, 1)]), lambda = lambda),
  longley_lm      = ridge(lm(lform, data = longley), lambda = lambda),
  # (default interface: ridge.formula(svd = FALSE) failed in 0.8.1, fixed afterwards)
  longley_nosvd   = ridge(longley[, "Employed"], data.matrix(longley[, c(2:6, 1)]), lambda = lambda, svd = FALSE),
  prostate_df     = ridge(prostate[, "lpsa"], data.matrix(prostate[, 1:8]), df = 8:1),
  manpower_df     = ridge(Hours ~ ., data = Manpower, df = seq(5, 3.75, -.25)),
  acetylene       = ridge(yield ~ temp + ratio + time + I(time^2) + temp:time, data = Acetylene,
                          lambda = c(0, 0.0005, 0.001, 0.002, 0.005, 0.01))
)

# downstream methods, on each fit
methods <- lapply(fits, function(f) list(
  coef      = coef(f),
  vcov      = vcov(f),
  best      = genridge:::best(f),   # best() generic is not exported in 0.8.1
  print     = capture.output(print(f)),
  vif       = vif(f),
  precision = precision(f),
  pca       = if (!is.null(f$svd.V)) pca(f) else NULL
))

ref <- list(version = as.character(packageVersion("genridge")),
            R = R.version.string, date = Sys.Date(),
            fits = fits, methods = methods)
saveRDS(ref, "dev/ridge-ref-0.8.1.rds")
cat("saved", length(fits), "fits to dev/ridge-ref-0.8.1.rds\n")
