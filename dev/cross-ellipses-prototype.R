# Prototype: cross-response coefficient ellipses for multivariate ridge ("mlmridge"),
# plan §5.3. Not yet in the package.
#
# For predictor j, row j of B(lambda) is the q-vector of its coefficients for all responses.
# Its covariance matrix is  Var(b_j.) = [V_lambda]_jj * Sigma  (Kronecker structure), so in the
# plane of two responses the ellipse for predictor j at lambda has
#   center = (B[j, r1], B[j, r2]),  shape = [V_lambda]_jj * Sigma[c(r1, r2), c(r1, r2)].
# The shape (orientation, eccentricity) is that of the residual covariance Sigma for every
# predictor and every lambda; only the centre and the size (sqrt([V_lambda]_jj)) change.
#
# Run from the package root:  Rscript dev/cross-ellipses-prototype.R [outdir]

suppressMessages(devtools::load_all(".", quiet = TRUE))
args <- commandArgs(TRUE)
outdir <- if (length(args)) args[1] else tempdir()

ellipse_xy <- function(center, shape, radius = 1, segments = 60) {
  a <- seq(0, 2 * pi, length.out = segments + 1)
  sweep(radius * cbind(cos(a), sin(a)) %*% chol(shape), 2, center, "+")
}

# coefficient ellipses for one or more predictors in the space of two responses, over lambda
coefplot_ridge <- function(x, predictor, responses = 1:2, which.lambda = seq_along(x$lambda),
                           radius = 1, col = NULL, fill = FALSE, fill.alpha = 0.2,
                           labels = TRUE, xlim = NULL, ylim = NULL, main = NULL, ...) {
  ynames <- colnames(x$Sigma)
  if (is.numeric(responses)) responses <- ynames[responses]
  S <- x$Sigma[responses, responses]
  lam <- x$lambda[which.lambda]
  nl <- length(lam)
  if (is.null(col)) col <- grDevices::hcl.colors(nl, "viridis", rev = TRUE)
  col <- rep_len(col, nl)

  # all ellipses, for the axis limits
  ells <- lapply(predictor, function(j) lapply(which.lambda, function(i) {
    ellipse_xy(x$coef[[i]][j, responses], x$cov[[i]][j, j] * S, radius)
  }))
  allxy <- do.call(rbind, unlist(ells, recursive = FALSE))
  if (is.null(xlim)) xlim <- range(allxy[, 1], 0)
  if (is.null(ylim)) ylim <- range(allxy[, 2], 0)
  if (is.null(main)) main <- if (length(predictor) == 1) predictor else ""

  plot(NA, xlim = xlim, ylim = ylim, xlab = paste("coefficient for", responses[1]),
       ylab = paste("coefficient for", responses[2]), main = main, ...)
  abline(h = 0, v = 0, col = "gray")
  for (p in seq_along(predictor)) {
    ctr <- t(sapply(which.lambda, function(i) x$coef[[i]][predictor[p], responses]))
    for (i in seq_len(nl))
      polygon(ells[[p]][[i]], border = col[i],
              col = if (fill) grDevices::adjustcolor(col[i], fill.alpha) else NA)
    lines(ctr, lty = 1)
    points(ctr, pch = 16, col = col)
    if (labels && length(predictor) > 1) text(ctr[1, , drop = FALSE], predictor[p], pos = 3, cex = 0.9)
  }
  invisible(NULL)
}

data(schooldata, package = "heplots")
lambda <- c(0, 0.1, 0.5, 1, 2, 5)
sridge <- ridge(cbind(reading, mathematics, selfesteem) ~
                education + occupation + visit + counseling + teacher,
                data = schooldata, lambda = lambda)
legend_lambda <- function(pos = "topleft")
  legend(pos, legend = lambda, title = expression(lambda), pch = 16, bty = "n", cex = 0.8,
         col = grDevices::hcl.colors(length(lambda), "viridis", rev = TRUE))

# A. one predictor per panel, all lambda: reading x mathematics
png(file.path(outdir, "ce-A-one-predictor.png"), width = 1200, height = 450, res = 110)
op <- par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
for (j in c("occupation", "visit", "counseling")) {
  coefplot_ridge(sridge, j, responses = c("reading", "mathematics"), radius = 1, fill = TRUE)
  if (j == "occupation") legend_lambda("bottomright")
}
par(op); dev.off()

# B. same, reading x selfesteem: lower residual correlation (0.47) -> rounder ellipses
png(file.path(outdir, "ce-B-other-pair.png"), width = 1200, height = 450, res = 110)
op <- par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
for (j in c("occupation", "visit", "counseling"))
  coefplot_ridge(sridge, j, responses = c("reading", "selfesteem"), radius = 1, fill = TRUE)
par(op); dev.off()

# C. all predictors in one panel, ellipses only at lambda = 0 and 5, paths between:
#    the ridge analogue of heplots::coefplot()
png(file.path(outdir, "ce-C-all-predictors.png"), width = 700, height = 650, res = 110)
coefplot_ridge(sridge, rownames(sridge$coef[[1]]), responses = c("reading", "mathematics"),
               which.lambda = c(1, 6), col = c("black", "red"), radius = 1,
               main = expression("All predictors, " * lambda * " = 0 (black) and 5 (red)"))
# full paths through all lambda values
for (j in rownames(sridge$coef[[1]]))
  lines(t(sapply(sridge$coef, function(B) B[j, c("reading", "mathematics")])), lty = 2, col = "gray40")
dev.off()
cat("wrote ce-*.png to", outdir, "\n")
