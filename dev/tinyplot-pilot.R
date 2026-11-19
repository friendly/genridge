# Pilot: genridge plots drawn with tinyplot (see dev/tinyplot-notes.md).
# Run from the package root:  Rscript dev/tinyplot-pilot.R [outdir]
suppressMessages(devtools::load_all(".", quiet = TRUE))
library(tinyplot)
args <- commandArgs(TRUE)
outdir <- if (length(args)) args[1] else tempdir()

data(schooldata, package = "heplots")
lambda <- c(0, 0.1, 0.5, 1, 2, 5)
sridge <- ridge(cbind(reading, mathematics, selfesteem) ~
                education + occupation + visit + counseling + teacher,
                data = schooldata, lambda = lambda)

# ---- 1. traceplot for all responses: one tinyplot call, facets = responses, by = predictor ----
arr <- coef(sridge, type = "array")                      # lambda x predictor x response
long <- as.data.frame.table(arr, responseName = "coef")
long$lambda <- as.numeric(as.character(long$lambda))
png(file.path(outdir, "tp-traceplot.png"), width = 1200, height = 450, res = 110)
tinytheme("clean")
tinyplot(coef ~ lambda | predictor, facet = ~ response, data = long, type = "b", pch = 16,
         facet.args = list(free = TRUE), xlab = "Ridge constant", ylab = "Coefficient")
tinytheme()                                              # reset
dev.off()

# ---- 2. bivariate ridge trace (plot.ridge) for one response: ellipses as polygons, by = lambda ----
r <- as.ridge(sridge, "reading")
vars <- c("occupation", "visit")
ell <- do.call(rbind, lapply(seq_along(lambda), function(i) {
  ctr <- r$coef[i, vars]; S <- r$cov[[i]][vars, vars]
  a <- seq(0, 2 * pi, length.out = 61)
  xy <- cbind(cos(a), sin(a)) %*% chol(S) * 0.5 + rep(ctr, each = 61)   # radius 0.5, as in gallery
  data.frame(x = xy[, 1], y = xy[, 2], lambda = factor(lambda[i]))
}))
ctrs <- data.frame(x = r$coef[, vars[1]], y = r$coef[, vars[2]], lambda = factor(lambda))
png(file.path(outdir, "tp-plot.png"), width = 700, height = 600, res = 110)
tinytheme("clean")
tinyplot(y ~ x | lambda, data = ell, type = "polygon", fill = 0.15,
         xlab = vars[1], ylab = vars[2], main = "reading")
tinyplot_add(data = ctrs, type = "p", pch = 19)
lines(ctrs$x, ctrs$y)                                    # base graphics still works on top
tinytheme()
dev.off()
cat("wrote tp-traceplot.png and tp-plot.png to", outdir, "\n")
