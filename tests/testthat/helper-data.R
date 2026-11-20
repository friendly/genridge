# Shared test data: mtcars, 3 responses, 5 collinear predictors (max VIF about 18).
# Base R only, so the tests run without the Suggests packages.

mt_ynames <- c("mpg", "qsec", "drat")
mt_xnames <- c("cyl", "disp", "hp", "wt", "carb")
mt_Y <- as.matrix(mtcars[, mt_ynames])
mt_X <- as.matrix(mtcars[, mt_xnames])
mt_lambda <- c(0, 0.5, 1, 2, 5, 10, 20)
mt_formula <- cbind(mpg, qsec, drat) ~ cyl + disp + hp + wt + carb

mt_fit <- function(...) ridge(mt_Y, mt_X, lambda = mt_lambda, ...)

# centered X, scaled to unit mean square (divisor n), as ridge() does
std_X <- function(X) {
  Xc <- scale(X, scale = FALSE)
  Xc / rep(sqrt(colMeans(Xc^2)), each = nrow(X))
}
