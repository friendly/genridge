#' Extract a Univariate Ridge Regression from a Multivariate One
#'
#' @description
#' `ridge()` with a matrix response (e.g., `ridge(cbind(Y1, Y2) ~ x1 + x2, data)`) returns an
#' object of class `"mlmridge"`, fitting all responses with a common set of ridge constants.
#' `as.ridge()` extracts the results for one response as an ordinary `"ridge"` object, so that
#' all the methods for `"ridge"` objects ([traceplot()], [plot.ridge()], [pairs.ridge()],
#' [plot3d.ridge()], [precision()], [vif.ridge()], [pca()], ...) can be applied to it.
#'
#' The result is the same as fitting `ridge()` to that response alone, with the same
#' `lambda` or `df`.
#'
#' @param x An object of class `"mlmridge"`
#' @param response The response to extract: a name, or a number in `1:q`. The default is the
#'        first response.
#' @param ... Not used
#' @return An object of class `"ridge"`; see [ridge()] for its components.
#' @seealso [ridge()]
#' @export
#' @examples
#' if (requireNamespace("heplots", quietly = TRUE)) {
#'   data(schooldata, package = "heplots")
#'   sridge <- ridge(cbind(reading, mathematics, selfesteem) ~
#'                   education + occupation + visit + counseling + teacher,
#'                   data = schooldata, lambda = c(0, 1, 2, 5, 10, 20))
#'   # results for one response
#'   read.ridge <- as.ridge(sridge, "reading")
#'   read.ridge
#'   traceplot(read.ridge)
#' }
as.ridge <- function(x, ...) UseMethod("as.ridge")

#' @rdname as.ridge
#' @exportS3Method
as.ridge.ridge <- function(x, ...) x

#' @rdname as.ridge
#' @exportS3Method
as.ridge.mlmridge <- function(x, response = 1, ...) {
  ynames <- colnames(x$Sigma)
  k <- if (is.character(response)) match(response, ynames) else response
  if (length(k) != 1 || is.na(k) || k < 1 || k > length(ynames))
    stop("response must be one of: ", paste(ynames, collapse = ", "), ", or a number in 1:",
         length(ynames))

  coef <- t(sapply(x$coef, function(B) B[, k]))
  if (nrow(x$coef[[1]]) == 1) coef <- t(coef)       # single predictor
  dimnames(coef) <- list(format(x$lambda), rownames(x$coef[[1]]))
  s2 <- x$Sigma[k, k]
  cov <- lapply(x$cov, function(V) s2 * V)
  names(cov) <- NULL
  mse <- x$mse[, k]
  names(mse) <- NULL
  GCV <- x$GCV[, k]

  # same components, in the same order, as ridge.default()
  result <- list(lambda = x$lambda, df = x$df,
                 coef = coef, cov = cov,
                 mse = mse, scales = x$scales,
                 kHKB = x$kHKB[[k]], kLW = x$kLW[[k]],
                 GCV = GCV, kGCV = x$kGCV[[k]],
                 criteria = x$criteria[k, ])
  if (!is.null(x$svd.D))
    result <- c(result, list(svd.D = x$svd.D, svd.U = x$svd.U, svd.V = x$svd.V))
  class(result) <- "ridge"
  result
}
