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


#' Methods for Multivariate Ridge Regression Objects
#'
#' @description
#' Basic methods for objects of class `"mlmridge"`, returned by [ridge()] when the response is
#' a matrix.
#'
#' * `coef()` extracts the coefficients: by default a list with one \eqn{p \times q} matrix for
#'   each ridge constant; with `type = "array"`, a 3-d array indexed by `[lambda, predictor, response]`;
#'   with `response`, the \eqn{\lambda \times p} matrix for one response, as for a `"ridge"` object.
#' * `vcov()` extracts the covariance matrices of the coefficients: by default a list with, for each
#'   ridge constant, the \eqn{pq \times pq} covariance matrix of the stacked coefficients
#'   \eqn{\mathrm{vec}(\mathbf{B})}, which is `kronecker(Sigma, cov[[i]])`; with `response`, the list of
#'   \eqn{p \times p} matrices for one response.
#' * `print()` shows the coefficients for each response and the optimal ridge constants given by
#'   the HKB, LW and GCV criteria, for each response and pooled over responses.
#'
#' @name mlmridge-methods
#' @param object,x An object of class `"mlmridge"`
#' @param response If not `NULL`, a response (name or number) for which to return the results,
#'        in the same form as for a univariate `"ridge"` object; see [as.ridge()].
#' @param type For `coef()` with `response = NULL`: `"list"` (the default) or `"array"`.
#' @param digits Number of significant digits to print.
#' @param ... Other arguments, not used.
#' @return `coef()` and `vcov()` return a list, array or matrix, as described above;
#'        `print()` returns `x`, invisibly.
#' @seealso [ridge()], [as.ridge()]
#' @examples
#' if (requireNamespace("heplots", quietly = TRUE)) {
#'   data(schooldata, package = "heplots")
#'   sridge <- ridge(cbind(reading, mathematics, selfesteem) ~
#'                   education + occupation + visit + counseling + teacher,
#'                   data = schooldata, lambda = c(0, 0.1, 0.5, 1, 2, 5))
#'   sridge
#'   coef(sridge)[[4]]                       # coefficients for lambda = 1
#'   coef(sridge, response = "reading")      # lambda x p, for one response
#'   dim(coef(sridge, type = "array"))
#'   dim(vcov(sridge)[[1]])                  # pq x pq
#' }
NULL

#' @rdname mlmridge-methods
#' @exportS3Method coef mlmridge
coef.mlmridge <- function(object, response = NULL, type = c("list", "array"), ...) {
  if (!is.null(response)) return(as.ridge(object, response)$coef)
  type <- match.arg(type)
  if (type == "list") return(object$coef)
  B <- object$coef
  arr <- array(unlist(lapply(B, function(b) t(b))),          # each slice: response x predictor
               dim = c(ncol(B[[1]]), nrow(B[[1]]), length(B)))
  arr <- aperm(arr, c(3, 2, 1))                              # -> lambda x predictor x response
  dimnames(arr) <- list(lambda = names(B), predictor = rownames(B[[1]]),
                        response = colnames(B[[1]]))
  arr
}

#' @rdname mlmridge-methods
#' @exportS3Method vcov mlmridge
vcov.mlmridge <- function(object, response = NULL, ...) {
  if (!is.null(response)) return(as.ridge(object, response)$cov)
  S <- object$Sigma
  xnames <- rownames(object$coef[[1]])
  vnames <- paste(rep(colnames(S), each = length(xnames)), xnames, sep = ":")
  lapply(object$cov, function(V) {
    K <- kronecker(S, V)
    dimnames(K) <- list(vnames, vnames)
    K
  })
}

#' @rdname mlmridge-methods
#' @exportS3Method print mlmridge
print.mlmridge <- function(x, digits = max(5, getOption("digits") - 5), ...) {
  ynames <- colnames(x$Sigma)
  B <- x$coef[[1]]
  cat("Multivariate ridge regression:", length(ynames), "responses,", nrow(B), "predictors,",
      length(x$lambda), "ridge constants\n")
  for (k in seq_along(ynames)) {
    cat("\nRidge Coefficients for ", ynames[k], ":\n", sep = "")
    print.default(format(coef(x, response = k), digits = digits), print.gap = 2, quote = FALSE)
  }
  cat("\nOptimal ridge constants:\n")
  print.default(format(x$criteria, digits = digits), print.gap = 2, quote = FALSE)
  cat("GCV pooled over responses (", x$pool, "): ", format(x$kGCV.pooled, digits = digits),
      "\n", sep = "")
  invisible(x)
}
