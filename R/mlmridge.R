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


#' Plot and Summary Methods for Multivariate Ridge Regression
#'
#' @description
#' These methods apply the plots and summaries for univariate ridge regression to one response
#' of a multivariate ridge regression (class `"mlmridge"`). Each extracts the results for that
#' response with [as.ridge()] and calls the corresponding method for `"ridge"` objects:
#'
#' * `traceplot()`: univariate ridge trace plot ([traceplot()])
#' * `plot()`: bivariate ridge trace plot of coefficients with covariance ellipses ([plot.ridge()])
#' * `pairs()`: scatterplot matrix of bivariate ridge trace plots ([pairs.ridge()])
#' * `plot3d()`: 3D ridge trace plot ([plot3d.ridge()])
#' * `biplot()`: ridge trace plot with variable vectors ([biplot.ridge()])
#' * `pca()`: transformation to PCA/SVD space ([pca.ridge()]), returning a `"pcaridge"` object
#' * `precision()`: measures of precision and shrinkage ([precision()])
#' * `vif()`: variance inflation factors ([vif.ridge()])
#'
#' Because all responses share the same ridge constants and predictors, the shrinkage factors,
#' effective degrees of freedom, PCA/SVD of the predictors and the variance inflation factors
#' are the same for every response; the coefficients, their covariance matrices (scaled by each
#' response's residual variance) and the HKB, LW and GCV criteria differ.
#'
#' @name mlmridge-plots
#' @param x,object,mod An object of class `"mlmridge"`
#' @param response The response to use: a name, or a number in `1:q`. If not given, the first
#'        response is used, with a message.
#' @param main Title for the plot. The default is the name of the response.
#' @param det.fun,normalize See [precision()].
#' @param ... Other arguments, passed to the method for `"ridge"` objects.
#' @return As for the corresponding method for `"ridge"` objects.
#' @seealso [ridge()], [as.ridge()], [mlmridge-methods]
#' @importFrom graphics pairs
#' @importFrom stats biplot
#' @examples
#' if (requireNamespace("heplots", quietly = TRUE)) {
#'   data(schooldata, package = "heplots")
#'   lambda <- c(0, 0.1, 0.5, 1, 2, 5)
#'   sridge <- ridge(cbind(reading, mathematics, selfesteem) ~
#'                   education + occupation + visit + counseling + teacher,
#'                   data = schooldata, lambda = lambda)
#'
#'   traceplot(sridge, response = "reading")
#'   traceplot(sridge, response = "selfesteem", X = "df")
#'
#'   # bivariate ridge trace plot for two predictors
#'   plot(sridge, response = "mathematics", variables = c("occupation", "visit"),
#'        radius = 0.5)
#'
#'   pairs(sridge, response = "reading", radius = 0.5)
#'
#'   vif(sridge)
#'   precision(sridge, response = "reading")
#' }
NULL

# select one response of an "mlmridge" object, as a "ridge" object; message if not specified
select_response <- function(x, response, missing_response) {
  ynames <- colnames(x$Sigma)
  r <- as.ridge(x, response)
  k <- if (is.character(response)) match(response, ynames) else response
  if (missing_response && length(ynames) > 1)
    message("Using response '", ynames[k], "' (1 of ", length(ynames),
            "); use response = to choose another")
  list(ridge = r, name = ynames[k])
}

#' @rdname mlmridge-plots
#' @exportS3Method traceplot mlmridge
traceplot.mlmridge <- function(x, response = 1, main, ...) {
  sel <- select_response(x, response, missing(response))
  if (missing(main)) main <- sel$name
  traceplot(sel$ridge, main = main, ...)
}

#' @rdname mlmridge-plots
#' @exportS3Method plot mlmridge
plot.mlmridge <- function(x, response = 1, main, ...) {
  sel <- select_response(x, response, missing(response))
  if (missing(main)) main <- sel$name
  plot(sel$ridge, main = main, ...)
}

#' @rdname mlmridge-plots
#' @exportS3Method pairs mlmridge
pairs.mlmridge <- function(x, response = 1, ...) {
  sel <- select_response(x, response, missing(response))
  pairs(sel$ridge, ...)
}

#' @rdname mlmridge-plots
#' @exportS3Method plot3d mlmridge
plot3d.mlmridge <- function(x, response = 1, ...) {
  sel <- select_response(x, response, missing(response))
  plot3d(sel$ridge, ...)
}

#' @rdname mlmridge-plots
#' @exportS3Method biplot mlmridge
biplot.mlmridge <- function(x, response = 1, main, ...) {
  sel <- select_response(x, response, missing(response))
  if (missing(main)) main <- sel$name
  biplot(sel$ridge, main = main, ...)
}

#' @rdname mlmridge-plots
#' @exportS3Method pca mlmridge
pca.mlmridge <- function(x, response = 1, ...) {
  sel <- select_response(x, response, missing(response))
  pca(sel$ridge, ...)
}

#' @rdname mlmridge-plots
#' @exportS3Method precision mlmridge
precision.mlmridge <- function(object, det.fun = c("log", "root"), normalize = TRUE,
                               response = 1, ...) {
  sel <- select_response(object, response, missing(response))
  precision(sel$ridge, det.fun = det.fun, normalize = normalize, ...)
}

#' @rdname mlmridge-plots
#' @exportS3Method vif mlmridge
vif.mlmridge <- function(mod, response = 1, ...) {
  # VIFs depend only on X, so they are the same for every response; the criteria
  # carried along for plot.vif.ridge() are those of the chosen response
  vif(as.ridge(mod, response), ...)
}
