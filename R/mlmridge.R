#' Extract a Univariate Ridge Regression from a Multivariate One
#'
#' @description
#' **Experimental.** Multivariate ridge regression (class `"mlmridge"`) is new in this
#' version; its functions, arguments and results may change in future releases.
#'
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
#' **Experimental.** Multivariate ridge regression (class `"mlmridge"`) is new in this
#' version; its functions, arguments and results may change in future releases.
#'
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
#' **Experimental.** Multivariate ridge regression (class `"mlmridge"`) is new in this
#' version; its functions, arguments and results may change in future releases.
#'
#' These methods apply the plots and summaries for univariate ridge regression to one response
#' of a multivariate ridge regression (class `"mlmridge"`). Each extracts the results for that
#' response with [as.ridge()] and calls the corresponding method for `"ridge"` objects:
#'
#' * `traceplot()`: univariate ridge trace plots ([traceplot()]), by default for all responses,
#'   one panel each
#' * `plot()`: bivariate ridge trace plot of coefficients with covariance ellipses ([plot.ridge()]);
#'   with `predictor =`, coefficient ellipses for one or more predictors in the space of two
#'   responses (see **Details**)
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
#' @details
#' **Two kinds of plot from `plot()`.** `plot()` for an `"mlmridge"` object draws one of two
#' different displays, depending on whether `predictor` is given:
#'
#' * **Within one response** (the default, `response =`): as for a univariate ridge regression
#'   ([plot.ridge()]), the axes are two *predictors*, and each ellipse shows the joint
#'   variability of their coefficients for that response, at one value of \eqn{\lambda}.
#' * **Across responses** (`predictor =`): the axes are two *responses* (`responses =`), and each
#'   ellipse shows the joint variability of one predictor's coefficients for those two responses,
#'   at one value of \eqn{\lambda}. The path of the centres traces how the predictor's effects on
#'   both responses change together as \eqn{\lambda} increases. With several predictors, each has
#'   its own path, labelled with its name.
#'
#' For predictor \eqn{j}, the coefficients for all responses (row \eqn{j} of \eqn{\mathbf{B}(\lambda)})
#' have covariance matrix \eqn{[\mathbf{V}_\lambda]_{jj} \, \boldsymbol{\Sigma}}, where
#' \eqn{\boldsymbol{\Sigma}} is the residual covariance matrix of the responses (`x$Sigma`) and
#' \eqn{\mathbf{V}_\lambda} depends only on the predictors (`x$cov`). So in the across-responses
#' plot, all the ellipses have the **same shape and orientation**, that of the residual
#' covariance of the two responses: only their centres and sizes change with \eqn{\lambda} and
#' between predictors.
#'
#' The across-responses display is the ridge analogue of [heplots::coefplot()] for a
#' multivariate linear model, which draws confidence ellipses for the coefficients of each
#' predictor in the space of two responses. There are two differences in meaning:
#'
#' * `coefplot()` shows each predictor once, for the OLS fit; here each predictor is traced over
#'   the values of \eqn{\lambda}.
#' * `coefplot()` scales its ellipses as confidence regions (by `level`, using an \eqn{F} or
#'   \eqn{\chi^2} quantile), so whether an ellipse covers the origin is a test of the
#'   coefficients. Here, as in [plot.ridge()], the ellipses are scaled by `radius` (1 = standard
#'   "unit" ellipse) and describe the sampling variability around a *biased* estimate; they are
#'   not confidence regions for the true coefficients, and covering the origin is not a test.
#'
#' @name mlmridge-plots
#' @param x,object,mod An object of class `"mlmridge"`
#' @param response The response to use: a name, or a number in `1:q`. If not given, the first
#'        response is used, with a message. For `traceplot()`, `response` can name several
#'        responses, and the default (`NULL`) plots all of them, one panel per response.
#' @param main Title for the plot. The default is the name of the response. For `traceplot()`
#'        with several responses, recycled over the panels.
#' @param layout For `traceplot()` with several responses: if `TRUE` (the default), set up a
#'        multi-panel layout with one panel per response, side by side where possible, and
#'        restore the previous layout afterwards. Use `layout = FALSE` to draw into a layout you
#'        have set up yourself, e.g., with `par(mfrow = )`.
#' @param predictor For `plot()`: if given, plot coefficient ellipses across responses for this
#'        predictor, or several predictors (names or numbers), instead of the within-response
#'        plot; see **Details**. Further arguments to control the plot (`radius`, `which.lambda`,
#'        `labels`, `col`, `fill`, `fill.alpha`, `lwd`, `lty`, `center.pch`, `center.cex`, `ref`,
#'        `ref.col`, `xlim`, `ylim`, `xlab`, `ylab`) have the same meaning as in [plot.ridge()].
#'        By default, the centres are labelled with the values of \eqn{\lambda} for one predictor,
#'        and with the predictor names for several.
#' @param responses For `plot()` with `predictor`: the two responses (names or numbers) giving
#'        the horizontal and vertical axes. Default: the first two.
#' @param det.fun,normalize See [precision()].
#' @param ... Other arguments, passed to the method for `"ridge"` objects.
#' @return As for the corresponding method for `"ridge"` objects.
#' @seealso [ridge()], [as.ridge()], [mlmridge-methods]
#' @importFrom graphics pairs lines polygon points text abline
#' @importFrom stats biplot
#' @examples
#' if (requireNamespace("heplots", quietly = TRUE)) {
#'   data(schooldata, package = "heplots")
#'   lambda <- c(0, 0.1, 0.5, 1, 2, 5)
#'   sridge <- ridge(cbind(reading, mathematics, selfesteem) ~
#'                   education + occupation + visit + counseling + teacher,
#'                   data = schooldata, lambda = lambda)
#'
#'   # all responses, one panel each
#'   traceplot(sridge)
#'   traceplot(sridge, X = "df")
#'   # one response
#'   traceplot(sridge, response = "selfesteem")
#'
#'   # bivariate ridge trace plot for two predictors
#'   plot(sridge, response = "mathematics", variables = c("occupation", "visit"),
#'        radius = 0.5)
#'
#'   # coefficient ellipses across two responses, for one predictor ...
#'   plot(sridge, predictor = "visit", responses = c("reading", "mathematics"), fill = TRUE)
#'   # ... and for all predictors, at lambda = 0 and 5
#'   plot(sridge, predictor = 1:5, which.lambda = c(1, 6), col = c("black", "red"))
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
traceplot.mlmridge <- function(x, response = NULL, main, layout = TRUE, ...) {
  ynames <- colnames(x$Sigma)
  if (is.null(response)) response <- ynames
  if (is.numeric(response)) response <- ynames[response]
  bad <- setdiff(response, ynames)
  if (length(bad) || anyNA(response))
    stop("response must be in: ", paste(ynames, collapse = ", "), ", or numbers in 1:", length(ynames))
  nr <- length(response)
  if (missing(main)) main <- response
  main <- rep_len(main, nr)
  if (layout && nr > 1) {
    op <- par(mfrow = rev(grDevices::n2mfrow(nr)))   # side by side: 1 x 3, 2 x 2, 2 x 3, ...
    on.exit(par(op))
  }
  for (i in seq_len(nr))
    traceplot(as.ridge(x, response[i]), main = main[i], ...)
  invisible(NULL)
}

#' @rdname mlmridge-plots
#' @exportS3Method plot mlmridge
plot.mlmridge <- function(x, response = 1, predictor = NULL, responses = 1:2, main, ...) {
  if (!is.null(predictor)) {
    if (!missing(response))
      warning("'response' is ignored when 'predictor' is given; use 'responses' to choose the two responses")
    if (missing(main)) main <- if (length(predictor) == 1) predictor else ""
    return(invisible(plot_cross_response(x, predictor = predictor, responses = responses,
                                         main = main, ...)))
  }
  sel <- select_response(x, response, missing(response))
  if (missing(main)) main <- sel$name
  plot(sel$ridge, main = main, ...)
}

# Cross-response coefficient ellipses: for each predictor, its coefficients for two responses
# across lambda, with covariance ellipses [V_lambda]_jj * Sigma[r, r]. Arguments follow plot.ridge().
plot_cross_response <- function(x, predictor, responses = 1:2,
                                radius = 1, which.lambda = seq_along(x$lambda),
                                labels = NULL, pos = 3, cex = 1.2,
                                lwd = 2, lty = 1, xlim, ylim,
                                col = c("black", "red", "darkgreen", "blue", "darkcyan", "magenta",
                                        "brown", "darkgray"),
                                center.pch = 16, center.cex = 1.5,
                                fill = FALSE, fill.alpha = 0.3, ref = TRUE, ref.col = gray(.70),
                                xlab, ylab, main = "", ...) {
  ynames <- colnames(x$Sigma)
  xnames <- rownames(x$coef[[1]])
  if (is.numeric(responses)) responses <- ynames[responses]
  if (length(responses) != 2 || !all(responses %in% ynames))
    stop("'responses' must be two of: ", paste(ynames, collapse = ", "))
  if (is.numeric(predictor)) predictor <- xnames[predictor]
  if (!all(predictor %in% xnames))
    stop("'predictor' must be among: ", paste(xnames, collapse = ", "))

  lambda <- x$lambda[which.lambda]
  nl <- length(lambda)
  col <- rep_len(col, nl)
  lwd <- rep_len(lwd, nl)
  lty <- rep_len(lty, nl)
  fill <- rep_len(fill, nl)
  fill.alpha <- rep_len(fill.alpha, nl)
  fill.col <- ifelse(fill, mapply(function(cl, a) grDevices::adjustcolor(cl, alpha.f = a),
                                  col, fill.alpha), NA)
  # labels at the ellipse centres: lambda values for one predictor, predictor names for several
  if (is.null(labels)) labels <- if (length(predictor) == 1) lambda else predictor

  S <- x$Sigma[responses, responses]
  angles <- seq(0, 2 * pi, length.out = 61)
  circle <- radius * cbind(cos(angles), sin(angles))
  centers <- lapply(predictor, function(j)
    t(sapply(which.lambda, function(i) x$coef[[i]][j, responses])))
  ells <- lapply(predictor, function(j) lapply(which.lambda, function(i) {
    sweep(circle %*% chol(x$cov[[i]][j, j] * S), 2, x$coef[[i]][j, responses], "+")
  }))

  allxy <- do.call(rbind, unlist(ells, recursive = FALSE))
  if (missing(xlim)) xlim <- range(allxy[, 1], if (ref) 0)
  if (missing(ylim)) ylim <- range(allxy[, 2], if (ref) 0)
  if (missing(xlab)) xlab <- paste("Coefficient for", responses[1])
  if (missing(ylab)) ylab <- paste("Coefficient for", responses[2])

  plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main, ...)
  if (ref) abline(h = 0, v = 0, col = ref.col)
  for (p in seq_along(predictor)) {
    for (i in seq_len(nl))
      polygon(ells[[p]][[i]], col = fill.col[i], border = col[i], lwd = lwd[i], lty = lty[i])
    ctr <- centers[[p]]
    lines(ctr)
    points(ctr, pch = center.pch, cex = center.cex, col = col)
    if (!is.null(labels)) {
      if (length(predictor) == 1) text(ctr, labels = labels, pos = pos, cex = cex)
      else text(ctr[1, 1], ctr[1, 2], labels = labels[p], pos = pos, cex = cex)
    }
  }
  invisible(list(centers = stats::setNames(centers, predictor),
                 ellipses = stats::setNames(ells, predictor)))
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
