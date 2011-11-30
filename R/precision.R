## Measures of precision and bias for ridge regression

precision <- function(object, ...) {
	UseMethod("precision")
}

# TODO:  should we use log.det or det()^{1/p} ??
precision.ridge <- function(object, ...) {
	tr <- function(x) sum(diag(x))
	maxeig <- function(x) max(eigen(x)$values)

	V <- object$cov
	p <- ncol(coef(object))
	ldet <- log(unlist(lapply(V, det)))
	trace <- unlist(lapply(V, tr))
	meig <- unlist(lapply(V, maxeig))	
	norm <- sqrt(rowMeans(coef(object)^2))
	data.frame(lambda=object$lambda, df=object$df, log.det=ldet, trace=trace, max.eig=meig, norm.beta=norm)
}

precision.lm <- function(object, ...) {
	V <- vcov(object)
	beta <- coefficients(object)
	if (names(beta[1]) == "(Intercept)") {
		V <- V[-1, -1]
		beta <- beta[-1]
	}
#	else warning("No intercept: precision may not be sensible.")
	ldet <- log(det(V))
	trace <- sum(diag(V))
	meig <- max(eigen(V)$values)
	norm <- sqrt(mean(beta^2))
	res <- list(df=length(beta), log.det=ldet, trace=trace, max.eig=meig, norm.beta=norm)
	unlist(res)
}
