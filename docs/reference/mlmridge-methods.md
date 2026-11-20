# Methods for Multivariate Ridge Regression Objects

**Experimental.** Multivariate ridge regression (class `"mlmridge"`) is
new in this version; its functions, arguments and results may change in
future releases.

Basic methods for objects of class `"mlmridge"`, returned by
[`ridge()`](https://friendly.github.io/genridge/reference/ridge.md) when
the response is a matrix.

- [`coef()`](https://rdrr.io/r/stats/coef.html) extracts the
  coefficients: by default a list with one \\p \times q\\ matrix for
  each ridge constant; with `type = "array"`, a 3-d array indexed by
  `[lambda, predictor, response]`; with `response`, the \\\lambda \times
  p\\ matrix for one response, as for a `"ridge"` object.

- [`vcov()`](https://rdrr.io/r/stats/vcov.html) extracts the covariance
  matrices of the coefficients: by default a list with, for each ridge
  constant, the \\pq \times pq\\ covariance matrix of the stacked
  coefficients \\\mathrm{vec}(\mathbf{B})\\, which is
  `kronecker(Sigma, cov[[i]])`; with `response`, the list of \\p \times
  p\\ matrices for one response.

- [`print()`](https://rdrr.io/r/base/print.html) shows the coefficients
  for each response and the optimal ridge constants given by the HKB, LW
  and GCV criteria, for each response and pooled over responses.

- [`best()`](https://friendly.github.io/genridge/reference/ridge.md)
  returns the optimal ridge constants as a data frame, with one row for
  each response and criterion (HKB, LW, GCV), followed by a row for GCV
  pooled over responses (`response = "(pooled)"`; the pooling method is
  `object$pool`, see
  [`ridge()`](https://friendly.github.io/genridge/reference/ridge.md));
  with `response`, the univariate result for one response, as from
  [`best()`](https://friendly.github.io/genridge/reference/ridge.md) for
  a `"ridge"` object.

As for a `"ridge"` object, the HKB and LW values are closed-form
estimates, while the GCV values (per response and pooled) are the best
among the ridge constants supplied in `lambda` (or `df`).

## Usage

``` r
# S3 method for class 'mlmridge'
coef(object, response = NULL, type = c("list", "array"), ...)

# S3 method for class 'mlmridge'
vcov(object, response = NULL, ...)

# S3 method for class 'mlmridge'
best(object, response = NULL, ...)

# S3 method for class 'mlmridge'
print(x, digits = max(5, getOption("digits") - 5), ...)
```

## Arguments

- object, x:

  An object of class `"mlmridge"`

- response:

  If not `NULL`, a response (name or number) for which to return the
  results, in the same form as for a univariate `"ridge"` object; see
  [`as.ridge()`](https://friendly.github.io/genridge/reference/as.ridge.md).

- type:

  For [`coef()`](https://rdrr.io/r/stats/coef.html) with
  `response = NULL`: `"list"` (the default) or `"array"`.

- ...:

  Other arguments, not used.

- digits:

  Number of significant digits to print.

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) return a list, array or
matrix, as described above;
[`print()`](https://rdrr.io/r/base/print.html) returns `x`, invisibly.

       `best()` returns a data frame with columns `response` and `criterion` (factors), `k`
       (the optimal ridge constant) and, if the object contains the singular values of the
       predictors (`svd = TRUE`, the default), `df` (the equivalent degrees of freedom).

## See also

[`ridge()`](https://friendly.github.io/genridge/reference/ridge.md),
[`as.ridge()`](https://friendly.github.io/genridge/reference/as.ridge.md)

## Examples

``` r
if (requireNamespace("heplots", quietly = TRUE)) {
  data(schooldata, package = "heplots")
  sridge <- ridge(cbind(reading, mathematics, selfesteem) ~
                  education + occupation + visit + counseling + teacher,
                  data = schooldata, lambda = c(0, 0.1, 0.5, 1, 2, 5))
  sridge
  coef(sridge)[[4]]                       # coefficients for lambda = 1
  coef(sridge, response = "reading")      # lambda x p, for one response
  dim(coef(sridge, type = "array"))
  dim(vcov(sridge)[[1]])                  # pq x pq
  best(sridge)                            # optimal ridge constants, all responses
  best(sridge, response = "selfesteem")   # one response
}
#>              k       df
#> HKB 0.02679366 4.922545
#> LW  0.02268283 4.934046
#> GCV 0.00000000 5.000000
```
