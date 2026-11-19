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

## Usage

``` r
# S3 method for class 'mlmridge'
coef(object, response = NULL, type = c("list", "array"), ...)

# S3 method for class 'mlmridge'
vcov(object, response = NULL, ...)

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
}
#> [1] 15 15
```
