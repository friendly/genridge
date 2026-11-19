# Extract a Univariate Ridge Regression from a Multivariate One

**Experimental.** Multivariate ridge regression (class `"mlmridge"`) is
new in this version; its functions, arguments and results may change in
future releases.

[`ridge()`](https://friendly.github.io/genridge/reference/ridge.md) with
a matrix response (e.g., `ridge(cbind(Y1, Y2) ~ x1 + x2, data)`) returns
an object of class `"mlmridge"`, fitting all responses with a common set
of ridge constants. `as.ridge()` extracts the results for one response
as an ordinary `"ridge"` object, so that all the methods for `"ridge"`
objects
([`traceplot()`](https://friendly.github.io/genridge/reference/traceplot.md),
[`plot.ridge()`](https://friendly.github.io/genridge/reference/plot.ridge.md),
[`pairs.ridge()`](https://friendly.github.io/genridge/reference/pairs.ridge.md),
[`plot3d.ridge()`](https://friendly.github.io/genridge/reference/plot3d.md),
[`precision()`](https://friendly.github.io/genridge/reference/precision.md),
[`vif.ridge()`](https://friendly.github.io/genridge/reference/vif.ridge.md),
[`pca()`](https://friendly.github.io/genridge/reference/pca.md), ...)
can be applied to it.

The result is the same as fitting
[`ridge()`](https://friendly.github.io/genridge/reference/ridge.md) to
that response alone, with the same `lambda` or `df`.

## Usage

``` r
as.ridge(x, ...)

# S3 method for class 'ridge'
as.ridge(x, ...)

# S3 method for class 'mlmridge'
as.ridge(x, response = 1, ...)
```

## Arguments

- x:

  An object of class `"mlmridge"`

- ...:

  Not used

- response:

  The response to extract: a name, or a number in `1:q`. The default is
  the first response.

## Value

An object of class `"ridge"`; see
[`ridge()`](https://friendly.github.io/genridge/reference/ridge.md) for
its components.

## See also

[`ridge()`](https://friendly.github.io/genridge/reference/ridge.md)

## Examples

``` r
if (requireNamespace("heplots", quietly = TRUE)) {
  data(schooldata, package = "heplots")
  sridge <- ridge(cbind(reading, mathematics, selfesteem) ~
                  education + occupation + visit + counseling + teacher,
                  data = schooldata, lambda = c(0, 1, 2, 5, 10, 20))
  # results for one response
  read.ridge <- as.ridge(sridge, "reading")
  read.ridge
  traceplot(read.ridge)
}
```
