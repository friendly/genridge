# Enhanced Contour Plots

This is an enhancement to
[`contour`](https://rdrr.io/r/graphics/contour.html), written as a
wrapper for that function. It creates a contour plot, or adds contour
lines to an existing plot, allowing the contours to be filled and
returning the list of contour lines.

## Usage

``` r
contourf(
  x = seq(0, 1, length.out = nrow(z)),
  y = seq(0, 1, length.out = ncol(z)),
  z,
  nlevels = 10,
  levels = pretty(zlim, nlevels),
  zlim = range(z, finite = TRUE),
  col = par("fg"),
  color.palette = colorRampPalette(c("white", col)),
  fill.col = color.palette(nlevels + 1),
  fill.alpha = 0.5,
  add = FALSE,
  ...
)
```

## Arguments

- x, y:

  locations of grid lines at which the values in `z` are measured. These
  must be in ascending order. By default, equally spaced values from 0
  to 1 are used. If `x` is a list, its components `x$x` and `x$y` are
  used for x and y, respectively. If the list has component `x$z` this
  is used for `z`.

- z:

  a matrix containing the values to be plotted (NAs are allowed). Note
  that `x` can be used instead of `z` for convenience.

- nlevels:

  number of contour levels desired **iff** levels is not supplied

- levels:

  numeric vector of levels at which to draw contour lines

- zlim:

  z-limits for the plot. x-limits and y-limits can be passed through ...

- col:

  color for the lines drawn

- color.palette:

  a color palette function to be used to assign fill colors in the plot

- fill.col:

  a call to the `color.palette` function or an an explicit set of colors
  to be used in the plot. Use `fill.col=NULL` to suppress the filled
  polygons. a vector of fill colors corresponding to levels. By default,
  a set of possibly transparent colors is calculated ranging from white
  to `col`, using transparency given by `fill.alpha`

- fill.alpha:

  transparency value for `fill.col`, either a hex character string, or a
  numeric value between 0 and 1. Use `fill.alpha=NA` to suppress
  transparency.

- add:

  logical. If `TRUE`, add to a current plot.

- ...:

  additional arguments passed to
  [`contour`](https://rdrr.io/r/graphics/contour.html), including all
  arguments of
  [`contour.default`](https://rdrr.io/r/graphics/contour.html) not
  mentioned above, as well as additional graphical parameters passed by
  [`contour.default`](https://rdrr.io/r/graphics/contour.html) to more
  basic functions.

## Value

Returns invisibly the list of contours lines, with components `levels`,
`x`, `y`. See
[`contourLines`](https://rdrr.io/r/grDevices/contourLines.html).

## See also

[`contour`](https://rdrr.io/r/graphics/contour.html),
[`contourLines`](https://rdrr.io/r/grDevices/contourLines.html)

[`contourplot`](https://rdrr.io/pkg/lattice/man/levelplot.html) from
package lattice.

## Author

Michael Friendly

## Examples

``` r

x <- 10*1:nrow(volcano)
y <- 10*1:ncol(volcano)
contourf(x,y,volcano, col="blue")

contourf(x,y,volcano, col="blue", nlevels=6)


# return value, unfilled, other graphic parameters
res <- contourf(x,y,volcano, col="blue", fill.col=NULL, lwd=2)

# levels used in the plot
sapply(res, function(x) x[[1]])
#>  [1] 100 100 100 110 110 110 110 120 130 140 150 150 160 160 170 170 170 180 180
#> [20] 190

```
