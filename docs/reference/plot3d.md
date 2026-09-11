# 3D Ridge Trace Plots

The 3D ridge trace plot displays 3D projections of the covariance
ellipsoids for a set of ridge regression estimates indexed by a ridge
tuning constant.

The centers of these ellipses show the bias induced for each parameter,
and also how the change in the ridge estimate for one parameter is
related to changes for other parameters.

The size and shapes of the covariance ellipsoids show directly the
effect on precision of the estimates as a function of the ridge tuning
constant.

`plot3d.ridge` and `plot3d.pcaridge` differ only in the defaults for the
variables plotted.

## Usage

``` r
plot3d(x, ...)

# S3 method for class 'pcaridge'
plot3d(x, variables = (p - 2):p, ...)

# S3 method for class 'ridge'
plot3d(
  x,
  variables = 1:3,
  radius = 1,
  which.lambda = 1:length(x$lambda),
  lwd = 1,
  lty = 1,
  xlim,
  ylim,
  zlim,
  xlab,
  ylab,
  zlab,
  col = c("black", "red", "darkgreen", "blue", "darkcyan", "magenta", "brown",
    "darkgray"),
  labels = lambda,
  ref = TRUE,
  ref.col = gray(0.7),
  segments = 40,
  shade = TRUE,
  shade.alpha = 0.1,
  wire = FALSE,
  aspect = 1,
  add = FALSE,
  ...
)
```

## Arguments

- x:

  A `ridge` object, as fit by
  [`ridge`](https://friendly.github.io/genridge/reference/ridge.md) or a
  `pcaridge` object as transformed by
  [`pca.ridge`](https://friendly.github.io/genridge/reference/pca.md)

- ...:

  Other arguments passed down

- variables:

  Predictors in the model to be displayed in the plot: an integer or
  character vector of length 3, giving the indices or names of the
  variables. Defaults to the first three predictors for `ridge` objects
  or the *last* three dimensions for `pcaridge` objects.

- radius:

  Radius of the ellipse-generating circle for the covariance ellipsoids.
  The default, `radius=1` gives a standard “unit” ellipsoid. Typically,
  `radius<1` gives less cluttered displays.

- which.lambda:

  A vector of indices used to select the values of `lambda` for which
  ellipsoids are plotted. The default is to plot ellipsoids for all
  values of `lambda` in the `ridge` object.

- lwd, lty:

  Line width and line type for the covariance ellipsoids. Recycled as
  necessary.

- xlim, ylim, zlim:

  X, Y, Z limits for the plot, each a vector of length 2. If missing,
  the range of the covariance ellipsoids is used.

- xlab, ylab, zlab:

  Labels for the X, Y, Z variables in the plot. If missing, the names of
  the predictors given in `variables` is used.

- col:

  A numeric or character vector giving the colors used to plot the
  covariance ellipsoids. Recycled as necessary.

- labels:

  A numeric or character vector giving the labels to be drawn at the
  centers of the covariance ellipsoids.

- ref:

  Logical: whether to draw horizontal and vertical reference lines at 0.
  This is not yet implemented.

- ref.col:

  Color of reference lines.

- segments:

  Number of line segments used in drawing each dimension of a covariance
  ellipsoid.

- shade:

  a logical scalar or vector, indicating whether the ellipsoids should
  be rendered with
  [`shade3d`](https://dmurdoch.github.io/rgl/dev/reference/shade3d.html).
  Recycled as necessary.

- shade.alpha:

  a numeric value in the range \[0,1\], or a vector of such values,
  giving the alpha transparency for ellipsoids rendered with
  `shade=TRUE`.

- wire:

  a logical scalar or vector, indicating whether the ellipsoids should
  be rendered with
  [`wire3d`](https://dmurdoch.github.io/rgl/dev/reference/shade3d.html).
  Recycled as necessary.

- aspect:

  a scalar or vector of length 3, or the character string "iso",
  indicating the ratios of the x, y, and z axes of the bounding box. The
  default, `aspect=1` makes the bounding box display as a cube
  approximately filling the display. See
  [`aspect3d`](https://dmurdoch.github.io/rgl/dev/reference/aspect3d.html)
  for details.

- add:

  if `TRUE`, add to the current `rgl` plot; the default is `FALSE`.

## Value

None. Used for its side-effect of plotting

## Note

This is an initial implementation. The details and arguments are subject
to change.

## References

Friendly, M. (2013). The Generalized Ridge Trace Plot: Visualizing Bias
*and* Precision. *Journal of Computational and Graphical Statistics*,
**22**(1), 50-68, doi:10.1080/10618600.2012.681237,
<https://www.datavis.ca/papers/genridge-jcgs.pdf>

## See also

[`plot.ridge`](https://friendly.github.io/genridge/reference/plot.ridge.md),
[`pairs.ridge`](https://friendly.github.io/genridge/reference/pairs.ridge.md),
[`pca.ridge`](https://friendly.github.io/genridge/reference/pca.md)

## Author

Michael Friendly

## Examples

``` r

lmod <- lm(Employed ~ GNP + Unemployed + Armed.Forces + Population + 
                      Year + GNP.deflator, data=longley)
longley.y <- longley[, "Employed"]
longley.X <- model.matrix(lmod)[,-1]

lambda <- c(0, 0.005, 0.01, 0.02, 0.04, 0.08)
lambdaf <- c("0", ".005", ".01", ".02", ".04", ".08")
lridge <- ridge(longley.y, longley.X, lambda=lambda)

# \donttest{
plot3d(lridge, var=c(1,4,5), radius=0.5)

# view in SVD/PCA space
plridge <- pca(lridge)
plot3d(plridge, radius=0.5)
# }


```
