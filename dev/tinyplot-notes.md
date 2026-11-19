# Should genridge's base-graphics plot methods use tinyplot?

*Notes started 2026-11-19, during the 0.8.1 → mlmridge development cycle. Pilot code:
`dev/tinyplot-pilot.R`.*

## What tinyplot is (checked 2026-11-19)

- **tinyplot 0.8.0** (CRAN, 2026-09-28), by Grant McDermott, Vincent Arel-Bundock and
  others.
- It's an extension of **base graphics**, not a new graphics system. `tinyplot()` and `plt()`
  draw with base R, so `lines()`, `text()`, `abline()`, `par()` and friends still work on top.
- **Dependencies:** none outside base R. Imports are graphics, grDevices, stats, tools and
  utils; it needs R ≥ 4.0. So it's a very light dependency.
- **Features relevant to genridge:**
  - a formula interface: `y ~ x | by`;
  - grouping (`by =`) with automatic legends and palettes;
  - **facets** (`facet = ~ f`, with free scales through `facet.args`);
  - **themes** (`tinytheme("clean")`, `"classic"`, `"minimal"`, … 22 built in);
  - `tinyplot_add()` to add layers;
  - transparency (`alpha`, `fill`);
  - 37 plot types: `type_polygon`, `type_segments`, `type_text` (with a `repel` option),
    `type_vline`, `type_ribbon` and others.
- **Not directly usable:**
  - `type_ellipse()` draws *data* (confidence) ellipses computed from grouped points. genridge
    draws ellipses from a **given** centre and covariance matrix, so it would still compute
    the ellipse coordinates itself and draw them with `type = "polygon"`.
  - `tinypairs()` is a pairs plot of the columns of a data frame. `pairs.ridge()` draws custom
    panels (ellipse traces), so tinypairs can't replace it.

## Pilot (`dev/tinyplot-pilot.R`, schooldata, 3 responses)

1. **`traceplot` for all responses in one call.** I reshaped `coef(sridge, type = "array")`
   to long form, then called
   `tinyplot(coef ~ lambda | predictor, facet = ~ response, type = "b", facet.args = list(free = TRUE))`.
   - This gives three panels with free y-scales, one legend for the predictors, and a
     consistent theme. It's simpler than the `par(mfrow)` loop in `traceplot.mlmridge()`, and
     the layout doesn't need restoring afterwards.
   - It has a *legend* instead of genridge's **direct labels** at the end of each trace. Direct
     labels are a deliberate genridge design choice, and they avoid the clipping problem in
     `traceplot.ridge()`. `type_text(repel = TRUE)` might reproduce them.
   - HKB/LW reference lines differ per facet. They'd need `tinyplot_add(type = "vline")` with
     per-facet data; not tried yet.
2. **`plot.ridge`-style bivariate trace for one response.** I computed the ellipse coordinates
   (radius 0.5) as genridge does, then called
   `tinyplot(y ~ x | lambda, type = "polygon", fill = 0.15)`, and added the centres with
   `tinyplot_add(type = "p")` and the path with base `lines()`.
   - It looks good: translucent filled ellipses and a λ legend.
   - λ is *ordered*, so a sequential palette would suit it better than the default
     qualitative one.
   - genridge's labels at the ellipse centres would be `type = "text"`.

## Options

| Option | What | For | Against |
|---|---|---|---|
| A. Replace | Rewrite all base plot methods with tinyplot | One consistent look; facets and themes everywhere | Large rewrite; changes the look of existing plots and possibly behaviour that users' code relies on (e.g., adding to a plot with `lines()`, `par()` settings); new Imports dependency |
| B. Engine switch | Keep base graphics as the default and add `engine = c("base", "tinyplot")` | No breakage; people can opt in | Two code paths to maintain and test |
| C. New plots only | Use tinyplot for the *new* multivariate plots (faceted `traceplot.mlmridge`, cross-response ellipses §5.3, maybe ridge HE plots) and leave the univariate methods as base | Biggest gain where facets matter; limited scope | Two looks in one package |
| D. No change | Stay with base | No work | Misses facets and themes |

**Leaning: C, then perhaps B.**
- Facets are the main gain, and they matter most for the multivariate plots: one panel per
  response.
- The univariate methods work, and they are what the JCGS paper and the book show.
- tinyplot is base graphics underneath, so a tinyplot `traceplot.mlmridge()` still combines
  with base annotations.

**Decision (MF, 2026-11-19): option C is the direction**, using tinyplot for the new
multivariate plots, but **not yet**. Finish the multivariate work in base graphics first, as
now. Take up the tinyplot versions as a separate step, after reviewing the pilot images.

## Questions to settle

1. **Dependency:** should tinyplot be in Imports, which is cheap since it only uses base R, or
   in Suggests with a fallback to base? Imports is simpler if any default plot uses it.
2. **Direct labels vs legends:** keep genridge's end-of-trace labels (`type_text`, maybe with
   `repel`), or accept legends?
3. **Theme:** use the tinyplot default, or set a genridge theme (`tinytheme_register()`)?
   Don't change the user's global theme; set and restore it inside each method.
4. **Palettes:** sequential for λ (ellipse plots), qualitative for predictors (trace plots).
   genridge currently uses fixed colour vectors: `col = c("black", "red", "darkgreen", ...)`.
5. **Return values:** base methods return `NULL` or invisible data. tinyplot can record the
   plot (`record = TRUE`). Is that useful?
6. **pkgdown and vignette:** tinyplot output works in knitr like base graphics. Check figure
   sizing with facets.
7. **Maturity:** tinyplot is pre-1.0 (0.8.0), so its API may still change. Pin a minimum
   version and watch its NEWS.

## Next steps (when this is taken up)

- Decide between options C and B.
- Prototype `traceplot.mlmridge()` with tinyplot: facets, direct labels with `type_text`,
  and HKB/LW lines per facet. Compare it side by side with the current base version in the
  gallery (`dev/mlmridge-plot-examples.R`).
- If the cross-response ellipses (§5.3 of `ridge_mlm-plan.md`) are built next, consider
  writing them with tinyplot from the start.
