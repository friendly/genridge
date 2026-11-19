# Multivariate ridge regression in genridge: planning notes

*Started 2026-11-19. Source draft: `dev/ridge_mlm.R` (adapted from `VARshrink::lm_multiv_ridge()`).*

The goal is to extend `ridge()` from a single response `y` to a matrix of responses `Y`
(a multivariate linear model, MLM), in a way that fits genridge's existing design: the `ridge`
object, its methods (`coef`, `vcov`, `print`, `best`) and, above all, its graphical methods
(`traceplot`, `plot`, `pairs`, `plot3d`, `pca`/`biplot`, `precision`, `vif`). This also links
genridge to the other MLM packages (heplots, candisc, mvinfluence) and to the Vis-MLM book,
chapter 9 (`09-collinearity-ridge.qmd`).

---

## 1. What the draft does

`ridge_mlm(Y, X, lambda, do_scale)` computes, via the SVD of `X`,

$$\widehat{\mathbf{B}}(\lambda) = (\mathbf{X}^\top \mathbf{X} + n\lambda \mathbf{I})^{-1} \mathbf{X}^\top \mathbf{Y}$$

for a vector of `lambda`, and returns `list(B, lambda, GCV)`: a list of p × q coefficient
matrices and a GCV score per `lambda`, pooled over all responses.

What it lacks compared with genridge's `ridge.default()`:

- no covariance matrices of the estimates (which all of genridge's ellipse plots need);
- no `df` parameterization;
- no HKB / LW criteria;
- no intercept handling, unless `do_scale = TRUE`; `do_scale = FALSE` is the default, and then
  `X` and `Y` aren't even centered;
- no formula or model interface, no class or methods.

**Licence.** VARshrink is GPL (>= 3); genridge is GPL (>= 2). The arithmetic is short and
standard, so I recommend **reimplementing it in genridge's own style** (reusing the
`ridge.default()` code), not copying the code. Cite VARshrink and Golub, Heath & Wahba (1979)
for the GCV formula in the docs.

## 2. Facts that shape the design (verified numerically, 2026-11-19)

Checked on `longley` with two responses (scratch script, not kept):

1. **Columns are just univariate ridges.** With the same λ for every response, column *k* of
   $\widehat{\mathbf{B}}(\lambda)$ is exactly `ridge(Y[, k], X, lambda)`: the max |difference|
   was ~1e-12. The penalty involves only `X`, so a shared-λ multivariate ridge doesn't "borrow
   strength" across responses in the estimates themselves.

2. **λ scales differ.** genridge's λ equals n × VARshrink's λ when `X` is scaled the same way.
   - genridge scales `X` to unit mean square (divisor n), so $\mathbf{X}^\top\mathbf{X}$ has
     diagonal n.
   - VARshrink's `do_scale` uses `scale()` (divisor n − 1) and puts nλ in the penalty.
   - So genridge's λ must stay as it is, and any examples taken from VARshrink need converting.

3. **The covariance has Kronecker structure.** With $\mathbf{A}_\lambda = (\mathbf{X}^\top
   \mathbf{X} + \lambda \mathbf{I})^{-1}\mathbf{X}^\top$,

   $$\mathrm{Var}(\mathrm{vec}\,\widehat{\mathbf{B}}) = \boldsymbol{\Sigma} \otimes \mathbf{V}_\lambda,
   \qquad \mathbf{V}_\lambda = \mathbf{A}_\lambda \mathbf{A}_\lambda^\top
   = (\mathbf{X}^\top \mathbf{X} + \lambda \mathbf{I})^{-1}\mathbf{X}^\top\mathbf{X}
     (\mathbf{X}^\top \mathbf{X} + \lambda \mathbf{I})^{-1}.$$

   - Σ is the q × q residual covariance matrix; estimate it from the OLS fit as
     **E** / (n − p − 1).
   - genridge's univariate `cov[[i]]` is already `MSE * V_λ`. I checked that `cov / MSE` is the
     same matrix for two different responses, and that `MSE` equals the diagonal element
     `E[k, k]`.
   - So the univariate machinery carries over unchanged: store **Σ** once, plus the list of
     `V_λ`.

4. **Consequences:**
   - `df(λ)`, the shrinkage factors, the SVD, `pca()` rotations and VIFs depend **only on X**.
     They are the same for every response, and identical to the univariate `ridge` object.
   - `precision()`: $\log|\mathbf{\Sigma} \otimes \mathbf{V}_\lambda| = p \log|\boldsymbol{\Sigma}|
     + q \log|\mathbf{V}_\lambda|$. It is the univariate measure shifted and rescaled, so it
     orders the λ values the same way. Trace and max-eigenvalue factor similarly.
   - **What's genuinely new in the multivariate case:**
     - the *cross-response* covariance of the coefficients for a given predictor,
       $\mathrm{Var}(\widehat{\mathbf{b}}_{j\cdot}) = [\mathbf{V}_\lambda]_{jj}\,\boldsymbol{\Sigma}$;
     - a pooled choice of λ, from a GCV summed over responses;
     - the behaviour of the hypothesis (**H**) and error (**E**) SSP matrices as λ varies, i.e.,
       the HE-plot view.

## 3. Interface (recommendation)

`ridge()` now has `.formula`, `.default` and (since 2026-11-19) `.lm` methods. The `.lm` method
currently *stops* for `mlm` objects. The natural extension:

| Call | Dispatch |
|---|---|
| `ridge(cbind(Y1, Y2) ~ x1 + x2, data, lambda)` | `ridge.formula()`: detect a matrix response and return the multivariate result |
| `ridge(lm(cbind(Y1, Y2) ~ ...), lambda)` | new `ridge.mlm()` (class `mlm` is dispatched before `lm`) |
| `ridge(Ymatrix, X, lambda)` | `ridge.default()`: branch on `is.matrix(y) && ncol(y) > 1` |

All three call one internal workhorse. The simplest version is to generalize `ridge.default()`
to accept a matrix `y` and have the univariate case be q = 1. Alternatively, keep it as is and
add `ridge_mlm_fit()`; that is lower risk for existing users.

**λ:** one shared λ (or `df`) for all responses, as in the draft. Response-specific λ would
just be q separate univariate ridges, and users can already do that.

## 4. The result object

Proposed class: **`c("mlmridge", "ridge")`**? Or a separate `"mlmridge"` only? Decide after
looking at which existing methods would silently do the wrong thing on a 3-d `coef`.

Likely components, mirroring `ridge`:

- `lambda`, `df`: as now.
- `coef`: a 3-d array `[lambda, predictor, response]`. A list of p × q `B` matrices, as in
  VARshrink, is the alternative; the array makes `coef[, , "Y1"]` give the univariate matrix
  directly.
- `cov`: a list of p × p `V_λ` (the shared X part, without Σ).
- `Sigma`: the q × q residual covariance from the OLS fit.
- `mse`: a matrix `[lambda, response]`.
- `scales`, `svd.D`, `svd.U`, `svd.V`: as now; they depend only on X.
- Criteria: HKB and LW per response (a matrix); GCV both pooled (as in VARshrink) and per
  response; `kGCV` from the pooled GCV.

**Key helper: extract a univariate `ridge` object for one response**, e.g.
`as.ridge(x, response = "Y1")` or `x[["Y1"]]`. It sets `coef = coef[, , k]`,
`cov = Sigma[k, k] * V_λ` and so on, so **every existing method works per response for free**:
`traceplot()`, `plot()`, `pairs()`, `plot3d()`, `precision()`, `vif()`, `pca()`, `biplot()`.
This should be step 1 of the implementation.

`coef()`, `vcov()` and `print()` methods:

- `coef(x, response = NULL)`: the array, or one response's matrix.
- `vcov(x, lambda = )`: the full pq × pq `kronecker(Sigma, V_λ)`, or per response.

## 5. Plot methods: ideas, roughly in order of value

1. **Per-response delegation (cheap, do first).** `plot(x, response = "Y1")` etc. call the
   univariate methods through `as.ridge()`. A `response = NULL` default could loop over the
   responses, or make a multi-panel `traceplot`.

2. **`traceplot.mlmridge()`.** Coefficient paths vs λ or df:
   - either one panel per response, which is the easiest to read, sharing the λ axis and the
     HKB/LW/GCV reference lines;
   - or one panel with colour = predictor and line type = response.

3. **Cross-response coefficient ellipses (new, and genuinely multivariate).** For a single
   predictor *j*, plot the q-vector $\widehat{\mathbf{b}}_{j\cdot}(\lambda)$ for two responses,
   with covariance ellipses $[\mathbf{V}_\lambda]_{jj}\,\boldsymbol{\Sigma}$ across λ.
   - Because of the Kronecker structure, the ellipses all have the **same shape and
     orientation** (that of Σ). Only their size and centre change with λ. That's a nice
     teaching point.
   - A `pairs()` version would take all response pairs for one predictor.

4. **Ridge HE plots (new; links to heplots).**
   - **H**(λ) = $\widehat{\mathbf{B}}(\lambda)^\top \mathbf{X}^\top\mathbf{X}\,\widehat{\mathbf{B}}(\lambda)$
     (the fitted SSP);
   - **E**(λ) = the residual SSP.
   - Show how the **H** ellipse shrinks and **E** grows as λ increases, in the response space.
   - Open questions:
     - Is the overall regression **H** the right object, or one per predictor (or term), as in
       `heplots::heplot()` for an `mlm`?
     - What does "significance" scaling mean for a biased estimator? It's probably best to
       show size only, without the test scaling.
     - Could this reuse `heplots::heplot()` by building a pseudo-`mlm`, or should it be
       implemented directly?

5. **PCA / biplot.**
   - `pca()` rotates rows of **B** by **V** (the SVD of X), the same for all responses.
   - Possible "coefficient biplot": rows of $\widehat{\mathbf{B}}(\lambda)$ (predictors) as
     points or vectors in the space of the responses, traced over λ. To be explored.

6. **`precision()` and `vif()`.** Little that's new (see §2.4). `vif()` is identical to the
   univariate case. `precision()` could report the multivariate measures
   (`log|Σ ⊗ V_λ|`, etc.) for completeness, but plotting them adds nothing beyond the univariate
   `plot.precision()`. Document this, rather than building new plots.

## 6. Beyond shared-λ ridge (scope: probably *not* in the first version)

Real multivariate shrinkage, which borrows strength across responses:

- **Brown & Zidek (1980)**, adaptive multivariate ridge regression: a ridge *matrix* in the
  response space, and a multivariate Hoerl–Kennard rule. In `statistics.bib` as
  `Brown-Zidek-1980` (*Annals of Statistics*, 8(1), 64–74). This is the most natural "next
  step" for genridge.
- **Curds and whey**: Breiman & Friedman (1997), "Predicting multivariate responses in multiple
  linear regression", *JRSS B*, 59(1), 3–54. *(Citation verified 2026-11-19, but not in
  statistics.bib.)* Their eq. (6.3) is exactly the shared-λ estimator
  $(\mathbf{X}^\top\mathbf{X} + \lambda\mathbf{I})^{-1}\mathbf{X}^\top\mathbf{Y}$, so it is the
  natural reference for the method itself.
- **Reduced-rank ridge regression**: Mukherjee & Zhu (2011), *Statistical Analysis and Data
  Mining*, 4(6). *(Not in statistics.bib; verify.)*
- **glmnet** with `family = "mgaussian"` (a group-lasso / elastic-net penalty across
  responses): a comparison point, like the existing `extra/glmnet-test.R` and
  `dev/boot-glmnet.R`.

Recommendation: ship the shared-λ version with good graphics first, and mention these in the docs
as extensions. The new plots (§5.3, §5.4) are where genridge can add something.

## 7. Example data

The example data need several responses and *collinear* predictors.

- `heplots::Rohwer`: 3 responses (`SAT`, `PPVT`, `Raven`), 5 paired-associate predictors.
  Collinearity is only modest: the largest VIF is 3.0 (`na`), and the correlations are up to
  0.72. It's familiar from heplots and the book, but the shrinkage may look unexciting.
- Other candidates to look for:
  - a multi-response version of a classic ridge dataset, e.g., Longley with several of its
    variables as responses (artificial);
  - an economics or chemometrics dataset with several outcomes;
  - simulated data with controlled collinearity, for the vignette.

Open: is there a natural dataset in heplots, candisc or mvinfluence with high VIFs?

**Search results (2026-11-19).** A web search plus checks in R found these. The full report,
which will be lost if the temp folder is cleared, is
`scratchpad/mlm-dataset-candidates.md`. The key numbers for `schooldata` were re-checked with
genridge's `ridge()`.

1. **`heplots::schooldata` (recommended main example).**
   - n = 70.
   - Responses: `reading`, `mathematics`, `selfesteem`, correlated 0.94–0.99.
   - Predictors: `education`, `occupation`, `visit`, `counseling`, `teacher`.
   - **VIFs: 7.4, 85.9, 61.6, 53.9, 1.6.** `occupation`, `visit` and `counseling` correlate
     at r = 0.98–0.99.
   - Clear wrong-sign OLS coefficients that ridge corrects:
     - `visit` is −5.3 for reading and −9.9 for math at λ = 0, and positive by λ = 5;
     - `occupation` is inflated (≈ 20–25) and shrinks to ≈ 4–7.
   - heplots already uses it for MLM and HE-plot examples, so it suits ridge HE plots
     (§5.4).
   - Caveats:
     - reading and math are almost collinear (r = 0.99), so ellipses in that plane will be
       thin;
     - several high-leverage or outlying cases (`robmlm()` downweights 12);
     - the variables are indices that are hard to interpret.
2. **`robustbase::pulpfiber`.**
   - n = 62; 4 responses (paper properties), 4 predictors (fibre properties).
   - VIFs up to 6.2; the responses correlate at 0.87–0.99.
   - Wrong-sign OLS coefficients for X1 and X3 that ridge corrects.
   - A milder second example.
3. **`mtcars`.**
   - Responses `mpg` and `qsec`; predictors `cyl`, `disp`, `hp`, `wt`, `drat`, `carb`.
   - Max VIF is 17.8.
   - Fine for short `@examples` with no extra dependency, though `qsec` as a response is a
     stretch.
4. Rejected, mostly because collinearity was too mild or the responses weren't a natural MLM:
   `Rohwer`, `Sake`, `NLSY`, `Hernior`, linnerud, `pls::oliveoil`, `rrr::tobacco` (Izenman's
   example, VIF only 4.1), `CASchools`, `Credit` and others.

## 8. Suggested implementation steps

1. Write `ridge_mlm_fit(Y, X, lambda, df, svd)` in `R/`, reusing `ridge.default()`'s
   centering, scaling, df ↔ λ solving, `V_λ` and criteria code. Add pooled GCV from the draft,
   converted to genridge's λ scale.
2. Choose the object structure (§4), with `print`, `coef` and `vcov` methods, and `as.ridge()`.
3. Interface: `ridge.mlm()`, plus a matrix-response branch in `ridge.formula()` and
   `ridge.default()` (§3). Remove the `mlm` stop in `ridge.lm()`.
4. Tests (the package has no `tests/` yet; add testthat?):
   - each response column equals the univariate `ridge()`;
   - `kronecker(Sigma, V_λ)` diagonal blocks equal the univariate `cov`;
   - λ ↔ VARshrink mapping (λ_g = n λ_v) on a fixed example;
   - the formula, `mlm` and matrix interfaces agree.
5. Plots: delegation (§5.1), `traceplot` (§5.2), cross-response ellipses (§5.3). Then HE-style
   (§5.4) as a separate step.
6. Docs:
   - roxygen in markdown;
   - an example using Rohwer or a better dataset;
   - possibly a short vignette, *Multivariate ridge regression*;
   - NEWS; version 0.9.0, since this is a new feature set.
7. Optionally, a section in Vis-MLM-book ch. 9.

## 9. Decisions (MF, 2026-11-19)

These replace the open questions that were here.

1. **Class:** `c("mlmridge", "ridge")`, inheriting from `ridge` for now.
2. **`coef`:** a list of p × q `B` matrices, one per λ, as in VARshrink. Add a helper to
   convert it to a 3-d array `[lambda, predictor, response]`.
3. **Fitting:** try to generalize `ridge.default()` to accept a matrix `y`, rather than
   writing a separate fit function.
4. **HE-style ridge plots (§5.4):** a great and novel contribution, but they need further
   thought and design. Not in the first version; treat them as a separate design task.
5. **Example dataset (§7):** MF will look for one that shows multivariate shrinkage. This may
   need a web search.
6. **Tests:** deferred.

**What follows from decisions 1 and 2.** Because the class inherits from `ridge`, every
existing `*.ridge` method will be dispatched for an `mlmridge` object. Those methods expect
`coef` to be a λ × p *matrix*, so with a list they will fail or give wrong answers. So:

- `mlmridge` needs its own method, or an explicit stop, for each inherited generic that
  touches `coef` or `cov`:
  - `coef`, `vcov`, `print`, `best`;
  - `plot`, `pairs`, `plot3d`, `traceplot`;
  - `pca`, `biplot`, `precision`, `vif`.
- The simplest first version: each of these methods takes `response =` and delegates to the
  univariate method through `as.ridge(x, response)` (§4, §5.1). Methods that don't depend on
  the response, such as `vif()`, can use any response.
- Audit which `*.ridge` methods read `x$coef` and `x$cov` directly before deciding; some may
  only need `x$lambda`, `x$df` or `x$svd.*`.

**What follows from decision 3.** `ridge.default()` is currently univariate throughout:
- `y` is centered as a vector;
- `MSE` and the HKB/LW criteria are computed with `lm.fit()`;
- `coef[i, ]`, `mse[i]` and GCV are per-λ scalars or vectors.

Generalizing it means:
- treating `y` as an n × q matrix, with q = 1 the current case;
- computing `Sigma` (q × q) in place of `MSE`;
- storing `cov` as `V_λ`, *without* the `MSE` factor, in the multivariate case;
- computing criteria per response plus a pooled GCV.

**The univariate output must stay byte-for-byte the same.** Before changing anything, save
reference results from the current `ridge()` for longley, prostate, Manpower and Acetylene, and
compare after the change. This stands in for the deferred tests.

## 10. Audit: what inherited `*.ridge` methods do with an `mlmridge` object (2026-11-19)

**Method.** I grepped `R/` for uses of `x$coef`, `x$cov`, `coef()`, `vcov()` and the other
components. Then I built a **mock** `mlmridge` object with the decided structure and called
every inherited method on it:
- class `c("mlmridge", "ridge")`;
- `coef` = a list of p × q `B` matrices, one per λ;
- `cov` = a list of p × p `V_λ`, *without* the Σ or MSE factor;
- `Sigma` = q × q;
- data: longley, 5 predictors, 2 responses.

| Method | Reads | Result on mock | Needed for `mlmridge` |
|---|---|---|---|
| `coef.ridge()` | `x$coef` | returns the list | fine as is (decision 2); add `coef(x, response=)` and the 3-d array helper |
| `vcov.ridge()` | `x$cov` | returns the list of `V_λ` | its meaning changes (no Σ factor). Own method: `vcov(x, response=)` → `Sigma[k,k] * V_λ`; `vcov(x)` → the list of `kronecker(Sigma, V_λ)`? |
| `print.ridge()` | `coef(x)` | runs; prints the raw list | own method (per-response coefficient tables, or one stacked table) |
| `best.ridge()` | `criteria`, `svd.D` | runs | depends on how criteria are stored (per response + pooled GCV). The `df` column is right, since `svd.D` comes from X only |
| `traceplot()` | `coef(x)` as a λ × p matrix; `kHKB`, `kLW` | **ERROR**: incorrect number of dimensions | own method (§5.2); `kHKB`/`kLW` become per response |
| `plot.ridge()` | `x$coef[which.lambda, variables]`, `x$cov`, `dimnames(x$coef)[[2]]` | **ERROR**: "There are only 0 predictor variables" | delegate per response (§5.1). Ellipses need `Sigma[k,k] * V_λ` |
| `pairs.ridge()` | `x$coef[, i]`, `apply(x$coef, 2, range)` | **ERROR**: dim(X) must have a positive length | delegate per response |
| `plot3d.ridge()` | same pattern as `plot.ridge()` | not run (rgl); will fail the same way | delegate per response |
| `precision.ridge()` | `x$cov`, `coef()` as a matrix, `ncol(b)` | **ERROR**: non-numeric argument to binary operator | delegate per response. A multivariate version would use `kronecker(Sigma, V_λ)` (§2.4) |
| `vif.ridge()` | `vcov()`, `colnames(coef())`, `rownames(coef())` | runs, but the **row and column names are lost** | values are correct: VIF uses `cov2cor()`, so it doesn't depend on Σ, and it equals the univariate VIF. Only the names need fixing; no `response=` argument needed |
| `pca.ridge()` | `x$coef %*% x$svd.V`, `x$cov` | **ERROR**: requires numeric matrix | own method: rotate each `B` (`t(V) %*% B`, i.e., rows of B) and each `V_λ` |
| `biplot.ridge()` / `biplot.pcaridge()` | `coef(x)`, `svd.V`, `svd.D` | **ERROR**: "There are only 0 predictor variables" | delegate per response, or a new coefficient biplot (§5.5) |

**Summary:**
- 7 of the 12 methods give an error; 5 run.
- Of the 5 that run:
  - `coef` and `vcov` return lists, which is by design;
  - `print` is untidy;
  - `best` depends on how the criteria are stored;
  - `vif` silently loses its names.
- There are no silently *wrong numbers*, as long as `cov` holds `V_λ` and every method that
  uses its scale (`plot`, `pairs`, `plot3d`, `precision`) goes through a per-response
  accessor that multiplies by `Sigma[k, k]`.

**Points to watch:**
- **Class order after `pca()`.** `pca.ridge()` sets `class(result) <- c("pcaridge", class(x))`,
  which gives `c("pcaridge", "mlmridge", "ridge")`. So the `*.pcaridge` methods (`plot`,
  `plot3d`, `biplot`) would be dispatched *before* any `*.mlmridge` method. `pca.mlmridge()`
  needs to produce something like `c("mlmpcaridge", "mlmridge", "pcaridge", "ridge")`, or
  `pca()` could be supported only per response.
- **`cov` changes meaning.** For `ridge` it is the full covariance (MSE × `V_λ`); for
  `mlmridge` it is the X part only. Every method that reads `x$cov` directly (`plot.ridge`,
  `plot3d.ridge`, `precision.ridge`, `pca.ridge`) must reach it through
  `as.ridge(x, response)`, never directly. Alternatively, store it under a different name
  (e.g., `covX`) and leave `cov` absent, so a stray inherited call fails loudly instead of
  quietly using the wrong scale. **Recommendation: use `covX`.**
- **`kHKB`, `kLW`, `kGCV` and `criteria`** are read directly by `traceplot`, `best`,
  `plot.vif.ridge` and `vif.ridge`. Decide their multivariate form (per-response vectors plus a
  pooled GCV) before writing those methods.

**Minimal first version** (enough that nothing errors):
1. `as.ridge(x, response)`.
2. `coef`, `vcov` and `print` methods.
3. A `response=` argument, defaulting to the first response with a message, on `traceplot`,
   `plot`, `pairs`, `plot3d`, `precision`, `pca` and `biplot`, delegating through `as.ridge()`.
4. A names fix for `vif`.

`best` comes once the criteria are settled. Multi-response `traceplot` (§5.2) and the
cross-response ellipses (§5.3) come after that.
