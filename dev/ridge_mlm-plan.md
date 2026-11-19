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
  linear regression", *JRSS B*, 59(1), 3–54. *(Not in statistics.bib; verify the citation
  before using.)*
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

## 9. Open questions for MF

- Class design: should `mlmridge` inherit from `ridge` or not?
- `coef` as a 3-d array or as a list of B matrices?
- Generalize `ridge.default()` in place, or add a separate fit function?
- Are HE-style ridge plots (§5.4) in scope for the first version, or a later one?
- Example dataset (§7).
- Add a `tests/` folder with testthat as part of this work?
