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

Class: **`"mlmridge"`**, a separate class that does **not** inherit from `ridge`. See §9
decision 1, revised after the audit in §10.

`ridge()` stays the single entry point. With a matrix response it returns an `mlmridge`
object, much as `lm()` returns an `mlm` for a matrix response. Generics that have no
`mlmridge` method yet fail with R's standard "no applicable method for '...' applied to an
object of class 'mlmridge'". That is clearer than the misleading errors inherited `ridge`
methods would give.

Components. The names follow `ridge` where the meaning is the same; since the classes are
separate, a name can take a multivariate meaning without risk:

- `lambda`, `df`: as now. They depend only on X, so they are shared by all responses.
- `coef`: a **list** of p × q `B` matrices, one per λ (decision 2). A helper,
  e.g. `coef_array(x)`, returns the 3-d array `[lambda, predictor, response]`.
- `cov`: a list of p × p `V_λ`, the shared X part *without* Σ. There's no need to rename it
  `covX`, because no inherited `ridge` method will ever read it.
- `Sigma`: the q × q residual covariance from the OLS fit, **E** / (n − p − 1).
- `mse`: a matrix `[lambda, response]`.
- `scales`, `svd.D`, `svd.U`, `svd.V`: as now; they depend only on X.
- Criteria: see **§4a** below.

**Key helper: extract a univariate `ridge` object for one response**, e.g.
`as.ridge(x, response = "Y1")`.
- It returns a genuine `ridge` object, with:
  - `coef` = the λ × p matrix for that response;
  - `cov` = `Sigma[k, k] * V_λ`;
  - `mse`, `kHKB`, `kLW`, `GCV`, `kGCV` and `criteria` for that response;
  - `lambda`, `df`, `scales` and `svd.*` unchanged.
- So **every existing `ridge` method applies per response**: `traceplot()`, `plot()`,
  `pairs()`, `plot3d()`, `precision()`, `vif()`, `pca()`, `biplot()`.
- This works the same with or without inheritance, and it is where the code reuse comes from.
- Step 1 of the implementation.

Basic methods:

- `coef(x, response = NULL)`: the list of B matrices, or one response's λ × p matrix.
- `vcov(x, response = NULL)`: a list of pq × pq `kronecker(Sigma, V_λ)`, or one response's list
  of `Sigma[k, k] * V_λ`.
- `print()`: λ and df, then a coefficient table per response.
- Optionally, `is_ridge_fit()` or a documentation note, if other code needs to recognize
  "either kind of ridge fit". (The one argument for inheriting was `inherits(x, "ridge")`.)

## 4a. Criteria layout (agreed 2026-11-19)

**The rule.** Per-response values go in fields with the **same names as in `ridge`**, each
widened by one dimension for responses. Pooled values get **new names**. Then
`as.ridge(x, k)` is just a subset, and the univariate methods (`traceplot`, `best`, `vif`,
`plot.vif.ridge`) work on its result unchanged.

| Field | Shape (q responses, m λ values) | Contents |
|---|---|---|
| `kHKB` | named vector, length q | HKB estimate per response, from that response's OLS fit (the univariate formula) |
| `kLW` | named vector, length q | Lawless–Wang estimate per response |
| `GCV` | m × q matrix | GCV per λ and response (the univariate formula, column by column) |
| `kGCV` | named vector, length q | the λ minimizing each response's GCV |
| `criteria` | q × 3 matrix | rows = responses; columns `kHKB`, `kLW`, `kGCV` |
| `GCV.pooled` | vector, length m | one GCV across responses (see below) |
| `kGCV.pooled` | number | the λ minimizing `GCV.pooled`: **the shared λ to recommend** |

`as.ridge(x, k)` sets:
- `kHKB = x$kHKB[k]`, `kLW = x$kLW[k]`, `kGCV = x$kGCV[k]`;
- `GCV = x$GCV[, k]`;
- `criteria = x$criteria[k, ]`.

That gives exactly the univariate structure.

**Pooling GCV.** In genridge the GCV denominator, $(n - \mathrm{df}(\lambda))^2$, depends only
on X, so it is the same for every response. Pooling therefore combines the residual sums of
squares $\mathrm{RSS}_k(\lambda)$, and their scale matters. There are three options, set with an
argument `pool = c("scaled", "raw", "mahalanobis")`:

- **`"scaled"` (default):**
  $\sum_k \mathrm{RSS}_k(\lambda) / \hat\sigma^2_{kk} \,/\, (n - \mathrm{df})^2$, where
  $\hat\sigma^2_{kk}$ is the OLS residual variance (the diagonal of `Sigma`).
  - Each response contributes relative to its own noise level.
  - The result doesn't depend on the units of the responses.
- **`"raw"`:** $\sum_k \mathrm{RSS}_k(\lambda) / (n - \mathrm{df})^2$, as in VARshrink.
  - The response with the largest variance dominates.
  - The result depends on units. Included for comparison with VARshrink.
- **`"mahalanobis"`:**
  $\mathrm{tr}(\hat{\boldsymbol\Sigma}^{-1} \mathbf{R}(\lambda)^\top \mathbf{R}(\lambda)) / (n - \mathrm{df})^2$.
  - It also uses the residual *covariances*, which is more principled.
  - It is unstable when responses are nearly collinear (in schooldata, reading and math have
    r = 0.99), so it is an option, not the default.

**HKB and LW:** per response only. There is no ad-hoc pooling (e.g., averaging), since it has
no theoretical basis. The principled multivariate analogue is Brown & Zidek's (1980)
multivariate Hoerl–Kennard rule, which belongs to §6 (later extensions).

**How the methods use the criteria:**
- **`best.mlmridge()`** returns a long data frame with columns `response`, `criterion`
  (HKB, LW, GCV), `k` and `df`:
  - one row per response and criterion;
  - a final `(pooled)` row for GCV.
  
  `df` comes from `svd.D`, which is shared, so it is computed exactly as in `best.ridge()`.
- **`traceplot.mlmridge()`** draws one panel per response. Each panel marks that response's
  HKB and LW, and every panel marks the pooled-GCV λ at the same place.
- **`print.mlmridge()`** shows the pooled-GCV λ and the q × 3 `criteria` table.

**To document:** `kGCV` and `kGCV.pooled` are the best values *among the λ values supplied*,
as now. HKB and LW are closed-form estimates.

## 5. Plot methods: ideas, roughly in order of value

1. **Per-response delegation (cheap, do first).** `plot(x, response = "Y1")` etc. call the
   univariate methods through `as.ridge()`. A `response = NULL` default could loop over the
   responses, or make a multi-panel `traceplot`.
   *Done 2026-11-19 (commit `3ed88a0`). Currently, omitting `response` uses the first
   response, with a message.*

2. **`traceplot.mlmridge()`.** Coefficient paths vs λ or df:
   - either one panel per response, which is the easiest to read, sharing the λ axis and the
     HKB/LW/GCV reference lines;
   - or one panel with colour = predictor and line type = response.

   **What `traceplot(sridge)` with no `response` should do (MF, 2026-11-19).** Rather than
   picking the first response with a message, it should show **all** responses.
   - Users can already do this by hand, because the default `main` is the response name, so
     each panel is labelled correctly:
     - a loop: `for (r in ynames) traceplot(sridge, response = r)`, as in
       `dev/mlmridge-plot-examples.R`;
     - a map: `invisible(lapply(ynames, \(r) traceplot(sridge, response = r)))`.
   - Proposal: make `response = NULL` the default for `traceplot.mlmridge()`. With `NULL`,
     loop over all responses:
     - set up a multi-panel layout (`par(mfrow = n2mfrow(q))`, restored on exit) unless the
       user has already set one, or add a `layout = TRUE/FALSE` argument;
     - draw one panel per response, each with its own HKB/LW lines;
     - optionally mark the pooled-GCV λ in every panel.
   - *Done 2026-11-19 (commit `981867e`):*
     - `traceplot.mlmridge(x, response = NULL, main, layout = TRUE, ...)`;
     - `response` can be several names or numbers;
     - layout `rev(n2mfrow(q))` (1 × 3, 2 × 2, 2 × 3, ...), restored on exit;
     - `layout = FALSE` uses the user's own `par(mfrow)`;
     - `main` is recycled over the panels.
     
     Not done: a pooled-GCV reference line in each panel. `traceplot.ridge()` draws only the
     HKB/LW lines, not GCV.
   - The same default may suit `plot()` and `precision()` (all responses, one panel each).
   - `pairs()` (already a matrix of panels) and `plot3d()` (rgl) should keep a single
     response.

3. **Cross-response coefficient ellipses (new, and genuinely multivariate).** For a single
   predictor *j*, plot the q-vector $\widehat{\mathbf{b}}_{j\cdot}(\lambda)$ for two responses,
   with covariance ellipses $[\mathbf{V}_\lambda]_{jj}\,\boldsymbol{\Sigma}$ across λ.
   - Because of the Kronecker structure, the ellipses all have the **same shape and
     orientation** (that of Σ). Only their size and centre change with λ. That's a nice
     teaching point.
   - A `pairs()` version would take all response pairs for one predictor.

   *Done 2026-11-19 (prototype `dev/cross-ellipses-prototype.R`, then in the package).*
   - **Interface (MF: option 2):** `plot(x, predictor = , responses = 1:2, ...)` on
     `plot.mlmridge()`, with the arguments of `plot.ridge()` (`radius`, `which.lambda`,
     `labels`, `col`, `fill`, `ref`, …), and internal `plot_cross_response()`. Without
     `predictor`, `plot()` keeps the within-response plot (`response =`).
   - The two modes have different meanings (axes = predictors vs axes = responses). This is
     explained under `@details` in `?mlmridge-plots`, including how the plot differs from
     `heplots::coefplot()`:
     - here each predictor is traced over λ;
     - the ellipses are `radius`-scaled sampling-variability ellipses around a biased estimate,
       not confidence regions, so covering the origin isn't a test.
   - Findings on schooldata:
     - `visit` crosses zero for both responses together;
     - `counseling` moves *away* from the origin (ridge shrinks the norm, not each
       coefficient);
     - at λ = 5 the collinear predictors converge into a cluster (ridge's grouping effect).
   - **Note: a `heplots::coefplot()` method for `mlmridge`** (option 1, not done). This would
     be a thin wrapper around `plot_cross_response()`, for users who know `coefplot()` from
     heplots.
     - Because heplots is only in Suggests, the method would use delayed registration
       (`@exportS3Method heplots::coefplot`, i.e. `S3method(heplots::coefplot, mlmridge)`).
       It is active only when heplots is loaded, so `coefplot(sridge)` needs
       `library(heplots)`.
     - It should take `coefplot()`'s arguments where they make sense: `variables` (the two
       responses), `parm` (predictors), `level` (convert to `radius`, perhaps with a note that
       it isn't a confidence level here), `fill`, `labels`.
     - Its documentation must repeat the difference in meaning.
     - Alternatively, it could live in heplots as `coefplot.mlmridge()`, with genridge in
       heplots' Suggests. Decide when doing it.
   - **Idea (held, 2026-11-19): make the several-predictor display easier to read.**
     - **The plot:** `plot(sridge, predictor = 1:5, which.lambda = c(1, 6))` (gallery
       `cross-all.png`).
     - **What it shows:** each predictor's coefficients (reading, math) at λ = 0 and λ = 5,
       with ellipses. The key message is the **grouping effect**: at λ = 5 the near-duplicate
       predictors (occupation, visit, counseling, education; r = 0.98–0.99) converge to a
       common value. Occupation comes *down* from ≈ 20 and visit comes *up* from negative
       values, rather than each shrinking independently towards 0.
     - **Why it's hard to read now:**
       1. only the selected λ values are joined, so the path through the intermediate λ is
          hidden;
       2. the small red ellipses at λ = 5 overprint each other in the cluster;
       3. labels appear only at the start (λ = 0), and they overlap (education/counseling);
       4. there is no legend for the λ colours.
     - **What could be done:**
       - draw each predictor's **full path** through all λ (a thin line, perhaps with small
         points), but **ellipses only at `which.lambda`**: a new argument such as
         `path = TRUE` / `path.lambda = "all"`;
       - optional **end labels** (`label.end = TRUE`) at the largest λ, and/or a small
         offset (`pos`) for label clashes;
       - a **λ legend** (`legend = TRUE`) when colours encode λ;
       - colour by *predictor* with line type or fill alpha by λ, which reads better when
         there are several predictors;
       - a "zoom" option or inset for the cluster region at large λ;
       - point to it in the vignette with an explanation, not as a default display.
     - It may be easier to redo in tinyplot (`by = predictor` gives automatic legends and
       palettes; see `dev/tinyplot-notes.md`).

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
   - **Criteria favour little shrinkage** (first run of the new code, 2026-11-19):
     - HKB is 0.15 for reading, 0.15 for math and 0.03 for self-esteem; LW is 0.25, 0.30 and
       0.02. GCV per response and pooled both have their minimum near λ = 0.1.
     - The wrong-sign `visit` coefficients only correct themselves at λ ≈ 1–5, well past these
       "optimal" values.
     - This is a teaching point for the example: the prediction-optimal λ and the λ that gives
       interpretable coefficients differ.
     - Use a fine grid at the low end, e.g. `lambda = c(0, 0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10)`.
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

1. **Class:** a separate class `"mlmridge"` that does **not** inherit from `ridge`.
   *Revised 2026-11-19.* The original decision was `c("mlmridge", "ridge")`, inheriting from
   `ridge`. The audit in §10 showed that inheriting gives almost no reuse: no `*.ridge`
   method works unchanged, and the real reuse comes from `as.ridge()`, which doesn't need
   inheritance. It also has costs:
   - misleading errors from inherited methods;
   - two traps: `cov` meaning different things, and the class order after `pca()`;
   - every future `*.ridge` method would silently apply to `mlmridge` too, the same problem
     stats has with `mlm` inheriting from `lm`.
2. **`coef`:** a list of p × q `B` matrices, one per λ, as in VARshrink. Add a helper to
   convert it to a 3-d array `[lambda, predictor, response]`.
3. **Fitting:** try to generalize `ridge.default()` to accept a matrix `y`, rather than
   writing a separate fit function.
4. **HE-style ridge plots (§5.4):** a great and novel contribution, but they need further
   thought and design. Not in the first version; treat them as a separate design task.
5. **Example dataset (§7):** MF will look for one that shows multivariate shrinkage. This may
   need a web search.
6. **Tests:** deferred.

**What follows from decisions 1 and 2.** No `*.ridge` method is dispatched for an `mlmridge`
object, so each generic we want to support needs an `*.mlmridge` method. Unsupported generics
fail with R's standard "no applicable method" error. The methods to provide, in order:

- `coef`, `vcov`, `print`; then `best`, once the criteria layout is settled;
- `traceplot`, `plot`, `pairs`, `plot3d`, `precision`, `vif`, `pca`, `biplot`.

The simplest first version: each method takes `response =` and delegates to the univariate
method through `as.ridge(x, response)` (§4, §5.1). Methods that don't depend on the response,
such as `vif()`, can use any response. Real multi-response methods (§5.2, §5.3) replace the
delegation later, one generic at a time.

The list form of `coef` (decision 2) fits a separate class better: no existing code expects a
matrix there.

**What follows from decision 3.** `ridge.default()` is currently univariate throughout:
- `y` is centered as a vector;
- `MSE` and the HKB/LW criteria are computed with `lm.fit()`;
- `coef[i, ]`, `mse[i]` and GCV are per-λ scalars or vectors.

Generalizing it means:
- treating `y` as an n × q matrix, with q = 1 the current case;
- computing `Sigma` (q × q) in place of `MSE`;
- storing `cov` as `V_λ`, *without* the `MSE` factor, in the multivariate case;
- computing criteria per response plus a pooled GCV;
- returning class `"ridge"` when q = 1, and `"mlmridge"` when q > 1.

**The univariate output must stay byte-for-byte the same.** Before changing anything, save
reference results from the current `ridge()` for longley, prostate, Manpower and Acetylene, and
compare after the change. This stands in for the deferred tests.

## 10. Audit: what inherited `*.ridge` methods do with an `mlmridge` object (2026-11-19)

*This audit led to revising decision 1 (§9): `mlmridge` is now a separate class, so none of
these inherited calls will happen. The table still records which `ridge` methods depend on
which components, and that is what each `*.mlmridge` method or the `as.ridge()` delegation has
to supply.*

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

**Points to watch.** These were problems under inheritance; with a separate class, the first
two go away.
- *Resolved:* **class order after `pca()`.** Under inheritance, `pca.ridge()` would have
  produced `c("pcaridge", "mlmridge", "ridge")`, sending `plot`, `plot3d` and `biplot` to the
  `*.pcaridge` methods first. With a separate class, `pca.mlmridge()` either delegates per
  response, returning an ordinary `pcaridge` object, or later returns its own
  `c("mlmpcaridge", "mlmridge")`.
- *Resolved:* **`cov` changing meaning.** Only `*.mlmridge` methods and `as.ridge()` ever
  read an `mlmridge` object's `cov`, so it can hold `V_λ` under the same name. The `covX`
  rename is no longer needed.
- *Resolved (§4a):* **`kHKB`, `kLW`, `kGCV` and `criteria`.** The univariate methods read
  these directly (`traceplot`, `best`, `plot.vif.ridge`, `vif.ridge`). The layout in §4a keeps
  the same names, widened per response, so `as.ridge()` fills them by subsetting. Pooled GCV
  goes in `GCV.pooled` and `kGCV.pooled`.

**Minimal first version** (enough that every supported generic works):
1. `as.ridge(x, response)`, returning a genuine `ridge` object.
2. `coef`, `vcov` and `print` methods for `mlmridge`.
3. `traceplot`, `plot`, `pairs`, `plot3d`, `precision`, `pca` and `biplot` methods for
   `mlmridge`, each with a `response=` argument (defaulting to the first response, with a
   message) and delegating through `as.ridge()`.
4. A `vif.mlmridge()` that returns the shared VIFs, with predictor and λ names.

`best` comes once the criteria are settled. Multi-response `traceplot` (§5.2) and the
cross-response ellipses (§5.3) come after that.

## 11. Vignette: multivariate ridge regression (planned 2026-11-19)

genridge has **no vignettes yet**. The multivariate extension is a natural first one, and the
example gallery `dev/mlmridge-plot-examples.R` (code next to each plot) is a good start: it
mainly needs descriptive text between the chunks.

**Working title:** *Multivariate ridge regression with genridge*
(file `vignettes/mlm-ridge.Rmd`).

**Outline:**

1. **Introduction.**
   - Collinearity when there are several responses.
   - Ridge with one λ shared by all responses:
     $\widehat{\mathbf{B}}(\lambda) = (\mathbf{X}^\top\mathbf{X} + \lambda\mathbf{I})^{-1}\mathbf{X}^\top\mathbf{Y}$
     (Breiman & Friedman, 1997, eq. 6.3).
   - Each column is the univariate ridge for that response. What's new is the shared λ, the
     pooled criteria and the cross-response covariance, $\boldsymbol\Sigma \otimes \mathbf{V}_\lambda$ (§2).
2. **Data: `heplots::schooldata`.**
   - VIFs (up to 86) and the near-duplicate predictors.
   - The `lm()` MLM fit and its wrong-sign coefficients (`visit`).
   - Possibly an HE plot of the OLS fit, to link to heplots.
3. **Fitting.**
   - `ridge(cbind(reading, mathematics, selfesteem) ~ ..., lambda = )`, using a fine grid at
     the low end.
   - `print()`.
   - Criteria per response, the pooled GCV and the `pool =` options.
   - The three interfaces: formula, `mlm` fit and matrices.
4. **Coefficient paths.** `traceplot()` for all responses, via the §5.2 default or a loop/map.
5. **Bivariate ridge traces.**
   - `plot()` for occupation × visit, showing the sign correction.
   - Ellipses for one response have the same shape across responses, scaled by
     $\sigma^2_{kk}$ (the Kronecker structure).
   - `pairs()` for one response.
6. **PCA/SVD view.** `pca()` with `plot()`/`biplot()`: shrinkage acts mostly on the smallest
   dimensions.
7. **VIF and precision.**
   - VIFs are the same for every response.
   - `precision()` curves have the same shape and differ only in level.
8. **Choosing λ: a teaching point (§7).**
   - HKB, LW and GCV favour small λ (≈ 0.1–0.3).
   - The wrong-sign coefficients only correct themselves at λ ≈ 1–5.
   - So the λ that's best for prediction and the λ that gives interpretable coefficients
     differ.
9. **Simulated example with known B.**
   - Longley X, as agreed. Choose a true **B** (2–3 responses with different coefficient
     patterns) and Σ with moderate correlation between responses; set a seed.
   - Show $\|\widehat{\mathbf{B}}(\lambda) - \mathbf{B}\|$ (bias, variance, MSE) against λ, and
     where GCV/HKB/LW fall.
   - This is the one thing a real dataset can't show.
10. **Further directions.**
    - Cross-response coefficient ellipses (§5.3).
    - Ridge HE plots (§5.4), once designed.
    - Brown & Zidek's multivariate ridge and the other methods in §6.

**Infrastructure** (none exists yet):

- `vignettes/` folder; `DESCRIPTION`: `VignetteBuilder: knitr`, and `knitr` and `rmarkdown`
  in Suggests. heplots is already in Suggests (added with `as.ridge()`).
- **Format:** R Markdown with `rmarkdown::html_vignette`, which CRAN builds with no extra
  system dependencies. Quarto vignettes are possible but need Quarto on the CRAN machines.
- **References:** a project-local `vignettes/references.bib`, copied from the canonical files
  (don't edit those):
  - Friendly (2013), from vita.bib;
  - Brown & Zidek (1980), from statistics.bib;
  - Breiman & Friedman (1997), not yet in statistics.bib: propose an entry to MF;
  - Hoerl & Kennard (1970), Golub, Heath & Wahba (1979), Charnes et al. (1981) for
    schooldata.
- **Speed and size:** keep figures modest. Skip `plot3d()`, which needs rgl, or show a static
  snapshot; CRAN vignette build time matters.
- **pkgdown:** the vignette appears under *Articles*; add it to `_pkgdown.yml` if an
  `articles:` section is used.
- **When:** after §5.2 (the all-responses `traceplot` default) and preferably §5.3
  (cross-response ellipses), so the vignette shows the multivariate plots, not only
  per-response delegation. Ridge HE plots can come in a later update of the vignette.

**Optional second vignette:** a univariate *Generalized ridge trace plots* article, largely
reusing the README (longley, prostate and Manpower examples). The README is long, so moving
its detail into a vignette would also shorten it.
