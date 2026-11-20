# genridge tasks

Moved from `C:\Dropbox\R\TASKS-all.md` (2026-11-20). Track and check off items here.

No `issues/` folder; this file is the tracker. As of 2026-11-20 `dev/` has 16 files (including
this one) and `extra/` has 5. `dev/` is in `.Rbuildignore`.

TODOs:

- [X] **Convert roxygen to markdown** — done 2026-11-19, commit `9691d47` (pushed). roxygen2md
  did `\code`/`\link`/`\emph`/`\bold`; `\describe`, `\itemize`, `\enumerate` and `@return \item`
  lists hand-converted to markdown lists. Rendered help diffed against the old version;
  `check()` clean (0/0/0); pkgdown reference pages checked. Done (2026-11-20). CRAN release is
  tracked below.

- [ ] **Multivariate ridge regression** — adapted from `VARshrink::lm_multiv_ridge`, extends ridge
  regression to a multivariate `Y` (matrix response) rather than the current univariate-only
  `ridge()`. Genuinely unshipped, and thematically connects to the other MLM-focused packages
  (heplots, candisc, mvinfluence).
  File: `dev/ridge_mlm.R`. Planning notes (2026-11-19): `dev/ridge_mlm-plan.md` — design,
  plot-method ideas, open questions.
  **On master in 0.8.1, marked experimental** (merged 2026-11-19; branch `mlmridge` done).
  - Done: `ridge()` with a matrix response → class `mlmridge`; `as.ridge()`;
    `coef`/`vcov`/`print`; per-response plot and summary methods; `traceplot()` showing all
    responses; cross-response ellipses, `plot(x, predictor = )` (§5.3); pkgdown article
    "Multivariate ridge regression" (`vignettes/articles/mlm-ridge.Rmd`, with simulation).
  - Next:
    - [X] `best.mlmridge()` — done 2026-11-19 (long data frame + pooled GCV row; `response =`).
    - [X] Pre-existing bug: `ridge.formula()` failed with a single predictor (`X[, -Inter]` lacked
      `drop = FALSE`; inherited from `MASS::lm.ridge()`, which still has it); `vif.ridge()` too.
      Fixed 2026-11-19.
    - [X] HKB and LW use a factor (p - 2), so they were negative for p = 1 and 0 for p = 2.
      Done 2026-11-20: now `NA` when p < 3 (univariate and `mlmridge`), not marked on plots;
      documented in `?ridge` Details, NEWS, tests.
    - [X] pooled-GCV reference line in `traceplot.mlmridge()` (§5.2) — done 2026-11-19 (`GCV = TRUE`).
    - [X] cite VARshrink in the docs (§1) — done 2026-11-19 (`?ridge`, MLM article, NEWS).
    - [X] tests (testthat) — done 2026-11-19: `tests/testthat/`, 214 expectations; `R CMD check`
      clean. Found and fixed a bug (one-column matrix response with several λ).
    - [X] decide: should `plot()`/`precision()` default to all responses, like `traceplot()`?
      Decided 2026-11-20: no, keep them to one response (`response = 1`) for now; a faceted
      all-responses display may come later with tinyplot (see the tinyplot item).
  - Held / needs design: ridge HE plots (§5.4); clearer several-predictor ellipse display
    (§5.3); `heplots::coefplot` method; coefficient biplot (§5.5); Brown–Zidek (§6).

- [ ] **Consider `tinyplot` for the plot methods** (added 2026-11-19, this development cycle).
  - Notes and options in `dev/tinyplot-notes.md`; pilot in `dev/tinyplot-pilot.R`.
  - tinyplot 0.8.0 is base graphics underneath and uses only base R. It adds facets, themes
    and `by =` legends.
  - The pilot shows that a faceted `traceplot` for all responses works in one call, and that
    genridge-computed ellipses work with `type = "polygon"`.
  - **Decided (2026-11-19):** tinyplot for the *new* multivariate plots, keeping the univariate
    base methods (maybe with an `engine =` switch later). **Not yet:** finish the multivariate
    work in base graphics first.
  - To decide:
    - Imports or Suggests;
    - direct labels vs legends;
    - theme and palettes (sequential for λ).
  
- [X] **`ridge.lm()` method** — done 2026-11-19, commit `fd0ea3a` (pushed). Rewritten rather
  than taken from `dev/ridge-lm.R` (that draft re-evaluated `match.call()`, which can't work on a
  fitted model): uses `model.matrix()`/`model.response()`/offset of the fitted `lm`; errors for
  `glm`, `mlm`, weighted fits. Identical to `ridge.formula()` on 5 test cases (factors, `log()`,
  NAs, offset, no intercept). `dev/ridge-lm.R` now superseded.
  
- [ ] **Ideas from lmridge** (Imdad Ullah et al., 2018, R Journal; notes in `dev/lmridge-notes.md`).
  - [X] Cite it — done 2026-11-19: "Related packages" in the Getting started vignette;
    `?ridge` `@seealso`/`@references`. Entry `ImdadUllah-etal:2018` in `statistics.bib`.
  - [ ] Plug-in bias² / variance / MSE trace (compare with the true tradeoff in the MLM
    article's simulation).
  - [ ] More criteria in `best()`: PRESS, Cp, AIC/BIC; a single "criteria trace" plot.
  - [ ] Mark several k estimators on `traceplot()` (rug on the λ axis).
  - [ ] `predict()`, `fitted()`, `residuals()` methods for `ridge` / `mlmridge`.
  - [ ] Condition-number trace in `precision()`.

- [ ] **CRAN release of 0.8.1** — CRAN has 0.8.0 (2024-11). Status 2026-11-19: ready to submit.
  win-builder R-devel clean; local `--as-cran` 0/0/0; no revdeps; spelling clean; `cran-comments.md`
  written (commit `20238f8`, pushed). 2026-11-20: HKB/LW `NA` change; NEWS bug-fix list
  tidied. Re-checked after the HKB/LW change, all clean (2026-11-20): local `--as-cran` (R 4.6.2),
  win-builder release + devel, R-hub linux/macos-arm64/windows (new `.github/workflows/rhub.yaml`).
  DESCRIPTION Date bumped to 2026-11-20. Remaining: review `cran-comments.md`,
  `devtools::submit_cran()` (bump Date again if not submitted today). After acceptance: GitHub release tagged `0.8.1`.
  - Note for myutil: `release_cran_comments()` matches only `## Version x.y.z` NEWS headings
    (genridge uses `## genridge x.y.z`, so it took the whole NEWS) and `trimws()` flattens
    nested NEWS lists; Comments section was fixed by hand.

- [ ] **Bootstrap ridge regression** — exploring bootstrapped CIs/variability for ridge coefficients
  on the Longley data. No `boot`-related code anywhere in `R/`.
  File: `dev/boot-ridge.R`

Reference material — not cleanup, not TODO: `dev/lmridge-RJ-2018-060.pdf` (the R Journal paper on
`lmridge`; see the lmridge item above), `dev/lab10.Rmd` (an ISLR course-lab exercise,
not genridge development).

Clean-up candidates (confirmed superseded by shipped v0.7.1–0.8.0 features):

- [X] Done 2026-11-19, commit `db4c94b`: removed `dev/plot.vif.ridge.R` (shipped version is a
  superset, adds `Y = c("vif", "sqrt")`), `dev/test-plot-precision.R` (shipped `plot.precision()`
  covers it, including the spline-placed HKB/LW points) and `dev/ridge-lm.R` (superseded by
  `ridge.lm()`). The `diab` files were moved to `data-raw/` and kept.

Worked examples / comparisons, not classified: `extra/Longley-ex.R` (+ its rendered `.html`),
`extra/glmnet-test.R` (compares against `glmnet`), `extra/JHBener-ex.txt`, `extra/SO-question.txt`.
