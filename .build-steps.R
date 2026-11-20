# Release build steps for genridge (run interactively, a step at a time).
# Adapted from vcdExtra/heplots .build-steps.R.

# keep rgl from popping up windows during examples/vignettes
Sys.setenv(RGL_USE_NULL = TRUE)

# build the README.md if it is older than README.Rmd
if (!file.exists("README.md") || file.mtime("README.Rmd") > file.mtime("README.md")) {
  devtools::build_readme()
}

# spellcheck -- fix any issues before running remote checks
(words <- devtools::spell_check())

# Check package locally
devtools::check(args = "--as-cran")

# Check on win-builder (CRAN's Windows servers); results come by email
devtools::check_win_release()
devtools::check_win_devel()

# Check on R-hub (GitHub Actions; workflow in .github/workflows/rhub.yaml)
# rhub::rhub_doctor()          # one-time: verify the GitHub PAT and workflow
rhub::rhub_check(platforms = c("linux", "macos-arm64", "windows"))
# results: https://github.com/friendly/genridge/actions

# reverse dependencies (none on CRAN as of 0.8.1)
devtools::revdep()

# build the pkgdown site
pkgdown::build_site()

# submit to CRAN
devtools::submit_cran()
