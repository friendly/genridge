# Compare the current (work-in-progress) genridge against the saved univariate reference results.
# Run from the package root:  Rscript dev/check-ridge-ref.R
# Each fit and method result must be identical() to the reference; all.equal() is shown on failure.

suppressMessages(devtools::load_all(".", quiet = TRUE))
ref <- readRDS("dev/ridge-ref-0.8.1.rds")

# re-create the same fits with the current code
src <- readLines("dev/make-ridge-ref.R")
start <- grep("^data\\(prostate", src); end <- grep("^ref <- list", src) - 1
eval(parse(text = src[start:end]))

ok <- TRUE
for (nm in names(ref$fits)) {
  same <- identical(fits[[nm]], ref$fits[[nm]])
  cat(sprintf("%-18s fit %s", nm, if (same) "identical" else "DIFFERS"))
  if (!same) { ok <- FALSE; cat("\n"); print(all.equal(ref$fits[[nm]], fits[[nm]])) }
  for (m in names(ref$methods[[nm]])) {
    if (!identical(methods[[nm]][[m]], ref$methods[[nm]][[m]])) {
      ok <- FALSE
      cat(sprintf("\n   %s DIFFERS: ", m)); print(all.equal(ref$methods[[nm]][[m]], methods[[nm]][[m]]))
    }
  }
  cat("\n")
}
if (ok) cat("\nALL IDENTICAL to reference", ref$version, "\n") else cat("\nDIFFERENCES FOUND\n")
