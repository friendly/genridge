## Test environments
* local Windows10 install, R version 4.4.1 (2024-06-14 ucrt)
* win-builder R Under development (unstable) (2024-11-30 r87409 ucrt)

## R CMD check results
There were no ERRORs or WARNINGS or NOTEs

## Reverse dependencies
There are no reverse dependencies

## Comments
### genridge 0.8.0 (2024-11-30)

This is a major update to the package adding additional plotting methods

o Fleshed out and installed `plot.precision()`
o Illustrate `plot.precision()` in `README.Rmd`
o Added `diab` data, diabetes from CASI
o Fixed documentation for `plot.ridge()` so that `plot.pcaridge()` is documented explicitly.
o 'vif.ridge()` now returns a "vif.ridge" object list to prepare for a plot method.
o Added `plot.vif.ridge()` to plot VIFs vs shrinkage
o Added `norm.diff` = sqrt((b_OLS - b_ridge)^2) as a measure of shrinkage to `precision()`
o Added `plot.precision()` for plots of shrinkage vs. precision
o Added `best()` to display the optimal shrinkage criteria
