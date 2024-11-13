## Test environments
* local Windows107 install, R version 4.4.1 (2024-06-14 ucrt)
* win-builder R Under development (unstable) (2024-11-11 r87319 ucrt)

## R CMD check results
There were no ERRORs or WARNINGS or NOTEs

## Reverse dependencies
character(0)

## Comments
This was modest update release, improving documentation and adding new functionality.


### genridge 0.7.1 (2024-11-07)

o Added links to gentalk.pdf
o Improved README and fixed some examples
o `precision()` result gains a class "precision" in preparation for a plot method
o Implemented `plot.precision()` for plots shrinkage vs. variance, using various criteria
o Fixed warning from `ridge()` related to contrasts
o `ridge()` now collects the optimal `criteria` in a named list, which can be added in `plot.precision()`
