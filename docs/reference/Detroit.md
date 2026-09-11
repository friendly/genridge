# Detroit Homicide Data for 1961-1973

The data set `Detroit` was used extensively in the book by Miller (2002)
on subset regression. The data are unusual in that a subset of three
predictors can be found which gives a very much better fit to the data
than the subsets found from the Efroymson stepwise algorithm, or from
forward selection or backward elimination. They are also unusual in
that, as time series data, the assumption of independence is patently
violated, and the data suffer from problems of high collinearity.

As well, ridge regression reveals somewhat paradoxical paths of
shrinkage in univariate ridge trace plots, that are more comprehensible
in multivariate views.

## Format

A data frame with 13 observations on the following 14 variables.

- `Police`:

  Full-time police per 100,000 population

- `Unemp`:

  Percent unemployed in the population

- `MfgWrk`:

  Number of manufacturing workers in thousands

- `GunLic`:

  Number of handgun licences per 100,000 population

- `GunReg`:

  Number of handgun registrations per 100,000 population

- `HClear`:

  Percent of homicides cleared by arrests

- `WhMale`:

  Number of white males in the population

- `NmfgWrk`:

  Number of non-manufacturing workers in thousands

- `GovWrk`:

  Number of government workers in thousands

- `HrEarn`:

  Average hourly earnings

- `WkEarn`:

  Average weekly earnings

- `Accident`:

  Death rate in accidents per 100,000 population

- `Assaults`:

  Number of assaults per 100,000 population

- `Homicide`:

  Number of homicides per 100,000 of population

## Source

<https://lib.stat.cmu.edu/datasets/detroit>

## Details

The data were originally collected and discussed by Fisher (1976) but
the complete dataset first appeared in Gunst and Mason (1980, Appendix
A). Miller (2002) discusses this dataset throughout his book, but
doesn't state clearly which variables he used as predictors and which is
the dependent variable. (`Homicide` was the dependent variable, and the
predictors were `Police` ... `WkEarn`.) The data were obtained from
StatLib.

A similar version of this data set, with different variable names
appears in the `bestglm` package.

## References

Fisher, J.C. (1976). Homicide in Detroit: The Role of Firearms.
*Criminology*, **14**, 387–400.

Gunst, R.F. and Mason, R.L. (1980). *Regression analysis and its
application: A data-oriented approach*. Marcel Dekker.

Miller, A. J. (2002). *Subset Selection in Regression*. 2nd Ed. Chapman
& Hall/CRC. Boca Raton.

## Examples

``` r

data(Detroit)

# Work with a subset of predictors, from Miller (2002, Table 3.14),
# the "best" 6 variable model
#    Variables: Police, Unemp, GunLic, HClear, WhMale, WkEarn
# Scale these for comparison with other methods

Det <- as.data.frame(scale(Detroit[,c(1,2,4,6,7,11)]))
Det <- cbind(Det, Homicide=Detroit[,"Homicide"])

# use the formula interface; specify ridge constants in terms
# of equivalent degrees of freedom
dridge <- ridge(Homicide ~ ., data=Det, df=seq(6,4,-.5))

# univariate trace plots are seemingly paradoxical in that
# some coefficients "shrink" *away* from 0
traceplot(dridge, X="df")

vif(dridge)
#> Variance inflaction factors:
#>             Police  Unemp  GunLic  HClear  WhMale  WkEarn
#> 0.00000000  24.813  4.857  19.272  46.923  56.711  54.717
#> 0.03901941  20.119  3.870  11.362  33.072  35.390  31.883
#> 0.10082082  15.415  2.971   6.432  21.289  20.776  17.353
#> 0.20585744  10.776  2.208   3.597  12.171  11.258   8.750
#> 0.40331592   6.534  1.623   2.128   5.996   5.503   4.113
pairs(dridge, radius=0.5)


# \donttest{
plot3d(dridge, radius=0.5, labels=dridge$df)

# transform to PCA/SVD space
dpridge <- pca(dridge)

# not so paradoxical in PCA space
traceplot(dpridge, X="df")

biplot(dpridge, radius=0.5, labels=dpridge$df)

#> Vector scale factor set to  2.864948 

# show PCA vectors in variable space
biplot(dridge, radius=0.5, labels=dridge$df)

#> Vector scale factor set to  0.8511932 
# }

```
