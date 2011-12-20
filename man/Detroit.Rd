\name{Detroit}
\alias{Detroit}
\docType{data}
\title{
Detroit Homicide Data for 1961-1973
}
\description{
The data set \code{Detroit} was used extensively in the book by Miller (2002)
on subset regression.
The data are
unusual in that a subset of three predictors can be found which gives a
very much better fit to the data than the subsets found from the Efroymson
stepwise algorithm, or from forward selection or backward elimination.
They are also unusual in that, as time series data, the assumption of
independence is patently violated, and the data suffer from problems
of high collinearity.

As well, ridge regression reveals somewhat paradoxical paths of
shrinkage in univariate ridge trace plots, that are more comprehensible
in multivariate views.

}
\usage{data(Detroit)}
\format{
  A data frame with 13 observations on the following 14 variables.
  \describe{
    \item{\code{Police}}{Full-time police per 100,000 population}
    \item{\code{Unemp}}{Percent unemployed in the population}
    \item{\code{MfgWrk}}{Number of manufacturing workers in thousands}
    \item{\code{GunLic}}{Number of handgun licences per 100,000 population}
    \item{\code{GunReg}}{Number of handgun registrations per 100,000 population}
    \item{\code{HClear}}{Percent of homicides cleared by arrests}
    \item{\code{WhMale}}{Number of white males in the population}
    \item{\code{NmfgWrk}}{Number of non-manufacturing workers in thousands}
    \item{\code{GovWrk}}{Number of government workers in thousands}
    \item{\code{HrEarn}}{Average hourly earnings}
    \item{\code{WkEarn}}{Average weekly earnings}
    \item{\code{Accident}}{Death rate in accidents per 100,000 population}
    \item{\code{Assaults}}{Number of assaults per 100,000 population}
    \item{\code{Homicide}}{Number of homicides per 100,000 of population}
  }
}

\details{
The data were orginally collected and discussed by Fisher (1976) but the complete dataset first 
appeared in Gunst and Mason (1980, Appendix A). 
Miller (2002) discusses this dataset throughout his book, but doeesn't state clearly
which variables he used as predictors and
which is the dependent variable.   (\code{Homicide} was the dependent variable, and the
predictors were \code{Police} \dots \code{WkEarn}.) 
The data were obtained from StatLib.

A similar version of this data set, with different variable names appears
in the \code{bestglm} package.

}
\source{
\url{http://lib.stat.cmu.edu/datasets/detroit}
}
\references{
Fisher, J.C. (1976). Homicide in Detroit: The Role of Firearms. \emph{Criminology}, \bold{14}, 387--400.

Gunst, R.F. and Mason, R.L. (1980). \emph{Regression analysis and its application: A data-oriented approach}. 
Marcel Dekker.

Miller, A. J. (2002). \emph{Subset Selection in Regression}. 2nd Ed. Chapman & Hall/CRC. Boca Raton. 
}
\examples{
data(Detroit)

# Work with a subset of predictors, from Miller (2002, Table 3.14),
# the "best" 6 variable model
#    Variables: Police, Unemp, GunLic, HClear, WhMale, WkEarn
# Scale these for comparison with other methods

Det <- as.data.frame(scale(Detroit[,c(1,2,4,6,7,11)]))
Det <- cbind(Det, Homicide=Detroit[,"Homicide"])

# use the formula interface; specify ridge constants in terms
# of equivalent degrees of freedom
dridge <- ridge(Homicide~., data=Det, df=seq(6,4,-.5))

# univariate trace plots are seemingly paradoxical in that
# some coefficients "shrink" *away* from 0
traceplot(dridge, X="df")
vif(dridge)
pairs(dridge, radius=0.5)

plot3d(dridge, radius=0.5, labels=dridge$df)

# transform to PCA/SVD space
dpridge <- pca.ridge(dridge)
# not so paradoxical in PCA space
traceplot(dpridge, X="df")
biplot(dpridge, radius=0.5)

}
\keyword{datasets}
