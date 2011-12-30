\name{Manpower}
\alias{Manpower}
\docType{data}
\title{
Hospital manpower data
}
\description{
The hospital manpower data, taken from Myers (1990), table 3.8,
are a well-known example of highly collinear data to which
ridge regression and various shrinkage and selection methods
are often applied.

The data consist of measures taken at 17 U.S. Naval Hospitals
and the goal is to predict the required monthly man hours
for staffing purposes.
}
\usage{data(Manpower)}
\format{
  A data frame with 17 observations on the following 6 variables.
  \describe{
    \item{\code{Hours}}{monthly man hours (response variable)}
    \item{\code{Load}}{average daily patient load}
    \item{\code{Xray}}{monthly X-ray exposures}
    \item{\code{BedDays}}{monthly occupied bed days}
    \item{\code{AreaPop}}{eligible population in the area in thousands}
    \item{\code{Stay}}{average length of patient's stay in days}
  }
}
\details{
Myers (1990) indicates his source was 
"Procedures and Analysis for Staffing Standards Development: Data/Regression Analysis Handbook", 
Navy Manpower and Material Analysis Center, San Diego, 1979. 
}
\source{
Raymond H. Myers (1990). \emph{Classical and Modern Regression with Applications}, 2nd ed.,
PWS-Kent, pp. 130-133.
}
\references{
Donald R. Jensen and Donald E. Ramirez (2012).
Variations on Ridge Traces in Regression,
\emph{Communications in Statistics - Simulation and Computation},
41 (2), 265-278.

}
\seealso{
\code{\link[bestglm]{manpower}} for the same data, and other analyses
}

\examples{
data(Manpower)
mmod <- lm(Hours ~ ., data=Manpower)
vif(mmod)
# ridge regression models, specified in terms of equivalent df
mridge <- ridge(Hours ~ ., data=Manpower, df=seq(5, 3.75, -.25))
vif(mridge)

# univariate ridge trace plots
traceplot(mridge)
traceplot(mridge, X="df")

# bivariate ridge trace plots
plot(mridge, radius=0.25, labels=mridge$df)
pairs(mridge, radius=0.25)

# 3D views
# ellipsoids for Load, Xray & BedDays are nearly 2D
plot3d(mridge, radius=0.2, labels=mridge$df)
# variables in model selected by AIC & BIC
plot3d(mridge, variables=c(2,3,5), radius=0.2, labels=mridge$df)

# plots in PCA/SVD space
mpridge <- pca.ridge(mridge)
traceplot(mpridge, X="df")
biplot(mpridge, radius=0.25)

}
\keyword{datasets}
