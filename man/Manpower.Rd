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
    \item{\code{Load}}{average daily patient load}
    \item{\code{Xray}}{monthly X-ray exposures}
    \item{\code{BedDays}}{monthly occupied bed days}
    \item{\code{AreaPop}}{eligible population in the area in thousands}
    \item{\code{Stay}}{average length of patient's stay in days}
    \item{\code{Hours}}{monthly man hours (response variable)}
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
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(Manpower)
## maybe str(Manpower) ; plot(Manpower) ...
}
\keyword{datasets}
