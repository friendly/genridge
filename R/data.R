
#' Acetylene Data
#' 
#' The data consist of measures of \code{yield} of a chemical manufacturing
#' process for acetylene in relation to numeric parameters.
#' 
#' Marquardt and Snee (1975) used these data to illustrate ridge regression in
#' a model containing quadratic and interaction terms, particularly the need to
#' center and standardize variables appearing in high-order terms.
#' 
#' Typical models for these data include the interaction of \code{temp:ratio},
#' and a squared term in \code{temp}
#' 
#' @name Acetylene
#' @docType data
#' @format A data frame with 16 observations on the following 4 variables.
#' \describe{ 
#'  \item{\code{yield}}{conversion percentage yield of acetylene}
#'  \item{\code{temp}}{reactor temperature (celsius)} 
#'  \item{\code{ratio}}{H2 to N-heptone ratio} 
#'  \item{\code{time}}{contact time (sec)} 
#' }
#' @references 
#' Marquardt, D.W., and Snee, R.D. (1975), "Ridge Regression in
#' Practice," \emph{The American Statistician}, \bold{29}, 3-20.
#' 
#' Marquardt, D.W. (1980), "A Critique of Some Ridge Regression Methods:
#' Comment," \emph{Journal of the American Statistical Association}, Vol. 75,
#' No. 369 (Mar., 1980), pp. 87-91
#' 
#' @source SAS documentation example for \code{PROC REG}, \emph{Ridge
#' Regression for Acetylene Data}.
#' @keywords datasets
#' @examples
#' 
#' data(Acetylene)
#' 
#' # naive model, not using centering
#' amod0 <- lm(yield ~ temp + ratio + time + I(time^2) + temp:time, data=Acetylene)
#' 
#' y <- Acetylene[,"yield"]
#' X0 <- model.matrix(amod0)[,-1]
#' 
#' lambda <- c(0, 0.0005, 0.001, 0.002, 0.005, 0.01)
#' aridge0 <- ridge(y, X0, lambda=lambda)
#' 
#' traceplot(aridge0)
#' traceplot(aridge0, X="df")
#' pairs(aridge0, radius=0.2)
#' 
#' 
#' 
NULL





#' Detroit Homicide Data for 1961-1973
#' 
#' @description 
#' The data set \code{Detroit} was used extensively in the book by Miller
#' (2002) on subset regression. The data are unusual in that a subset of three
#' predictors can be found which gives a very much better fit to the data than
#' the subsets found from the Efroymson stepwise algorithm, or from forward
#' selection or backward elimination. They are also unusual in that, as time
#' series data, the assumption of independence is patently violated, and the
#' data suffer from problems of high collinearity.
#' 
#' As well, ridge regression reveals somewhat paradoxical paths of shrinkage in
#' univariate ridge trace plots, that are more comprehensible in multivariate
#' views.
#' 
#' @details 
#' The data were originally collected and discussed by Fisher (1976) but the
#' complete dataset first appeared in Gunst and Mason (1980, Appendix A).
#' Miller (2002) discusses this dataset throughout his book, but doesn't state
#' clearly which variables he used as predictors and which is the dependent
#' variable.  (\code{Homicide} was the dependent variable, and the predictors
#' were \code{Police} \dots{} \code{WkEarn}.)  The data were obtained from
#' StatLib.
#' 
#' A similar version of this data set, with different variable names appears in
#' the \code{bestglm} package.
#' 
#' @name Detroit
#' @docType data
#' @format A data frame with 13 observations on the following 14 variables.
#' \describe{ 
#'  \item{\code{Police}}{Full-time police per 100,000 population}
#'  \item{\code{Unemp}}{Percent unemployed in the population}
#'  \item{\code{MfgWrk}}{Number of manufacturing workers in thousands}
#'  \item{\code{GunLic}}{Number of handgun licences per 100,000 population}
#'  \item{\code{GunReg}}{Number of handgun registrations per 100,000 population} 
#'  \item{\code{HClear}}{Percent of homicides cleared by arrests}
#'  \item{\code{WhMale}}{Number of white males in the population}
#'  \item{\code{NmfgWrk}}{Number of non-manufacturing workers in thousands}
#'  \item{\code{GovWrk}}{Number of government workers in thousands}
#'  \item{\code{HrEarn}}{Average hourly earnings} 
#'  \item{\code{WkEarn}}{Average weekly earnings} 
#'  \item{\code{Accident}}{Death rate in accidents per 100,000 population} 
#'  \item{\code{Assaults}}{Number of assaults per 100,000 population} 
#'  \item{\code{Homicide}}{Number of homicides per 100,000 of population} 
#' }
#' 
#' @references 
#' Fisher, J.C. (1976). Homicide in Detroit: The Role of Firearms.
#' \emph{Criminology}, \bold{14}, 387--400.
#' 
#' Gunst, R.F. and Mason, R.L. (1980). \emph{Regression analysis and its
#' application: A data-oriented approach}.  Marcel Dekker.
#' 
#' Miller, A. J. (2002). \emph{Subset Selection in Regression}. 2nd Ed. Chapman
#' & Hall/CRC. Boca Raton.
#' 
#' @source \url{https://lib.stat.cmu.edu/datasets/detroit}
#' @keywords datasets
#' @examples
#' 
#' data(Detroit)
#' 
#' # Work with a subset of predictors, from Miller (2002, Table 3.14),
#' # the "best" 6 variable model
#' #    Variables: Police, Unemp, GunLic, HClear, WhMale, WkEarn
#' # Scale these for comparison with other methods
#' 
#' Det <- as.data.frame(scale(Detroit[,c(1,2,4,6,7,11)]))
#' Det <- cbind(Det, Homicide=Detroit[,"Homicide"])
#' 
#' # use the formula interface; specify ridge constants in terms
#' # of equivalent degrees of freedom
#' dridge <- ridge(Homicide ~ ., data=Det, df=seq(6,4,-.5))
#' 
#' # univariate trace plots are seemingly paradoxical in that
#' # some coefficients "shrink" *away* from 0
#' traceplot(dridge, X="df")
#' vif(dridge)
#' pairs(dridge, radius=0.5)
#' 
#' \donttest{
#' plot3d(dridge, radius=0.5, labels=dridge$df)
#' 
#' # transform to PCA/SVD space
#' dpridge <- pca(dridge)
#'
#' # not so paradoxical in PCA space
#' traceplot(dpridge, X="df")
#' biplot(dpridge, radius=0.5, labels=dpridge$df)
#' 
#' # show PCA vectors in variable space
#' biplot(dridge, radius=0.5, labels=dridge$df)
#' }
#' 
#' 
NULL


#' Hospital manpower data
#' 
#' @description 
#' The hospital manpower data, taken from Myers (1990), table 3.8, are a
#' well-known example of highly collinear data to which ridge regression and
#' various shrinkage and selection methods are often applied.
#' 
#' The data consist of measures taken at 17 U.S. Naval Hospitals and the goal
#' is to predict the required monthly man hours for staffing purposes.
#' 
#' @details
#' Myers (1990) indicates his source was "Procedures and Analysis for Staffing
#' Standards Development: Data/Regression Analysis Handbook", Navy Manpower and
#' Material Analysis Center, San Diego, 1979.
#' 
#' @name Manpower
#' @docType data
#' @format A data frame with 17 observations on the following 6 variables.
#' \describe{ 
#'  \item{\code{Hours}}{monthly man hours (response variable)}
#'  \item{\code{Load}}{average daily patient load} 
#'  \item{\code{Xray}}{monthly X-ray exposures} 
#'  \item{\code{BedDays}}{monthly occupied bed days}
#'  \item{\code{AreaPop}}{eligible population in the area in thousands}
#'  \item{\code{Stay}}{average length of patient's stay in days} 
#' }
#' @seealso \code{\link[bestglm]{manpower}} for the same data, and other
#' analyses
#' 
#' @references 
#' Donald R. Jensen and Donald E. Ramirez (2012). Variations on
#' Ridge Traces in Regression, \emph{Communications in Statistics - Simulation
#' and Computation}, 41 (2), 265-278.
#' 
#' @source 
#' Raymond H. Myers (1990). \emph{Classical and Modern Regression with
#' Applications}, 2nd ed., PWS-Kent, pp. 130-133.
#' 
#' @keywords datasets
#' @examples
#' 
#' data(Manpower)
#' mmod <- lm(Hours ~ ., data=Manpower)
#' vif(mmod)
#' # ridge regression models, specified in terms of equivalent df
#' mridge <- ridge(Hours ~ ., data=Manpower, df=seq(5, 3.75, -.25))
#' vif(mridge)
#' 
#' # univariate ridge trace plots
#' traceplot(mridge)
#' traceplot(mridge, X="df")
#' 
#' # bivariate ridge trace plots
#' plot(mridge, radius=0.25, labels=mridge$df)
#' pairs(mridge, radius=0.25)
#' 
#' \donttest{
#' # 3D views
#' # ellipsoids for Load, Xray & BedDays are nearly 2D
#' plot3d(mridge, radius=0.2, labels=mridge$df)
#' # variables in model selected by AIC & BIC
#' plot3d(mridge, variables=c(2,3,5), radius=0.2, labels=mridge$df)
#' 
#' # plots in PCA/SVD space
#' mpridge <- pca(mridge)
#' traceplot(mpridge, X="df")
#' biplot(mpridge, radius=0.25)
#' }
#' 
#' 
NULL


#' Prostate Cancer Data
#' 
#' @description
#' Data to examine the correlation between the level of prostate-specific
#' antigen and a number of clinical measures in men who were about to receive a
#' radical prostatectomy.
#' 
#' @details
#' This data set came originally from the (now defunct) ElemStatLearn package.
#' 
#' The last column indicates which 67 observations were used as the "training
#' set" and which 30 as the test set, as described on page 48 in the book.
#' 
#' @name prostate
#' @docType data
#' @format A data frame with 97 observations on the following 10 variables.
#' \describe{ 
#'  \item{lcavol}{log cancer volume} 
#'  \item{lweight}{log prostate weight} 
#'  \item{age}{in years} 
#'  \item{lbph}{log of the amount of benign prostatic hyperplasia} 
#'  \item{svi}{seminal vesicle invasion} 
#'  \item{lcp}{log of capsular penetration} 
#'  \item{gleason}{a numeric vector}
#'  \item{pgg45}{percent of Gleason score 4 or 5} 
#'  \item{lpsa}{response}
#'  \item{train}{a logical vector} 
#' }
#' @note There was an error in this dataset in earlier versions of the package,
#' as indicated in a footnote on page 3 of the second edition of the book. As
#' of version 2012.04-0 this was corrected.
#' 
#' @source 
#' Stamey, T., Kabalin, J., McNeal, J., Johnstone, I., Freiha, F.,
#' Redwine, E. and Yang, N (1989) Prostate specific antigen in the diagnosis
#' and treatment of adenocarcinoma of the prostate II. Radical prostatectomy
#' treated patients, \emph{Journal of Urology}, \bold{16}: 1076--1083.
#' @keywords datasets
#' @examples
#' 
#' data(prostate)
#' str( prostate )
#' cor( prostate[,1:8] )
#' prostate <- prostate[, -10]
#' 
#' prostate.mod <- lm(lpsa ~ ., data=prostate)
#' vif(prostate.mod)
#' 
#' py <- prostate[, "lpsa"]
#' pX <- data.matrix(prostate[, 1:8])
#' pridge <- ridge(py, pX, df=8:1)
#' pridge
#' 
#' # univariate ridge trace plots
#' traceplot(pridge)
#' traceplot(pridge, X="df")
#' 
#' # bivariate ridge trace plots
#' plot(pridge)
#' pairs(pridge)
#' 
#' 
#' 
NULL



