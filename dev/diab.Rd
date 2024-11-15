#' @name diab
#' @aliases diab
#' @docType data
#' @title
#' Diabetes Progression
#'
#' @description
#' These data consist of observations on 442 patients, with the response of interest being 
#' a quantitative measure of disease progression one year after baseline.
#'
#' There are ten baseline variables: age, sex, body-mass index (\code{bmi}), average blood pressure (\code{map})
#' and six blood serum measurements.
#'
#' @usage data("diab")
#' @format
#'  A data frame with 442 observations on the following 11 variables.
#'  \describe{
#'    \item{\code{prog}}{integer, a numeric vector}
#'    \item{\code{age}}{integer, a numeric vector}
#'    \item{\code{sex}}{integer, a numeric vector}
#'    \item{\code{bmi}}{numeric, a numeric vector}
#'    \item{\code{map}}{numeric, a numeric vector}
#'    \item{\code{tc}}{integer, a numeric vector}
#'    \item{\code{ldl}}{numeric, a numeric vector}
#'    \item{\code{hdl}}{numeric, a numeric vector}
#'    \item{\code{tch}}{numeric, a numeric vector}
#'    \item{\code{ltg}}{numeric, a numeric vector}
#'    \item{\code{glu}}{integer, a numeric vector}
#'  }
#'
#' @details 
#' %%  ~~ If necessary, more details than the __description__ above ~~
#' @source 
#' The dataset was taken from the web site for Efron & Hastie (2021), 
#' \url{http://hastie.su.domains/CASI_files/DATA/diabetes.csv}.
#' @references
#' %%  ~~ possibly secondary sources and usages ~~
#'
#' @concept %% non-keyword concepts, one per line
#' @examples
#' data(diab)
#' ## maybe str(diab) ; plot(diab) ...
#'
#' @keywords datasets
NULL
