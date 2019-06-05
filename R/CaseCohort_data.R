#' Toy Example for Case-Cohort Design with Interval-Censored Data
#' 
#' This data set gives a simple toy example of case-cohort design with interval-censored data.
#' It was generated following the simulation study used to illustrate the method
#' in the original manuscript referenced below. 
#' This dataset is not meaningful and is intended for demonstration purposes only.
#' 
#' The data can be understood as follow.
#'   There are 
#'   K follow-up examinations at times TE = (T1, T2, ..., TK), and the failure
#'   time is denoted as TF. 
#'   For left-censored data, the failure occurred prior to the first 
#'   follow-up examination (TF < T1); therefore, U = T1, V = tau, and 
#'   (del1,del2)=(1,0). For right-censored data, the
#'   failure had not yet occurred at the last follow-up examination 
#'   (TF > TK); therefore, U = 0, V = TK, 
#'   and (del1,del2)=(0,0). For interval-censored data, the failure occurred
#'   between two follow-up examinations, e.g. T2 < TF < T3; therefore, 
#'   U and V to be the two consecutive follow-up examination times 
#'   bracketing the failure time TF and (del1,del2)=(0,1). 
#'
#' @usage data(ccData)
#'
#' @format A data.frame containing 500 observations with 6 columns: 
#' \describe{
#' \item{U}{examination time.}
#' \item{V}{examination time.}
#' \item{del1}{indicator of a left-censored observation I(T<=U).} 
#' \item{del2}{indicator of an interval-censored observation I(U<T<=V).} 
#' \item{xi}{indicating membership of the case-cohort sample.}
#' \item{z}{covariates.}
#' }
#' See Details for further information.
#'
#' @references Zhou, Q., Zhou, H., and Cai, J. (2017). Case-cohort studies 
#'   with interval-censored failure time data. Biometrika, 104(1): 17--29. 
#'   <doi:10.1093/biomet/asw067>
#'
#' @keywords datasets
"ccData"
