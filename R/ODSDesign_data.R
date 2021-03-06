#' Toy Example for ODS Design with Interval-Censored Data
#' 
#' This data set gives a simple toy example of ODS design with interval-censored data.
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
#' @usage data(odsData)
#'
#' @format A data.frame containing 501 observations with 6 columns:
#' \describe{
#' \item{U}{examination time; see Details.}
#' \item{V}{examination time; see Details.}
#' \item{del1}{indicator of a left-censored observation I(T<=U).} 
#' \item{del2}{indicator of an interval-censored observation I(U<T<=V).}
#' \item{z}{covariates.}
#' \item{ind}{indicating membership of the simple random sample (0), lower-tail 
#'   supplemental sample (1), or upper-tail supplemental sample (2).}
#' }
#'
#' @references Zhou, Q., Cai, J., and Zhou, H. (2018). Outcome-dependent 
#'   sampling with interval-censored failure time data. Biometrics, 
#'   74(1): 58--67.  <doi:10.1111/biom.12744>
#'
#' @keywords datasets
"odsData"
