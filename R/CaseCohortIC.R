#' Case-Cohort Studies with Interval-Censored Failure Time Data
#'
#' Provides a sieve semiparametric likelihood approach under the proportional 
#'   hazards model for analyzing data from a case-cohort design with failure 
#'   times subject to interval-censoring. The likelihood function is constructed
#'   using inverse probability weighting, and the sieves are built with 
#'   Bernstein polynomials. A weighted bootstrap procedure is implemented for 
#'   variance estimation. 
#'
#' The implementation uses stats::optim() to minimize the likelihood. The
#'   hard-coded method is "BFGS". Users are able to make changes to the
#'   'control' input of optim() by passing named inputs through the ellipses.
#'   If a call to optim() returns convergence = 1, i.e., optim() reached its
#'   internal maximum number of iterations before convergence was attained,
#'   the software automatically repeats the call to optim() with input 
#'   variable par set to the last parameter values. This procedure is 
#'   repeated at most maxit times.
#'
#' Input parameters U, V, del1, and del2 are defined as follows. 
#'   Suppose there are 
#'   K follow-up examinations at times TE = (T1, T2, ..., TK), and the failure
#'   time is denoted as TF. 
#'   For left-censored data, the failure occurs prior to the first 
#'   follow-up examination (TF < T1); therefore, define U = T1, V = tau, and 
#'   (del1,del2)=(1,0). For right-censored data, the
#'   failure has not yet occurred at the last follow-up examination 
#'   (TF > TK); therefore, define U = 0, V = TK, 
#'   and (del1,del2)=(0,0). For interval-censored data, the failure occurs
#'   between two follow-up examinations, e.g. T2 < TF < T3; therefore, 
#'   define U and V to be the two consecutive follow-up examination times 
#'   bracketing the failure time TF and (del1,del2)=(0,1). 
#'
#' @references Zhou, Q., Zhou, H., and Cai, J. (2017). Case-cohort studies 
#'   with interval-censored failure time data. Biometrika, 104(1): 17--29. 
#'   <doi:10.1093/biomet/asw067>
#'
#' @param U numeric vector (n); examination time. 
#'   See Details for further information.
#' @param V numeric vector (n); examination time.
#'   See Details for further information.
#' @param del1 integer vector (n); indicator of a left-censored observation I(T<=U).
#'   See Details for further information.
#' @param del2 integer vector (n); indicator of an interval-censored observation I(U<T<=V).
#'   See Details for further information.
#' @param xi integer vector (n); indicating membership of the case-cohort sample.
#' @param z matrix (nxp); covariates.
#' @param sp numeric (1); sampling probability 0 < sp < 1.
#' @param mVal integer vector (m); one or more options for the degree of 
#'   the Bernstein polynomials. If more than one option provided, the value 
#'   resulting in the lowest
#'   AIC is selected. The results returned are for only that m-value.
#' @param B integer (1); number of bootstrap samples used to calculate the variance
#'   estimate.
#' @param beta numeric vector (p); initial values for beta. If NULL, initial
#'   guess set to 0.5 for each of the p parameters.
#' @param maxit integer(1); maximum number of calls to optimization method.
#' @param verbose logical; TRUE generates progress screen prints.
#' @param ... optional inputs to "control" of function optim().
#'
#' @return an object of class CaseCohort (inheriting from ICODS) containing
#' \item{optim}{a list of the results returned by optim().}
#' \item{beta}{the estimated beta parameters.}
#' \item{se}{the standard error of the estimated beta parameters.}
#' \item{pValue}{the p-value of the estimated beta parameters.}
#' \item{m}{the selected degree of the Bernstein polynomials.}
#' \item{AIC}{the AIC value for the selected degree of the Bernstein polynomials.}
#'
#' @include bernstein.R CaseCohort_class.R
#'
#' @export
#'
#' @examples
#'
#' data(ccData)
#'
#' result <- CaseCohortIC(U = ccData$U, 
#'                        V = ccData$V,  
#'                        del1 = ccData$del1,  
#'                        del2 = ccData$del2, 
#'                        xi = ccData$xi,
#'                        z = ccData$z, 
#'                        sp = 0.2, 
#'                        mVal = 1L,
#'                        B = 10L, 
#'                        beta = NULL, 
#'                        maxit = 10L,
#'                        verbose = TRUE)
#'
#' print(result)
#' mVal(result)
#' estimate(result)
#' optimObj(result)
#' minAIC(result)
#' summary(result)
#'
CaseCohortIC <- function(U, 
                         V,  
                         del1,  
                         del2,  
                         xi,  
                         z,  
                         sp,  
                         mVal,  
                         B,  
                         beta = NULL, 
                         maxit = 10L, 
                         verbose = TRUE, ...) {

  dt <- .testInputData(U = U,
                       V = V,
                       z = z,
                       del1 = del1,
                       del2 = del2,
                       mVal = mVal,
                       beta = beta,
                       verbose = verbose)

  rm( U, V, z, del1, del2, mVal, beta )

  n <- length(x = dt$V)

  # ensure that xi is integer 0,1
  if (!is.numeric(x = xi)) stop("xi must be integer 0,1")
  if (!is.integer(x = xi)) xi <- as.integer(x = round(x = xi, digits = 0L))
  if (any(!{xi %in% c(0L,1L)})) stop("xi must be integer 0,1")
  if (sum(xi) == 0L) stop("no cases are in cohort")
  if (length(x = xi) != n) stop("xi is not of appropriate length")

  # ensure sampling probability is 0 < sp < 1
  if (sp <= 0.0) stop("sampling probability must be 0 < sp < 1")
  if (sp >= 1.0) stop("sampling probability must be 0 < sp < 1")

  # ensure that number of bootstrap samples is non-negative integer
  if (!is.numeric(x = B)) stop("B must be a positive integer")
  if (!is.integer(x = B)) B <- as.integer(x = round(x = B, digits = 0L))
  if (B <= 0L) stop("B must be a positive integer")

  # probability of observing the covariate Z
  piq <- dt$del1 + dt$del2 + {1L-dt$del1-dt$del2}*sp

  minML <- NULL
  mAIC <- Inf

  for (i in 1L:length(x = dt$mVal)) {

    if (verbose) message(paste("degree of Bernstein Poly", dt$mVal[i]))

    bernstein_u <- .bernsteinBasis(em = dt$mVal[i], 
                                   sigma = dt$sigma, 
                                   tau = dt$tau, 
                                   tee = dt$U)

    bernstein_v <- .bernsteinBasis(em = dt$mVal[i], 
                                   sigma = dt$sigma, 
                                   tau = dt$tau, 
                                   tee = dt$V)

    tmp <- .newCaseCohort(wg = xi / piq,  
                          bernstein_u = bernstein_u,
                          bernstein_v = bernstein_v,
                          z = dt$z,  
                          del1 = dt$del1,  
                          del2 = dt$del2,  
                          B = B,
                          beta = dt$beta, 
                          maxit = maxit,
                          verbose = verbose, ...)

    if (is.null(x = tmp)) {
      warning(paste("unable to obtain estimate for m=", dt$mVal[i]))
      next
    }

    if (minAIC(object = tmp) < mAIC) {
      minML <- tmp
      mAIC <- minAIC(object = tmp)
    }

  }

  if (is.infinite(x = mAIC)) {
    stop("unable to obtain an estimate for any of the specified mVal")
  }

  if (verbose) {
    cat("\nEstimates with minimum AIC\n")
    print(x = minML)
  }

  return( minML )
}
