#' Outcome-Dependent Sampling with Interval-Censored Failure Time Data
#'
#' Provides an outcome-dependent sampling (ODS) design with interval-censored
#'   failure time data, where the observed sample is enriched by selectively 
#'   including certain more informative failure subjects. The method is a 
#'   sieve semiparametric maximum empirical likelihood approach for 
#'   fitting the proportional hazards model to data from the interval-
#'   censoring ODS design. 
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
#'   For left-censored data, the failure occurred prior to the first 
#'   follow-up examination (TF < T1); therefore, define U = T1, V = tau, and 
#'   (del1,del2)=(1,0). For right-censored data, the
#'   failure had not yet occurred at the last follow-up examination 
#'   (TF > TK); therefore, define U = 0, V = TK, 
#'   and (del1,del2)=(0,0). For interval-censored data, the failure occurred
#'   between two follow-up examinations, e.g. T2 < TF < T3; therefore, 
#'   define U and V to be the two consecutive follow-up examination times 
#'   bracketing the failure time TF and (del1,del2)=(0,1). 
#'
#' @references Zhou, Q., Cai, J., and Zhou, H. (2018). Outcome-dependent 
#'   sampling with interval-censored failure time data. Biometrics, 
#'   74(1): 58--67.  <doi:10.1111/biom.12744>
#'
#' @param U numeric vector (n); examination time. 
#'   See Details for further information.
#' @param V numeric vector (n); examination time.
#'   See Details for further information.
#' @param del1 integer vector (n); indicator of a left-censored observation I(T<=U).
#'   See Details for further information.
#' @param del2 integer vector (n); indicator of an interval-censored observation I(U<T<=V).
#'   See Details for further information.
#' @param z matrix (nxp); covariates.
#' @param mVal integer vector (m); one or more options for the degree of 
#'   the Bernstein polynomials. If more than one option provided, the value 
#'   resulting in the lowest
#'   AIC is selected. The results returned are for only that m-value.
#' @param ind integer vector (n); indicating membership of the simple random 
#'   sample (0), lower-tail supplemental sample (1), or upper-tail 
#'   supplemental sample (2).
#' @param a1 numeric (1); lower cut-off point for selecting ODS sample 
#'   (0 < a1 < a2 < tau).
#' @param a2 numeric (1); upper cut-off point for selecting ODS sample 
#'   (0 < a1 < a2 < tau).
#' @param beta numeric vector (p); initial values for beta. If NULL, initial
#'   guess set to 0.5 for each of the p parameters.
#' @param maxit integer(1); maximum number of calls to optimization method.
#' @param verbose logical; TRUE generates progress screen prints.
#' @param ... optional inputs to "control" of function optim().
#'
#' @return an object of class ODSDesign (inheriting from ICODS) containing
#' \item{optim}{a list of the results returned by optim().}
#' \item{beta}{the estimated beta parameters.}
#' \item{se}{the standard error of the estimated beta parameters.}
#' \item{pValue}{the p-value of the estimated beta parameters.}
#' \item{m}{the selected degree of the Bernstein polynomials.}
#' \item{AIC}{the AIC value for the selected degree of the Bernstein polynomials.}
#'
#' @include bernstein.R ODSDesign_class.R
#'
#' @export
#'
#' @examples
#'
#' data(odsData)
#'
#' result <- ODSDesignIC(U = odsData$U, 
#'                       V = odsData$V,  
#'                       del1 = odsData$del1,  
#'                       del2 = odsData$del2, 
#'                       z = odsData$z, 
#'                       mVal = 1L, 
#'                       ind = odsData$ind, 
#'                       a1 = 0.43, 
#'                       a2 = 0.45, 
#'                       beta = NULL, 
#'                       maxit = 10L,
#'                       verbose = TRUE)
#'
#' print(result)
#' mVal(result)
#' estimate(result)
#' optimObj(result)
#' minAIC(result)
#' summary(result)
#'
ODSDesignIC <- function(U, 
                        V, 
                        del1, 
                        del2, 
                        z, 
                        mVal, 
                        ind, 
                        a1, 
                        a2, 
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

  # ensure that 0 < a1 < a2 < tau
  if (a2 < a1) {
    message("reodered a1, a2\n")
    n <- a1
    a1 <- a2
    a2 <- n
  }

  if (a2 > dt$tau) stop("a2 must be less than tau")

  n <- length(x = dt$V)

  # ensure that ind is integer 0,1,2
  if (!is.numeric(x = ind)) stop("ind must be integer 0,1,2")
  if (!is.integer(x = ind)) ind <- as.integer(x = round(x = ind, digits = 0L))
  if (any(!{ind %in% c(0L,1L,2L)})) stop("ind must be integer 0,1,2")
  if (!any(ind == 0L)) stop("there are no members in SRS")
  if (!any(ind == 1L)) stop("there are no members in lower-tail")
  if (!any(ind == 2L)) stop("there are no members in upper-tail")
  if (length(x = ind) != n) stop("ind is not of appropriate length")

  maxInd1 <- max(c(dt$U[ind == 1L & dt$del1 == 1L],
                   dt$V[ind == 1L & dt$del2 == 1L]))

  if (a1 < maxInd1) stop("a1 does not agree with sample data")

  minInd2 <- min(c(dt$U[ind == 2L & dt$del2 == 1L]))

  if (a2 > minInd2) stop("a2 does not agree with sample data")

  # identify SRS
  ind0 <- ind == 0L
  n0 <- sum(ind0)
  n1 <- sum(ind == 1L)
  n2 <- sum(ind == 2L)

  if (verbose) {
    cat("memberships:\n")
    cat("\tSRS:", n0, "\n",
        "\tlower-tail:", n1, "\n",
        "\tupper-tail:", n2, "\n")
  }

  # order data according to membership
  oind <- order(ind)
  Uods <- dt$U[oind]
  Vods <- dt$V[oind]
  del1ods <- dt$del1[oind]
  del2ods <- dt$del2[oind]
  zods <- dt$z[oind,,drop=FALSE]

  # indicators Gk
  I1 <- Uods < a1
  I2 <- Vods < a1
  I3 <- Uods > a2

  n01 <- sum( {{dt$del1[ind0] == 1L} & {dt$U[ind0] <= a1}} | 
              {{dt$del2[ind0] == 1L} & {dt$V[ind0] <= a1}} )
  n02 <- sum( {dt$del2[ind0] == 1L} & {dt$U[ind0] >= a2} )
  n0e <- sum( dt$del1[ind0] + dt$del2[ind0] )

  cutp <- c(n01,n02)/n0e
  pevent <- n0e/n0
  cutp <- cutp / pevent
  pis <- suppressWarnings(-log(1.0/cutp-1.0))

  tst <- is.infinite(pis) | is.nan(pis)

  pis[tst] <- 0.0

  if (verbose) cat("initial pi parameters:", pis, "\n")

  minML <- NULL
  mAIC <- Inf

  for (i in 1L:length(x = dt$mVal)) {

    if (verbose) cat("degree of Bernstein Poly", dt$mVal[i], "\n")

    bernstein_u <- .bernsteinBasis(em = dt$mVal[i], 
                                   sigma = dt$sigma, 
                                   tau = dt$tau, 
                                   tee = Uods)

    bernstein_v <- .bernsteinBasis(em = dt$mVal[i], 
                                   sigma = dt$sigma, 
                                   tau = dt$tau, 
                                   tee = Vods)

    tmp <- .newODSDesign(bernstein_u = bernstein_u,
                         bernstein_v = bernstein_v,
                         z = zods,  
                         del1 = del1ods,  
                         del2 = del2ods,  
                         I1 = I1,
                         I2 = I2,
                         I3 = I3,
                         n0 = n0,
                         n1 = n1,
                         n2 = n2,
                         beta = dt$beta, 
                         pis = pis,
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
