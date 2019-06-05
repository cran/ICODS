#' @include class_ICODS.R
setClass(Class = "CaseCohort",
         contains = c("ICODS"))


#' @include CaseCohort_Obj.R
#' @importFrom stats pnorm
#' @importFrom stats rexp
#' @importFrom stats sd
.newCaseCohort <- function(wg, 
                           bernstein_u, 
                           bernstein_v, 
                           z, 
                           del1, 
                           del2, 
                           B, 
                           beta, 
                           maxit,
                           verbose, ...) {

  np <- ncol(x = z)
  m <- ncol(x = bernstein_u)

  cbu <- apply(X = bernstein_u,
              MARGIN = 1L,
              FUN = function(x){ rev(x = cumsum(x = rev(x = x)))})

  cbv <- apply(X = bernstein_v,
              MARGIN = 1L,
              FUN = function(x){ rev(x = cumsum(x = rev(x = x)))})

  minData <- list("z" = z, 
                  "bernstein_u" = bernstein_u, 
                  "bernstein_v" = bernstein_v, 
                  "cbu" = cbu, 
                  "cbv" = cbv, 
                  "del1" = del1, 
                  "del2" = del2)

  ccInfo <- list("wg" = wg, 
                 "wb" = 1.0)

  # create initial parameter vector
  if (is.null(x = beta)) {
    initial <- numeric(length = np + m) + 0.5
  } else {
    initial <- c(beta, rep(x = -0.5, times = m))
  }

  res <- .optimStep_CaseCohort(initial = initial, 
                               infoObj = ccInfo,
                               minData = minData,
                               m = m, 
                               B = B, 
                               maxit = maxit,
                               verbose = verbose, ...)

  if (verbose) show(object = res)

  return( res )
}

#' @include CaseCohort_fn.R CaseCohort_gr.R myOptim.R
.optimStep_CaseCohort <- function(initial, 
                                  infoObj,  
                                  minData,  
                                  m,  
                                  B,  
                                  maxit,  
                                  verbose, ...) {

  # generate control input for optim
  controlList <- list(...)

  # if user did not specify a value for trace, base it on the verbose input
  if (is.null(x = controlList[[ "trace" ]])) {
    controlList[[ "trace" ]] <- 1L*verbose
  }

  # verify that user provided inputs are allowed to be modified
  optimControls <- c("trace", 
                     "fnscale",  
                     "parscale",  
                     "ndeps",  
                     "maxit", 
                     "abstol", 
                     "reltol", 
                     "REPORT")

  if (any(!{names(x = controlList) %in% optimControls})) {
    stop("unrecognized inputs for optim() control")
  }

  # create argument list of call to optim()
  argList <- list()

  argList[[ "par" ]] <- initial 
  argList[[ "fn" ]] <- .CaseCohort_fn
  argList[[ "gr" ]] <- .CaseCohort_gr
  argList[[ "method" ]] <- "BFGS"
  argList[[ "control" ]] <- controlList
  argList[[ "minData" ]] <- minData
  argList[[ "info" ]] <- infoObj

  if (verbose) cat("\nOptimizing\n")

  # optimize parameter vectpr
  optimResult <- .myOptim(argList = argList, maxit = maxit)

  # if optimization failed, a NULL value is returned by .myOptim()
  if (is.null(x = optimResult)) return( NULL )

  nopt <- length(x = optimResult)

  # Create an ODSObj with final estimated parameters
  mo <- .newMethodObj_CaseCohort(info = infoObj,
                                 par = optimResult[[ nopt ]]$par, 
                                 minData = minData)

  if (verbose) cat("\nEstimating Standard Error\n")

  se <- .se_CaseCohort(object = mo, B = B, argList = argList)

  pValue <- .pValue(object = mo, se = se)

  AIC <- .AIC(object = mo, value = optimResult[[ nopt ]]$value)

  res <- new(Class = "CaseCohort",
             "optim" = optimResult, 
             "beta" = mo$baseInfo$beta,
             "se" = se,
             "pValue" = pValue, 
             "AIC" = AIC,
             "m" = m - 1L)

  return( res )
}
