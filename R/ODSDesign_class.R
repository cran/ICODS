#' @include class_ICODS.R ODSDesign_Obj.R
setClass(Class = "ODSDesign",
         contains = c("ICODS"))

#' @importFrom stats pnorm
.newODSDesign <- function(bernstein_u, 
                          bernstein_v, 
                          z, 
                          del1, 
                          del2, 
                          I1,  
                          I2,  
                          I3,
                          n0,  
                          n1,  
                          n2, 
                          beta, 
                          pis,
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

  # create odsInfo object for compact inputs
  odsInfo <- list("I1" = I1, 
                  "I2" = I2, 
                  "I3" = I3, 
                  "n0" = as.integer(x = n0), 
                  "n1" = as.integer(x = n1), 
                  "n2" = as.integer(x = n2))

  # create initial guess for parameter vector
  if (is.null(x = beta)) {
    initial <- c(rep(x = 0.5, times = np), rep(x = -0.5, times = m), pis)
  } else {
    initial <- c(beta, rep(x = -0.5, times = m), pis)
  }

  res <- .optimStep_ODSDesign(initial = initial, 
                              infoObj = odsInfo, 
                              minData = minData, 
                              m = m, 
                              maxit = maxit,
                              verbose = verbose, ...)

  if (verbose) print(x = res)

  return( res )
}

#' @include CaseCohort_fn.R CaseCohort_gr.R myOptim.R
.optimStep_ODSDesign <- function(initial, 
                                 infoObj,  
                                 minData,  
                                 m,  
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
  argList[[ "fn" ]] <- .ODSDesign_fn
  argList[[ "gr" ]] <- .ODSDesign_gr
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
  mo <- .newMethodObj_ODSDesign(info = infoObj,
                                par = optimResult[[ nopt ]]$par, 
                                minData = minData)

  if (verbose) cat("\nEstimating Standard Error\n")

  se <- .se_ODSDesign(object = mo, argList = argList)

  pValue <- .pValue(object = mo, se = se)

  AIC <- .AIC(object = mo, value = optimResult[[ nopt ]]$value)

  res <- new(Class = "ODSDesign",
             "optim" = optimResult, 
             "beta" = mo$baseInfo$beta,
             "se" = se,
             "pValue" = pValue, 
             "AIC" = AIC,
             "m" = m - 1L)

  return( res )
}
