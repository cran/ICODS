.testInputData <- function(U, V, z, del1, del2, mVal, beta, verbose) {

  sigma <- min(c(U,V))
  tau <- max(c(U,V))

  if (verbose) cat("sigma =", sigma, "tau =", tau, "\n")

  n <- length(x = V)

  # ensure that del1 is integer 0,1
  if (!is.numeric(x = del1)) stop("del1 must be integer 0,1")
  if (!is.integer(x = del1)) del1 <- as.integer(x = round(x = del1, digits = 0L))
  if (any(!{del1 %in% c(0L,1L)})) stop("del1 must be integer 0,1")

  # ensure that del2 is integer 0,1
  if (!is.numeric(x = del2)) stop("del2 must be integer 0,1")
  if (!is.integer(x = del2)) del2 <- as.integer(x = round(x = del2, digits = 0L))
  if (any(!{del2 %in% c(0L,1L)})) stop("del2 must be integer 0,1")

  if (sum(del1+del2) == 0L) stop("all data are right-censored")
  if (any(!{{del1+del2} %in% c(0L,1L)})) stop("del1+del2 must be integer 0,1")

  # ensure that z is a matrix
  if (!is.matrix(x = z)) z <- as.matrix(x = z)
  if (is.null(x = colnames(x = z))) {
    colnames(z) <- paste0("Z", 1L:ncol(x = z), collapse="")
  }

  # ensure data is all of appropriate length
  if (length(x = U) != n) stop("U and V are not of equivalent length")
  if (length(x = del1) != n) stop("del1 is not of appropriate length")
  if (length(x = del2) != n) stop("del2 is not of appropriate length")
  if (nrow(x = z) != n) stop("z is not of appropriate dimension")

  # ensure ms are non-negative integers
  if (!is.numeric(x = mVal)) stop("mVal must be positive integer(s)")
  if (!is.integer(x = mVal)) mVal <- as.integer(x = round(x = mVal, digits = 0L))
  if (any(mVal <= 0L)) stop("mVal must be positive integer(s)")

  # ensure that if provided, beta is appropriate length
  if (!is.null(x = beta)) {
    if (length(x = beta) != ncol(x = z)) {
      stop("if provided, ncol(z) and length(beta) must be the same")
    }
  }

  return( list("U" = U, 
               "V" = V,  
               "z" = z,  
               "del1" = del1,  
               "del2" = del2,  
               "mVal" = mVal,  
               "beta" = beta,
               "sigma" = sigma,
               "tau" = tau) )

}
