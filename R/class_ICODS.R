#' Hidden methods
#'
#' @name ICODS-internal-api
#' @rdname ICODS-internal-api
#' @keywords internal
#' @import methods
tmp <- function(x){}

#' @import methods
.validity_ICODS <- function(object) {

  if ({length(x = object@beta) != length(x = object@se)} ||
      {length(x = object@beta) != length(x = object@pValue)}) {
    return( "beta, se, and pValue are not of same length" )
  }

  if (!any(is.na(object@beta))) {
    if (is.null(x = names(object@beta)) ||
        is.null(x = names(object@se)) ||
        is.null(x = names(object@pValue))) {
      return( "beta, se, and pValue must be named vectors" )
    }
    if (!all(names(object@beta) %in% names(object@se)) ||
        !all(names(object@beta) %in% names(object@pValue))) {
      return( "beta, se, and pValue must have same names" )
    }
  }

  return( TRUE )
}

setClass(Class = "ICODS",
         slots = c("optim" = "list", 
                   "beta" = "vector",
                   "se" = "vector",
                   "pValue" = "vector", 
                   "AIC" = "numeric",
                   "m" = "integer"),
         prototype = list("optim" = list(), 
                          "beta" = NA,  
                          "se" = NA,  
                          "pValue" = NA,  
                          "AIC" = 0L/0L,   
                          "m" = 0L),
         validity = .validity_ICODS)

#' Retrieve the Minimum AIC
#'
#' Retrieves the minimum AIC.
#'   
#' @name minAIC
#' @rdname minAIC
#'
#' @param object An object of class ICODS
#' @param ... ignored
#'
#' @return numeric
#'
#' @export minAIC
#'
#' @examples
#'
#' data(odsData)
#'
#' resultODS <- ODSDesignIC(U = odsData$U, 
#'                          V = odsData$V,  
#'                          del1 = odsData$del1,  
#'                          del2 = odsData$del2, 
#'                          z = odsData$z, 
#'                          mVal = 1L, 
#'                          ind = odsData$ind, 
#'                          a1 = 0.43, 
#'                          a2 = 0.45, 
#'                          beta = NULL, 
#'                          maxit = 10L,
#'                          verbose = TRUE)
#'
#' minAIC(resultODS)
#'
#' data(ccData)
#'
#' resultCC <- CaseCohortIC(U = ccData$U, 
#'                          V = ccData$V,  
#'                          del1 = ccData$del1,  
#'                          del2 = ccData$del2, 
#'                          xi = ccData$xi,
#'                          z = ccData$z, 
#'                          sp = 0.2, 
#'                          mVal = 1L,
#'                          B = 10L, 
#'                          beta = NULL, 
#'                          maxit = 10L,
#'                          verbose = TRUE)
#'
#' minAIC(resultCC)
#'
setGeneric(name = "minAIC",
           def = function(object, ...) { standardGeneric("minAIC") })

#' @rdname ICODS-internal-api
setMethod(f = "minAIC",
          signature = c(object = "ICODS"),
          definition = function(object, ...) {
              return( object@AIC )
            })

#' Retrieve the Optimization Results
#'
#' Retrieves the final optimization results for the m value that minimizes the 
#'   AIC.
#'   
#' @name optimObj
#' @rdname optimObj
#'
#' @param object An object of class ICODS
#' @param ... ignored
#'
#' @return the value object returned by stats::optim()
#' @examples
#'
#' data(odsData)
#'
#' resultODS <- ODSDesignIC(U = odsData$U, 
#'                          V = odsData$V,  
#'                          del1 = odsData$del1,  
#'                          del2 = odsData$del2, 
#'                          z = odsData$z, 
#'                          mVal = 1L, 
#'                          ind = odsData$ind, 
#'                          a1 = 0.43, 
#'                          a2 = 0.45, 
#'                          beta = NULL, 
#'                          maxit = 10L,
#'                          verbose = TRUE)
#'
#' optimObj(resultODS)
#'
#' data(ccData)
#'
#' resultCC <- CaseCohortIC(U = ccData$U, 
#'                          V = ccData$V,  
#'                          del1 = ccData$del1,  
#'                          del2 = ccData$del2, 
#'                          xi = ccData$xi,
#'                          z = ccData$z, 
#'                          sp = 0.2, 
#'                          mVal = 1L,
#'                          B = 10L, 
#'                          beta = NULL, 
#'                          maxit = 10L,
#'                          verbose = TRUE)
#'
#' optimObj(resultCC)
#'
setGeneric(name = "optimObj",
           def = function(object, ...) { standardGeneric("optimObj") })

#' @rdname ICODS-internal-api
#' @export optimObj
setMethod(f = "optimObj",
          signature = c(object = "ICODS"),
          definition = function(object, ...) {
              return( object@optim )
            })

#' Retrieve Degree of Optimal Bernstein Polynomial
#'
#' Retrieves the degree of the Bernstein polynomial basis provided as input
#'   that minimizes the AIC.
#'   
#' @name mVal
#' @rdname mVal
#'
#' @param object An object of class ICODS
#' @param ... ignored
#'
#' @return an integer
#' @examples
#'
#' data(odsData)
#'
#' resultODS <- ODSDesignIC(U = odsData$U, 
#'                          V = odsData$V,  
#'                          del1 = odsData$del1,  
#'                          del2 = odsData$del2, 
#'                          z = odsData$z, 
#'                          mVal = 1L, 
#'                          ind = odsData$ind, 
#'                          a1 = 0.43, 
#'                          a2 = 0.45, 
#'                          beta = NULL, 
#'                          maxit = 10L,
#'                          verbose = TRUE)
#'
#' mVal(resultODS)
#'
#' data(ccData)
#'
#' resultCC <- CaseCohortIC(U = ccData$U, 
#'                          V = ccData$V,  
#'                          del1 = ccData$del1,  
#'                          del2 = ccData$del2, 
#'                          xi = ccData$xi,
#'                          z = ccData$z, 
#'                          sp = 0.2, 
#'                          mVal = 1L,
#'                          B = 10L, 
#'                          beta = NULL, 
#'                          maxit = 10L,
#'                          verbose = TRUE)
#'
#' mVal(resultCC)
#'
setGeneric(name = "mVal",
           def = function(object, ...) { standardGeneric("mVal") })

#' @rdname ICODS-internal-api
#' @export mVal
setMethod(f = "mVal",
          signature = c(object = "ICODS"),
          definition = function(object, ...) {
              return( object@m )
            })

#' Retrieve the Estimated Beta Parameters
#'
#' Retrieves the estimated beta parameters for the m value that minimizes the 
#'   AIC.
#'   
#' @name estimate
#' @rdname estimate
#'
#' @param object An object of class ICODS
#' @param ... ignored
#'
#' @return A matrix containing the estimated parameter value, the standard
#'   error, and the p-value.
#' @examples
#'
#' data(odsData)
#'
#' resultODS <- ODSDesignIC(U = odsData$U, 
#'                          V = odsData$V,  
#'                          del1 = odsData$del1,  
#'                          del2 = odsData$del2, 
#'                          z = odsData$z, 
#'                          mVal = 1L, 
#'                          ind = odsData$ind, 
#'                          a1 = 0.43, 
#'                          a2 = 0.45, 
#'                          beta = NULL, 
#'                          maxit = 10L,
#'                          verbose = TRUE)
#'
#' estimate(resultODS)
#'
#' data(ccData)
#'
#' resultCC <- CaseCohortIC(U = ccData$U, 
#'                          V = ccData$V,  
#'                          del1 = ccData$del1,  
#'                          del2 = ccData$del2, 
#'                          xi = ccData$xi,
#'                          z = ccData$z, 
#'                          sp = 0.2, 
#'                          mVal = 1L,
#'                          B = 10L, 
#'                          beta = NULL, 
#'                          maxit = 10L,
#'                          verbose = TRUE)
#'
#' estimate(resultCC)
#'
setGeneric(name = "estimate",
           def = function(object, ...) { standardGeneric("estimate") })

#' @rdname ICODS-internal-api
#' @export estimate
setMethod(f = "estimate",
          signature = c(object = "ICODS"),
          definition = function(object, ...) {
              mat <- cbind(object@beta, object@se, object@pValue)
              colnames(x = mat) <- c("Estimate", "SE", "p-value")
              rownames(x = mat) <- names(x = object@beta)
              return( mat )
            })

#' @export print
#' @rdname ICODS-internal-api
setMethod(f = "print",
          signature = c(x = "ICODS"),
          definition = function(x, ...) {
              mat <- cbind(x@beta, x@se, x@pValue)
              colnames(x = mat) <- c("Estimate", "SE", "p-value")
              rownames(x = mat) <- names(x = x@beta)

              print(x = mat)
              cat("\n")
              cat("Degree of Bernstein polynomials:", x@m, "\n")
              cat("AIC:", x@AIC, "\n")
            })

#' @export show
#' @rdname ICODS-internal-api
setMethod(f = "show",
          signature = c(object = "ICODS"),
          definition = function(object) {
              mat <- cbind(object@beta, object@se, object@pValue)
              colnames(x = mat) <- c("Estimate", "SE", "p-value")
              rownames(x = mat) <- names(x = object@beta)

              show(object = mat)
              cat("\n")
              cat("Degree of Bernstein polynomials:", object@m, "\n")
              cat("AIC:", object@AIC, "\n")
            })

#' Retrieve the Key Results
#'
#' Retrieves the estimated beta parameters for the m value that minimizes the 
#'   AIC; the m value; and the AIC value.
#'   
#' @name summary
#' @rdname summary
#'
#' @param object An object of class ICODS
#' @param ... ignored
#'
#' @return A list containing
#' \item{par}{A matrix containing the estimated parameter value, the standard
#'   error, and the p-value.}
#' \item{m}{The selected m value.}
#' \item{AIC}{The AIC.}
#'
#' @examples
#'
#' data(odsData)
#'
#' resultODS <- ODSDesignIC(U = odsData$U, 
#'                          V = odsData$V,  
#'                          del1 = odsData$del1,  
#'                          del2 = odsData$del2, 
#'                          z = odsData$z, 
#'                          mVal = 1L, 
#'                          ind = odsData$ind, 
#'                          a1 = 0.43, 
#'                          a2 = 0.45, 
#'                          beta = NULL, 
#'                          maxit = 10L,
#'                          verbose = TRUE)
#'
#' summary(resultODS)
#'
#' data(ccData)
#'
#' resultCC <- CaseCohortIC(U = ccData$U, 
#'                          V = ccData$V,  
#'                          del1 = ccData$del1,  
#'                          del2 = ccData$del2, 
#'                          xi = ccData$xi,
#'                          z = ccData$z, 
#'                          sp = 0.2, 
#'                          mVal = 1L,
#'                          B = 10L, 
#'                          beta = NULL, 
#'                          maxit = 10L,
#'                          verbose = TRUE)
#'
#' summary(resultCC)
#'
NULL

#' @export summary
#' @rdname ICODS-internal-api
setMethod(f = "summary",
          signature = c(object = "ICODS"),
          definition = function(object, ...) {

              mat <- cbind(object@beta, object@se, object@pValue)
              colnames(x = mat) <- c("Estimate", "SE", "p-value")
              rownames(x = mat) <- names(x = object@beta)

              res <- list()

              res$par <- mat
              res$m <- object@m
              res$AIC <- object@AIC

              return( res )
            })
