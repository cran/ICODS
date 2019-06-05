#' @include minObj.R
.newBaseInfo <- function(par, minData) {

  np <- ncol(x = minData$z)
  if (is.null(x = np)) np <- 1L
  m <- ncol(x = minData$bernstein_u)

  if (length(x = par) != {np + m}) stop("dim error")

  beta <- par[1L:np]
  names(x = beta) <- colnames(x = minData$z)

  et <- exp(x = par[{np+1L}:{np+m}])
  eta <- cumsum(x = et)

  return( list("U" = .newMinObj(beta = beta, 
                                eta = eta,  
                                bernstein = minData$bernstein_u,  
                                cb = minData$cbu,
                                z = minData$z),  
               "V" = .newMinObj(beta = beta, 
                                eta = eta,  
                                bernstein = minData$bernstein_v,  
                                cb = minData$cbv,
                                z = minData$z),  
               "beta" = beta,  
               "et" = unname(obj = et),
               "del1" = minData$del1,
               "del2" = minData$del2) )

}
