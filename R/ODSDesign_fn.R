#' @include ODSDesign_Obj.R
.ODSDesign_fn <- function(par, minData, info, ...) {

  mo <- .newMethodObj_ODSDesign(info = info, par = par, minData = minData)

  return( suppressWarnings(.loglik_ODSDesign(object = mo, ...)) )

}
