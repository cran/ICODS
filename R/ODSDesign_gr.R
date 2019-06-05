#' @include ODSDesign_Obj.R
.ODSDesign_gr <- function(par, minData, info, ...) {

  mo <- .newMethodObj_ODSDesign(info = info, par = par, minData = minData)

  return( suppressWarnings(colSums(x = .dloglik_ODSDesign(object = mo, ...))) )

}
