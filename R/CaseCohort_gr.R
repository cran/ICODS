#' @include CaseCohort_Obj.R
.CaseCohort_gr <- function(par, minData, info, ...) {

  mo <- .newMethodObj_CaseCohort(info = info, par = par, minData = minData)

  return( suppressWarnings(colSums(x = .dloglik_CaseCohort(object = mo, ...))) )

}
