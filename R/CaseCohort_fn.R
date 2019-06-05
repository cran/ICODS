#' @include CaseCohort_Obj.R
.CaseCohort_fn <- function(par, minData, info, ...) {

  mo <- .newMethodObj_CaseCohort(info = info, par = par, minData = minData)

  return( suppressWarnings(.loglik_CaseCohort(object = mo, ...)) )

}
