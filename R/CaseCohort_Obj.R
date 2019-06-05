#' @include baseInfo.R
.newMethodObj_CaseCohort <- function(info, par, minData, ...) {

  base <- .newBaseInfo(par = par, minData = minData)

  return( list("wg" = info$wg,
               "wb" = info$wb,
               "np" = length(x = base$beta) + length(x = base$et),
               "baseInfo" = base) )
 
}

.loglik_CaseCohort <- function(object, ...) {

  Su <- object$baseInfo$U$S
  Sv <- object$baseInfo$V$S

  res <- -sum(object$wg*object$wb*
              {object$baseInfo$del1 * log(x = 1.0-Su) + 
               object$baseInfo$del2 * log(x = Su-Sv) + 
               {1.0-object$baseInfo$del1-object$baseInfo$del2} * log(x = Sv)})

  if (is.nan(x = res)) return( Inf )

  return( res )
}

.dloglik_CaseCohort <- function(object, ...) {

  Su <- object$baseInfo$U$S
  Sv <- object$baseInfo$V$S

  dSu <- .deriv1S(object = object$baseInfo$U, 
                  et = object$baseInfo$et, 
                  beta = object$baseInfo$beta)
  dSv <- .deriv1S(object = object$baseInfo$V, 
                  et = object$baseInfo$et, 
                  beta = object$baseInfo$beta)

  temp11 <- -dSu / {1.0-Su}
  temp12 <- {dSu - dSv} / {Su-Sv}
  temp13 <- dSv / Sv

  res <- object$baseInfo$del1*temp11 + 
         object$baseInfo$del2*temp12 + 
         {1.0-object$baseInfo$del1-object$baseInfo$del2}*temp13

  return( unname(-res*object$wb*object$wg) )

}

.ddloglik_CaseCohort <- function(object, ...) {

  n <- length(x = object$baseInfo$del1)
  np <- object$np

  res <- matrix(data = 0.0, nrow = np, ncol = np)

  Su <- object$baseInfo$U$S
  Sv <- object$baseInfo$V$S

  if (length(x = object$wb) == 1L) {
    object$wb <- rep(x = object$wb, times = n)
  }
  if (length(x = object$wg) == 1L) {
    object$wg <- rep(x = object$wg, times = n)
  }

  for (i in 1L:n) {

    dSu <- .derivS(object = object$baseInfo$U, 
                   i = i, 
                   et = object$baseInfo$et, 
                   beta = object$baseInfo$beta)

    dSv <- .derivS(object = object$baseInfo$V, 
                   i = i, 
                   et = object$baseInfo$et, 
                   beta = object$baseInfo$beta)

    temp11 <- -dSu$Stt / {1.0-Su[i]} - dSu$St %o% dSu$St / {{1.0-Su[i]}^2}
    temp12 <- {dSu$Stt - dSv$Stt} / {Su[i]-Sv[i]} -
              {dSu$St - dSv$St} %o% {dSu$St - dSv$St} / {{Su[i]-Sv[i]}^2}
    temp13 <- dSv$Stt / Sv[i] - dSv$St %o% dSv$St/{Sv[i]^2}

    res <- res + {object$baseInfo$del1[i]*temp11 + 
                  object$baseInfo$del2[i]*temp12 + 
                 {1.0-object$baseInfo$del1[i]-object$baseInfo$del2[i]}*temp13}*
                  object$wb[i]*object$wg[i]
  }

  return( unname(-res) )

}

#' @include myOptim.R
.se_CaseCohort <- function(object, B, argList, ...) {

  np <- length(x = object$baseInfo$beta)
  n <- length(x = object$baseInfo$del1)

  # estimate standard error via bootstrap
  boot <- matrix(data = NA, nrow = B, ncol = np)

  for (b in 1L:B) {

    argList[[ "info" ]] <- list("wg" = object$wg, 
                                "wb" = rexp(n = n, rate = 1.0))

    tmp <- .myOptim(argList = argList)

    if (is.null(x = tmp)) {
      warning(paste("optim did not converge for bootstrap iteration", b))
    } else {
      i <- length(x = tmp)
      boot[b,] <- tmp[[ i ]]$par[1L:np]
    }
  }

  # standard error across bootstrap samples for each parameter
  se <- drop(x = apply(X = boot, 
                       MARGIN = 2L, 
                       FUN = sd, 
                       na.rm = TRUE))
  names(x = se) <- names(x = object$baseInfo$beta)

  return( se )
}

#' @importFrom stats pnorm
.pValue <- function(object, se, ...) {
  pValue <- 2.0*{1.0 - pnorm(q = abs(x = object$baseInfo$beta / se), 
                             mean = 0.0, 
                             sd = 1.0)}
  names(x = pValue) <- names(x = object$baseInfo$beta)
  return( pValue )
}

.AIC <- function(object, value, ...) {
  n <- length(object$baseInfo$beta) + length(object$baseInfo$et)
  return( 2.0*{value + n} )
}
