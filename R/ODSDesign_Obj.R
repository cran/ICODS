.newMethodObj_ODSDesign <- function(info, par, minData) {

  npar <- length(x = par)

  ccObj <- .newMethodObj_CaseCohort(info = list("wg" = 1.0, "wb" = 1.0), 
                                    par = par[1L:{npar-2L}],  
                                    minData = minData)


  return( c(info, ccObj, 
            "pi1" = 1.0 / {1.0 + exp(x = -par[npar-1L])}, 
            "pi2" = 1.0 / {1.0 + exp(x = -par[npar])}) )

}

.G <- function(object, ...) {

  Su <- object$baseInfo$U$S
  Sv <- object$baseInfo$V$S

  G1 <- object$I1 * {1.0-Su} + object$I2 * {Su - Sv}
  G2 <- object$I3 * {Su - Sv}

  return( list("G1" = G1, "G2" = G2) )
}

# first and second derivatives wrt non-pi parameters
.dG <- function(object, i, ...) {

  dSu <- .derivS(object = object$baseInfo$U, 
                 i = i, 
                 beta = object$baseInfo$beta, 
                 et = object$baseInfo$et)

  dSv <- .derivS(object = object$baseInfo$V, 
                 i = i, 
                 beta = object$baseInfo$beta, 
                 et = object$baseInfo$et)

  G1t <- -object$I1[i]*dSu$St + object$I2[i]*{dSu$St-dSv$St}
  G2t <- object$I3[i]*{dSu$St - dSv$St}

  G1tt <- -object$I1[i]*dSu$Stt + object$I2[i]*{dSu$Stt - dSv$Stt}
  G2tt <- object$I3[i]*{dSu$Stt - dSv$Stt}
    
  return( list("d1g1" = unname(obj = G1t), 
               "d1g2" = unname(obj = G2t), 
               "d2g1" = unname(obj = G1tt), 
               "d2g2" = unname(obj = G2tt)) )
}

# first derivative wrt non-pi parameters
.d1G <- function(object, ...) {

  dSu <- .deriv1S(object = object$baseInfo$U, 
                  beta = object$baseInfo$beta, 
                  et = object$baseInfo$et)

  dSv <- .deriv1S(object = object$baseInfo$V, 
                  beta = object$baseInfo$beta, 
                  et = object$baseInfo$et)

  G1t <- -object$I1*dSu + object$I2*{dSu-dSv}
  G2t <- object$I3*{dSu - dSv}

  return( list("d1g1" = unname(obj = G1t), "d1g2" = unname(obj = G2t)) )
}

.piece <- function(object) {

  Gs <- .G(object = object)

  slz <- sum(log(x = {object$n0 + 
                      object$n1*Gs$G1/object$pi1 + 
                      object$n2*Gs$G2/object$pi2}))

  res <- slz + object$n1*log(x = object$pi1) + object$n2*log(x = object$pi2)

  return( res )

}

.dPiece <- function(object) {

  Gs <- .G(object = object)

  n1p1 <- object$n1 / object$pi1
  n2p2 <- object$n2 / object$pi2

  temp2 <- -n1p1*Gs$G1*{1.0-object$pi1}
  temp3 <- -n2p2*Gs$G2*{1.0-object$pi2}
  temp4 <- 1.0/{object$n0 + n1p1*Gs$G1 + n2p2*Gs$G2}

  temp5 <- c(rep(x = 0, times = object$np), 
             object$n1*{1.0-object$pi1},  
             object$n2*{1.0-object$pi2})

  dG1 <- .d1G(object = object)

  temp1 <- n1p1*dG1$d1g1 + n2p2*dG1$d1g2

  n <- nrow(x = temp1)

  l2xi <- cbind(temp1*temp4, temp2*temp4, temp3*temp4)

  return( unname(obj = t(x = {t(x = l2xi) + temp5/n})) )
}

.ddPiece <- function(object) {

  np <- object$np + 2L

  Gs <- .G(object = object)

  n1p1 <- object$n1 / object$pi1
  n2p2 <- object$n2 / object$pi2

  term <- object$n0 + n1p1*Gs$G1 + n2p2*Gs$G2

  # derivative wrt pi1/pi2 parameters
  hess2 <- matrix(data = 0.0, nrow = 2L, ncol = 2L)

  temp33 <- n1p1*Gs$G1*{1.0-object$pi1}/term
  temp34 <- n2p2*Gs$G2*{1.0-object$pi2}/term

  hess2[1L,1L] <- sum(temp33*{temp33 - 1.0})
  hess2[1L,2L] <- sum(temp33*temp34)
  hess2[2L,2L] <- sum(temp34*{temp34 - 1.0})

  # d/dpi d/dbeta, d/dpi d/deta terms
  res <- matrix(data = 0.0, nrow = {np-2L}, ncol = {np-2L})
  hess4 <- matrix(data = 0.0, nrow = 2L, ncol = {np-2L})

  for (i in 1L:length(x = term)) {

    dGs <- .dG(object = object, i = i)

    temp21 <- n1p1*dGs$d2g1 + n2p2*dGs$d2g2
    temp22 <- n1p1*dGs$d1g1 + n2p2*dGs$d1g2
    temp22 <- temp22 %o% temp22

    itemp23 <- 1.0/term[i]

    res <- res + {-temp21 + temp22*itemp23}*itemp23

    temp31 <- n1p1*dGs$d1g1*{1.0-object$pi1}*itemp23
    temp32 <- n2p2*dGs$d1g2*{1.0-object$pi2}*itemp23
    temp35 <- {n1p1*dGs$d1g1 + n2p2*dGs$d1g2}*itemp23

    hess4[1L,] <- hess4[1L,] + {temp31 - temp35*temp33[i]}
    hess4[2L,] <- hess4[2L,] + {temp32 - temp35*temp34[i]}
  }

  hess1 <- res
  hess3 <- t(x = hess4)
  hess2[1L,1L] <- hess2[1L,1L] - object$n1*object$pi1*(object$pi1-1L)
  hess2[2L,2L] <- hess2[2L,2L] - object$n2*object$pi2*(object$pi2-1L)
  hess2[2L,1L] <- hess2[1L,2L]

  return( -rbind(cbind(hess1, hess3),cbind(hess4,hess2)) )

}

.loglik_ODSDesign <- function(object, ...) {

  if (object$pi1 < 1e-8 || object$pi1 >= 1.0) return( Inf )
  if (object$pi2 < 1e-8 || object$pi2 >= 1.0) return( Inf )

  # calculate the weighted log likelihood
  logl <- - .loglik_CaseCohort(object = object, ...)

  res <- -{logl - .piece(object = object)}

  if (is.nan(x = res)) return( Inf )

  return( res )

}

.dloglik_ODSDesign <- function(object, ...) {

  # calculate the weighted log likelihood
  lxi <- .dloglik_CaseCohort(object = object, ...)

  # add columns for pi1 and pi2 parameters
  lxi <- cbind(-lxi, 0.0, 0.0)

  l2 <- .dPiece(object = object)

  lxi <- lxi - l2

  return( -lxi )

}

.ddloglik_ODSDesign <- function(object, ...) {

  hess1 <- .ddloglik_CaseCohort(object = object, ...)

  hess1 <- rbind(cbind(hess1, 0, 0),0,0)

  hess2 <- .ddPiece(object = object)

  return( unname(obj = hess1 + hess2) )

}

#' @importFrom MASS ginv
#' @importFrom stats cov
.se_ODSDesign <- function(object, ...) {

  np <- length(x = object$baseInfo$beta)

  # variance
  Jh <- .ddloglik_ODSDesign(object = object)
  res <- .dloglik_ODSDesign(object = object)

  V0 <- cov(x = res[1L:object$n0,,drop=FALSE])
  V1 <- cov(x = res[{object$n0+1L}:{object$n0+object$n1},,drop=FALSE])
  V2 <- cov(x = res[{object$n0+object$n1+1L}:{object$n0+object$n1+object$n2},,drop=FALSE])
  
  Sigh <- object$n0*V0 + object$n1*V1 + object$n2*V2

  Sigma <- ginv(X = Jh) %*% Sigh %*% ginv(X = Jh)

  se <- sqrt(x = diag(Sigma[1L:np,1L:np,drop=FALSE]))
  names(x = se) <- names(x = object$baseInfo$beta)

  return( se )

}
