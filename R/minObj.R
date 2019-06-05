.newMinObj <- function(beta, eta, cb, bernstein, z) {

  if (!is.matrix(x = z)) z <- matrix(data = z, 
                                     ncol = 1L,  
                                     dimnames = list(NULL, "z1"))

  L <- drop(x = bernstein %*% eta)

  S <- drop(x = exp(x = -L * drop(x = exp(x = z %*% beta))))

  return( list("S" = S, 
               "bernstein" = bernstein, 
               "cb" = cb,  
               "L" = L,  
               "z" = z) )

}

.deriv1S <- function(object, et, beta, ...) {

  ezb <- drop(x = exp(object$z %*% beta))

  cbet <- t(x = object$cb*et)

  Sb <- {-object$S*object$L*ezb}*object$z #{n x p}
  Se <- {-object$S*ezb}*cbet #{n x m}
  St <- cbind(Sb, Se) #{p+m}

  return( St )
}


.derivS <- function(object, i, et, beta, ...) {

  ezb <- drop(x = exp(object$z[i,] %*% beta))

  cbet <- object$cb[,i]*et

  Sb <- {-object$S[i]*object$L[i]*ezb}*object$z[i,] #{p}
  Se <- {-object$S[i]*ezb}*cbet #{m}
  St <- c(Sb, Se) #{p+m}

  Sbb <- {{-object$L[i]*ezb}*object$z[i,]} %o% {Sb + object$S[i]*object$z[i,]} #{pxp}

  See <- {-ezb*cbet} %o% Se #{mxm}
              # verify why this is different
  See <- -Se %o% {ezb*cbet} + diag(Se)

  Sbe <- {-ezb*object$z[i,]} %o% {Se*object$L[i] + object$S[i]*cbet} #{pxm}
  Stt <- rbind(cbind(Sbb,Sbe),
               cbind(t(x = Sbe),See)) #{p+m x p+m} 

  return( list("St" = St, "Stt" = Stt) )
}
