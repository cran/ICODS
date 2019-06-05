# Calculates the b_{k,en}(x) Bernstein polynomial
#
# function is not for export
.bernstein <- function(k, en, x) {

  if (any(x < 0.0 | x > 1.0)) stop("x must be [0,1] for Bernstein poly")
  if (k > en) stop("en must be greater than k")

  return( choose(n = en, k = k) * {x^k} * {{1.0 - x}^{en - k}} )

}

# Calculates the em+1 Bernstein basis polynomials of degree em
# @param em integer; degree of Bernstein polynomials
# @param sigma numeric; lower bound used to ensure that 0 <= x <= 1
# @param tau numeric; upper bound used to ensure that 0 <= x <= 1
# @param tee numeric; unscaled value of x
.bernsteinBasis <- function(em, sigma, tau, tee) {

    x <- {tee-sigma} / {tau-sigma}

    return( sapply(X = 0L:em, FUN = .bernstein, en = em, x = x) )

}
