mfun <- function(r, alpha, nu) {
  2^(1 - nu) / gamma(nu) * (r / alpha)^nu * besselK(r / alpha, nu)
}

mfunspec <- function(r, alpha, nu) {
  4 * pi * alpha^2 * nu * (1 + 4 * pi^2 * alpha^2 * r^2)^{-1-nu}
}

g <- function(r, alpha, nu) {
  1 - mfun(r, alpha, nu)^2
}

g12 <- function(r, tau, alpha12, nu) {
  1 - tau^2 * mfun(r, alpha12, nu)^2
}

alphamax <- function(rho, nu) {
  1 / sqrt(4 * pi * rho * nu)
}

Kspecdiag <- function(r, rho1 = 100, rho2 = 100, alpha1 = 0.03, alpha2 = 0.03, nu1 = 10, nu2 = 10) {
  diag(c(rho1 * mfunspec(r, alpha1, nu1), rho2 * mfunspec(r, alpha2, nu2)))
}

Kspecmatern <- function(r, rho1 = 100, rho2 = 100, alpha1 = 0.03, alpha2 = 0.03, alpha12 = 0.05, tau = 0.2, nu1 = 10, nu2 = 10, nu12 = 10){
  matrix(c(rho1 * mfunspec(r, alpha1, nu1), tau * sqrt(rho1 * rho2) * mfunspec(r, alpha12, nu12),
           tau * sqrt(rho1 * rho2) * mfunspec(r, alpha12, nu12), rho2 * mfunspec(r, alpha2, nu2)), 2, 2)
}
