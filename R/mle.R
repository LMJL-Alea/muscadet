#' @export
mle_dpp <- function(X, labels, lb, ub, rho1, rho2, alpha1, alpha2) {
  # optim(par = 0.5, fn = Evaluate, method = "L-BFGS-B", lower = 1e-4, upper = 1-1e-4,
  #       X = X, labels = labels, lb = lb, ub = ub, rho1 = rho1, rho2 = rho2, alpha1 = alpha1, alpha2 = alpha2)
  # optim(par = 0, fn = Evaluate, method = "L-BFGS-B", lower = -0.28, upper = 0.28,
  #       X = X, labels = labels, lb = lb, ub = ub, rho1 = rho1, rho2 = rho2, alpha1 = alpha1, alpha2 = alpha2)
  # optim(par = 0.02, fn = Evaluate, method = "L-BFGS-B", lower = 1e-4, upper = 0.053,
  #       X = X, labels = labels, lb = lb, ub = ub, rho1 = rho1, rho2 = rho2, alpha1 = alpha1, alpha2 = alpha2)
  optim(par = 0.035, fn = Evaluate, method = "L-BFGS-B", lower = 0.03, upper = 0.045,
        X = X, labels = labels, lb = lb, ub = ub, rho1 = rho1, rho2 = rho2, alpha1 = alpha1, alpha2 = alpha2)
}

#' @export
mle_dpp_bessel <- function(X, labels, lb, ub, a1 = 0, a2 = 0, a12 = 0, alpha1 = 0, alpha2 = 0, alpha12 = 0) {
  d <- ncol(X)
  # # rho_i
  # optimise(f = EvaluateBessel, interval = c(1e-4, 1-1e-4),
  #          X = X, labels = labels, lb = lb, ub = ub,
  #          amplitude1 = a1, amplitude2 = a2, amplitude12 = a12,
  #          alpha1 = alpha1, alpha2 = alpha2, alpha12 = alpha12)
  # # alpha_i
  # optimise(f = EvaluateBessel, interval = c(1e-4, 0.053),
  #          X = X, labels = labels, lb = lb, ub = ub,
  #          amplitude1 = a1, amplitude2 = a2, amplitude12 = a12,
  #          alpha1 = alpha1, alpha2 = alpha2, alpha12 = alpha12)
  # rho12
  optimise(f = EvaluateBessel, interval = c(-0.28, 0.28),
           X = X, labels = labels, lb = lb, ub = ub,
           amplitude1 = a1, amplitude2 = a2, amplitude12 = a12,
           alpha1 = alpha1, alpha2 = alpha2, alpha12 = alpha12)
}

#' @export
test_eval <- function(p, X, labels, lb, ub, a1 = 0, a2 = 0, a12 = 0, alpha1 = 0, alpha2 = 0, alpha12 = 0) {
  d <- ncol(X)
  EvaluateBessel(p,
        X = X, labels = labels, lb = lb, ub = ub,
        amplitude1 = a1, amplitude2 = a2, amplitude12 = a12,
        alpha1 = alpha1, alpha2 = alpha2, alpha12 = alpha12)
}
