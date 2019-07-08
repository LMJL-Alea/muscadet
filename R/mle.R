#' @export
mle_dpp <- function(X, labels, lb, ub, rho1, rho2, alpha1, alpha2) {
  # optim(par = 0.5, fn = Evaluate, method = "L-BFGS-B", lower = 1e-4, upper = 1-1e-4,
  #       X = X, labels = labels, lb = lb, ub = ub, rho1 = rho1, rho2 = rho2, alpha1 = alpha1, alpha2 = alpha2)
  # optim(par = 0, fn = Evaluate, method = "L-BFGS-B", lower = -0.28, upper = 0.28,
  #       X = X, labels = labels, lb = lb, ub = ub, rho1 = rho1, rho2 = rho2, alpha1 = alpha1, alpha2 = alpha2)
  # optim(par = 0.02, fn = Evaluate, method = "L-BFGS-B", lower = 1e-4, upper = 0.053,
  #       X = X, labels = labels, lb = lb, ub = ub, rho1 = rho1, rho2 = rho2, alpha1 = alpha1, alpha2 = alpha2)
  optim(par = 0.05, fn = Evaluate, method = "L-BFGS-B", lower = 0.03,
        X = X, labels = labels, lb = lb, ub = ub, rho1 = rho1, rho2 = rho2, alpha1 = alpha1, alpha2 = alpha2)
}
