AVAILABLE_MODELS <<- c("Gauss", "Bessel")

jinc <- function(x, alpha) {
  # J_alpha(x) / (x/2)^alpha * gamma(alpha + 1)
  if (x < sqrt(.Machine$double.eps))
    return(1)
  besselJ(x, alpha) / (x / 2)^alpha * gamma(alpha + 1)
}
