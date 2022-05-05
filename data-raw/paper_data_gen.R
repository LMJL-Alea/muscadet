library(muscadet)
library(future)
library(progressr)

data_gen <- function(tau_values) {
  p <- progressor(along = tau_values)
  sim3_full <- furrr::future_map(tau_values, function(tau) {
    p()
    purrr::map(1:3, ~ {
      rbidpp(
        n = 100,
        seed = 1234,
        rho1 = 100,
        rho2 = 100,
        tau = tau,
        alpha1 = 0.03,
        alpha2 = 0.03,
        alpha12 = 0.035,
        L = .x,
        Kspec = "Kspecgauss",
        testtau = "testtaugauss",
        progress = FALSE,
        ncores = 1
      )
    })
  })

}

tau_values <- c(0, 0.7, by = 0.1)
plan(multisession, workers = 4)
df <- data_gen(tau_values)
plan(sequential)
