library(muscadet)
library(future)
library(progressr)
progressr::handlers("rstudio")

R <- 1000

data_gen <- function(tau_values) {
  p <- progressor(along = tau_values)
  furrr::future_map(tau_values, function(tau) {
    p()
    rbidpp(
      n = R,
      seed = 1234,
      rho1 = 100,
      rho2 = 100,
      tau = tau,
      alpha1 = 0.03,
      alpha2 = 0.03,
      alpha12 = 0.035,
      L = 1,
      Kspec = "Kspecgauss",
      testtau = "testtaugauss",
      progress = FALSE,
      ncores = 1
    )
  })
}

tau_values <- seq(0, 0.7, by = 0.1)
plan(multisession, workers = 4)
with_progress({
  sim_L1 <- data_gen(tau_values)
})
plan(sequential)
save(sim_L1, file = "data-raw/sim_L1.RData")
