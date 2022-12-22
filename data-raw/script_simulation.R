library(muscadet)
library(future)

mod_gauss1 <- TwoMarkDPP$new(model = "Gauss")
mod_gauss1$first_intensity <- 100
mod_gauss1$second_intensity <- 100
mod_gauss1$first_repulsion_rate <- 0.03
mod_gauss1$second_repulsion_rate <- 0.03
mod_gauss1$cross_repulsion_rate <- 0.035
mod_gauss1

mod_gauss2 <- TwoMarkDPP$new(model = "Gauss")
mod_gauss2$first_intensity <- 100
mod_gauss2$second_intensity <- 100
mod_gauss2$first_repulsion_rate <- 0.04
mod_gauss2$second_repulsion_rate <- 0.04
mod_gauss2$cross_repulsion_rate <- 0.041
mod_gauss2

mod_bessel1 <- TwoMarkDPP$new(model = "Bessel")
mod_bessel1$first_intensity <- 100
mod_bessel1$second_intensity <- 100
mod_bessel1$first_repulsion_rate <- 0.03
mod_bessel1$second_repulsion_rate <- 0.03
mod_bessel1$cross_repulsion_rate <- 0.035
mod_bessel1

mod_bessel2 <- TwoMarkDPP$new(model = "Bessel")
mod_bessel2$first_intensity <- 100
mod_bessel2$second_intensity <- 100
mod_bessel2$first_repulsion_rate <- 0.04
mod_bessel2$second_repulsion_rate <- 0.04
mod_bessel2$cross_repulsion_rate <- 0.041
mod_bessel2

R <- 10008
L <- 1:2
tau1 <- seq(0, 0.7, by = 0.1)
tau2 <- seq(0, 0.9, by = 0.1)

# GAUSS1 ------------------------------------------------------------------

run_gauss1 <- function(n, ws, tau) {
  cli::cli_alert_info("Generating {n} data sets with correlation of {tau} between marks on a window size L = {ws}.")
  mod_gauss1$window_size <- ws
  mod_gauss1$between_mark_correlation <- tau
  withr::with_seed(1234, {
    progressr::with_progress(
      mod_gauss1$random(n = n)
    )
  })
}

progressr::handlers("rstudio")
plan(multisession, workers = parallelly::availableCores(omit = 1))
sims_gauss1 <- crossing(L, tau1) |>
  mutate(data = map2(L, tau1, run_gauss1, n = R))
plan(sequential)
write_rds(sims_gauss1, file = "~/Documents/Professionnel/Data/DPP/sims_gauss1.rds", compress = "xz", version = 3)

# BESSEL1 -----------------------------------------------------------------

run_bessel1 <- function(n, ws, tau) {
  cli::cli_alert_info("Generating {n} data sets with correlation of {tau} between marks on a window size L = {ws}.")
  mod_bessel1$window_size <- ws
  mod_bessel1$between_mark_correlation <- tau
  withr::with_seed(1234, {
    progressr::with_progress(
      mod_bessel1$random(n = n)
    )
  })
}

progressr::handlers("rstudio")
plan(multisession, workers = parallelly::availableCores(omit = 1))
sims_bessel1 <- crossing(L, tau1) |>
  mutate(data = map2(L, tau1, run_bessel1, n = R))
plan(sequential)
write_rds(sims_bessel1, file = "~/Documents/Professionnel/Data/DPP/sims_bessel1.rds", compress = "xz", version = 3)

# GAUSS2 ------------------------------------------------------------------

run_gauss2 <- function(n, ws, tau) {
  cli::cli_alert_info("Generating {n} data sets with correlation of {tau} between marks on a window size L = {ws}.")
  mod_gauss2$window_size <- ws
  mod_gauss2$between_mark_correlation <- tau
  withr::with_seed(1234, {
    progressr::with_progress(
      mod_gauss2$random(n = n)
    )
  })
}

progressr::handlers("rstudio")
plan(multisession, workers = parallelly::availableCores(omit = 1))
sims_gauss2 <- crossing(L, tau2) |>
  mutate(data = map2(L, tau2, run_gauss2, n = R))
plan(sequential)
write_rds(sims_gauss2, file = "~/Documents/Professionnel/Data/DPP/sims_gauss2.rds", compress = "xz", version = 3)

# BESSEL2 -----------------------------------------------------------------

run_bessel2 <- function(n, ws, tau) {
  cli::cli_alert_info("Generating {n} data sets with correlation of {tau} between marks on a window size L = {ws}.")
  mod_bessel2$window_size <- ws
  mod_bessel2$between_mark_correlation <- tau
  withr::with_seed(1234, {
    progressr::with_progress(
      mod_bessel2$random(n = n)
    )
  })
}

progressr::handlers("rstudio")
plan(multisession, workers = parallelly::availableCores(omit = 1))
sims_bessel2 <- crossing(L, tau2) |>
  mutate(data = map2(L, tau2, run_bessel2, n = R))
plan(sequential)
write_rds(sims_bessel2, file = "~/Documents/Professionnel/Data/DPP/sims_bessel2.rds", compress = "xz", version = 3)
