sims_gauss1 <- read_rds("~/Documents/Professionnel/Data/DPP/sims_gauss1.rds")
# sims_bessel1 <- read_rds("~/Documents/Professionnel/Data/DPP/sims_bessel1.rds")
library(muscadet)
library(future)

plan(multisession, workers = parallelly::availableCores(omit = 1))
system.time(
  withr::with_seed(1234, {
    test <- sims_gauss1$data[[9]][1:9] |>
      furrr::future_map(
        .f = fit_via_pcf,
        model = "Gauss",
        B = 19,
        full_bootstrap = FALSE,
        .options = furrr::furrr_options(seed = TRUE)
      )
  })
)
plan(sequential)

# 162 sec pour 1 ppp a 800 pts
# 237 pour 9 ppp 800 pts en parallelisant sur 9 cores
# 256s pour le 1 (400 pts)
# 265s pour le 9 (800 pts)

alpha <- 0.05
np_pvals <- map_dbl(test, "pvalue_np")
mean(np_pvals <= alpha)
p_pvals <- map_dbl(test, "pvalue_p")
mean(p_pvals <= alpha)


mod <- TwoMarkDPP$new()
mod$first_intensity <- 100
mod$second_intensity <- 100
mod$first_repulsion_rate <- 0.03
mod$second_repulsion_rate <- 0.03
mod$cross_repulsion_rate <- 0.035
mod$between_mark_correlation <- 0
mod$window_size <- 2
mod$random(1)
