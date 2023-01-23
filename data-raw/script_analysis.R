library(muscadet)
library(future)
progressr::handlers("rstudio")

.fit_data <- function(sims, model = "Gauss", full_bootstrap = TRUE) {
  sims <- sims[1:1008]
  pb <- progressr::progressor(along = sims)
  furrr::future_map(sims, ~ {
    pb()
    fit_via_pcf(.x, model = model, B = 19, full_bootstrap = full_bootstrap)
  }, .options = furrr::furrr_options(seed = TRUE))
}

# Gauss1 ------------------------------------------------------------------

sims_gauss1 <- read_rds("~/Documents/Professionnel/Data/DPP/sims_gauss1.rds")

plan(multisession, workers = parallelly::availableCores(omit = 1))
withr::with_seed(1234, {
  sims_gauss1 <- sims_gauss1 |>
    mutate(
      full_bootstrap = map(data, ~ {
        progressr::with_progress(
          res <- .fit_data(.x, model = "Gauss", full_bootstrap = TRUE)
        )
        res
      }),
      partial_bootstrap = map(data, ~ {
        progressr::with_progress(
          res <- .fit_data(.x, model = "Gauss", full_bootstrap = FALSE)
        )
        res
      })
    )
})
plan(sequential)

write_rds(sims_gauss1, file = "~/Documents/Professionnel/Data/DPP/sims_gauss1_R1008.rds")

# Bessel1 -----------------------------------------------------------------

sims_bessel1 <- read_rds("~/Documents/Professionnel/Data/DPP/sims_bessel1.rds")

plan(multisession, workers = parallelly::availableCores(omit = 1))
withr::with_seed(1234, {
  sims_bessel1 <- sims_bessel1 |>
    mutate(
      full_bootstrap = map(data, ~ {
        progressr::with_progress(
          res <- .fit_data(.x, model = "Bessel", full_bootstrap = TRUE)
        )
        res
      }),
      partial_bootstrap = map(data, ~ {
        progressr::with_progress(
          res <- .fit_data(.x, model = "Bessel", full_bootstrap = FALSE)
        )
        res
      })
    )
})
plan(sequential)

write_rds(sims_bessel1, file = "~/Documents/Professionnel/Data/DPP/sims_bessel1_R1008.rds")

# Gauss2 -----------------------------------------------------------------

sims_gauss2 <- read_rds("~/Documents/Professionnel/Data/DPP/sims_gauss2.rds")

plan(multisession, workers = parallelly::availableCores(omit = 1))
withr::with_seed(1234, {
  sims_gauss2 <- sims_gauss2 |>
    mutate(
      full_bootstrap = map(data, ~ {
        progressr::with_progress(
          res <- .fit_data(.x, model = "Gauss", full_bootstrap = TRUE)
        )
        res
      }),
      partial_bootstrap = map(data, ~ {
        progressr::with_progress(
          res <- .fit_data(.x, model = "Gauss", full_bootstrap = FALSE)
        )
        res
      })
    )
})
plan(sequential)

write_rds(sims_gauss2, file = "~/Documents/Professionnel/Data/DPP/sims_gauss2_R1008.rds")

# 1600 pour 90 ppp a 800 pts sur 9 cores
# 162 sec pour 1 ppp a 800 pts
# 195 pour 9 ppp 800 pts en parallelisant sur 9 cores
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




alpha <- 0.05

sim_gauss1 <- read_rds("~/Documents/Professionnel/Data/DPP/sims_gauss1_R1008.rds")
sim_gauss1 |>
  select(-data) |>
  mutate(
    np_full_bootstrap = purrr::map_dbl(full_bootstrap, ~ mean(purrr::map_dbl(.x, function(boot_data) {
      (1 + sum(boot_data$stat_np_boots >= boot_data$stat_np_obs)) / (1 + length(boot_data$stat_np_boots))
    }) <= alpha)),
    np_partial_bootstrap = purrr::map_dbl(partial_bootstrap, ~ mean(purrr::map_dbl(.x, function(boot_data) {
      (1 + sum(boot_data$stat_np_boots >= boot_data$stat_np_obs)) / (1 + length(boot_data$stat_np_boots))
    }) <= alpha)),
    p_full_bootstrap = purrr::map_dbl(full_bootstrap, ~ mean(purrr::map_dbl(.x, function(boot_data) {
      (1 + sum(boot_data$stat_p_boots >= boot_data$stat_p_obs)) / (1 + length(boot_data$stat_p_boots))
    }) <= alpha)),
    p_partial_bootstrap = purrr::map_dbl(partial_bootstrap, ~ mean(purrr::map_dbl(.x, function(boot_data) {
      (1 + sum(boot_data$stat_p_boots >= boot_data$stat_p_obs)) / (1 + length(boot_data$stat_p_boots))
    }) <= alpha))
  ) |>
  select(-full_bootstrap, -partial_bootstrap) |>
  pivot_longer(ends_with("bootstrap")) |>
  ggplot(aes(tau1, value, color = name)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = alpha, linetype = "dashed") +
  facet_wrap(vars(L), nrow = 1) +
  labs(title = "Gaussian Model", subtitle = "First set of parameters (alpha=0.03, alpha12=0.035)")

sim_bessel1 <- read_rds("~/Documents/Professionnel/Data/DPP/sims_bessel1_R1008.rds")
sim_bessel1 |>
  select(-data) |>
  mutate(
    np_full_bootstrap = purrr::map_dbl(full_bootstrap, ~ mean(purrr::map_dbl(.x, function(boot_data) {
      (1 + sum(boot_data$stat_np_boots >= boot_data$stat_np_obs)) / (1 + length(boot_data$stat_np_boots))
    }) <= alpha)),
    np_partial_bootstrap = purrr::map_dbl(partial_bootstrap, ~ mean(purrr::map_dbl(.x, function(boot_data) {
      (1 + sum(boot_data$stat_np_boots >= boot_data$stat_np_obs)) / (1 + length(boot_data$stat_np_boots))
    }) <= alpha)),
    p_full_bootstrap = purrr::map_dbl(full_bootstrap, ~ mean(purrr::map_dbl(.x, function(boot_data) {
      (1 + sum(boot_data$stat_p_boots >= boot_data$stat_p_obs)) / (1 + length(boot_data$stat_p_boots))
    }) <= alpha)),
    p_partial_bootstrap = purrr::map_dbl(partial_bootstrap, ~ mean(purrr::map_dbl(.x, function(boot_data) {
      (1 + sum(boot_data$stat_p_boots >= boot_data$stat_p_obs)) / (1 + length(boot_data$stat_p_boots))
    }) <= alpha))
  ) |>
  select(-full_bootstrap, -partial_bootstrap) |>
  pivot_longer(ends_with("bootstrap")) |>
  ggplot(aes(tau1, value, color = name)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = alpha, linetype = "dashed") +
  facet_wrap(vars(L), nrow = 1) +
  labs(title = "Bessel Model", subtitle = "First set of parameters (alpha=0.03, alpha12=0.035)")

sim_gauss2 <- read_rds("~/Documents/Professionnel/Data/DPP/sims_gauss2_R1008.rds")
sim_gauss2 |>
  select(-data) |>
  mutate(
    np_full_bootstrap = purrr::map_dbl(full_bootstrap, ~ mean(purrr::map_dbl(.x, function(boot_data) {
      (1 + sum(boot_data$stat_np_boots >= boot_data$stat_np_obs)) / (1 + length(boot_data$stat_np_boots))
    }) <= alpha)),
    np_partial_bootstrap = purrr::map_dbl(partial_bootstrap, ~ mean(purrr::map_dbl(.x, function(boot_data) {
      (1 + sum(boot_data$stat_np_boots >= boot_data$stat_np_obs)) / (1 + length(boot_data$stat_np_boots))
    }) <= alpha)),
    p_full_bootstrap = purrr::map_dbl(full_bootstrap, ~ mean(purrr::map_dbl(.x, function(boot_data) {
      (1 + sum(boot_data$stat_p_boots >= boot_data$stat_p_obs)) / (1 + length(boot_data$stat_p_boots))
    }) <= alpha)),
    p_partial_bootstrap = purrr::map_dbl(partial_bootstrap, ~ mean(purrr::map_dbl(.x, function(boot_data) {
      (1 + sum(boot_data$stat_p_boots >= boot_data$stat_p_obs)) / (1 + length(boot_data$stat_p_boots))
    }) <= alpha))
  ) |>
  select(-full_bootstrap, -partial_bootstrap) |>
  pivot_longer(ends_with("bootstrap")) |>
  ggplot(aes(tau2, value, color = name)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = alpha, linetype = "dashed") +
  facet_wrap(vars(L), nrow = 1) +
  labs(title = "Gaussian Model", subtitle = "Second set of parameters (alpha=0.04, alpha12=0.041)")

sim_bessel2 <- read_rds("~/Documents/Professionnel/Data/DPP/sims_bessel2_R1008.rds")
sim_bessel2 |>
  select(-data) |>
  mutate(
    np_full_bootstrap = purrr::map_dbl(full_bootstrap, ~ mean(purrr::map_dbl(.x, function(boot_data) {
      (1 + sum(boot_data$stat_np_boots > boot_data$stat_np_obs)) / (1 + length(boot_data$stat_np_boots))
    }) <= alpha)),
    np_partial_bootstrap = purrr::map_dbl(partial_bootstrap, ~ mean(purrr::map_dbl(.x, function(boot_data) {
      (1 + sum(boot_data$stat_np_boots > boot_data$stat_np_obs)) / (1 + length(boot_data$stat_np_boots))
    }) <= alpha)),
    p_full_bootstrap = purrr::map_dbl(full_bootstrap, ~ mean(purrr::map_dbl(.x, function(boot_data) {
      (1 + sum(boot_data$stat_p_boots > boot_data$stat_p_obs)) / (1 + length(boot_data$stat_p_boots))
    }) <= alpha)),
    p_partial_bootstrap = purrr::map_dbl(partial_bootstrap, ~ mean(purrr::map_dbl(.x, function(boot_data) {
      (1 + sum(boot_data$stat_p_boots > boot_data$stat_p_obs)) / (1 + length(boot_data$stat_p_boots))
    }) <= alpha))
  ) |>
  select(-full_bootstrap, -partial_bootstrap) |>
  pivot_longer(ends_with("bootstrap")) |>
  ggplot(aes(tau2, value, color = name)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = alpha, linetype = "dashed") +
  facet_wrap(vars(L), nrow = 1) +
  labs(title = "Bessel Model")
