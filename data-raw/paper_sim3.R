library(muscadet)
library(future)

rho1 <- rho2 <- 100
alpha1 <- alpha2 <- 0.03
alpha12 <- 0.035
tau <- 0.5
k1 <- rho1 * pi * alpha1^2
k2 <- rho2 * pi * alpha2^2
k12 <- tau * sqrt(rho1 * rho2) * pi * alpha12^2
beta <- 1 / alpha12^2

sim3_full <- map(c(0, 0.7, by = 0.1), function(tau) {
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

load("data-raw/sim3.RData")

plan(multisession, workers = 4)
df2 <- bind_rows(
  furrr::future_map_dfr(sim3[[1]], ~ {
    estimate(
      X = .x,
      model = "Gauss",
      rmin_alpha = 2,
      rmin_alpha12 = 2,
      rmin_tau = 2,
      q = 1,
      p = 2,
      divisor_marginal = "d",
      divisor_cross = "d",
      method = "profiling")
  }) %>% mutate(L = 1),
  furrr::future_map_dfr(sim3[[2]], ~ {
    estimate(
      X = .x,
      model = "Gauss",
      rmin_alpha = 2,
      rmin_alpha12 = 2,
      rmin_tau = 2,
      q = 1,
      p = 2,
      divisor_marginal = "d",
      divisor_cross = "d",
      method = "profiling")
  }) %>% mutate(L = 2),
  furrr::future_map_dfr(sim3[[3]], ~ {
    estimate(
      X = .x,
      model = "Gauss",
      rmin_alpha = 2,
      rmin_alpha12 = 2,
      rmin_tau = 2,
      q = 1,
      p = 2,
      divisor_marginal = "d",
      divisor_cross = "d",
      method = "profiling")
  }) %>% mutate(L = 3)
)
plan(sequential)

df2 %>%
  select(alpha1, alpha2, k12, tau, L) %>%
  group_by(L) %>%
  summarise(
    alpha1 = yardstick::rmse_vec(rep(0.03, 100), alpha1),
    alpha2 = yardstick::rmse_vec(rep(0.03, 100), alpha2),
    tau = yardstick::rmse_vec(rep(0.5, 100), tau),
    k12 = yardstick::rmse_vec(rep(0.5 * 100 * pi * 0.035^2, 100), k12)
  ) %>%
  pivot_longer(-L) %>%
  ggplot(aes(L, value, color = name)) +
  geom_line() +
  geom_point() +
  facet_wrap(vars(name), nrow = 1, scales = "free")

df2 %>%
  select(alpha1, alpha2, k12, tau, L) %>%
  pivot_longer(-L) %>%
  ggplot(aes(as.factor(L), value, fill = as.factor(L))) +
  geom_boxplot() +
  geom_hline(aes(yintercept = value), tibble(
    name = c("alpha1", "alpha2", "k12", "tau"),
    value = c(0.03, 0.03, 0.5 * 100 * pi * 0.035^2, 0.5)
  )) +
  facet_wrap(vars(name), nrow = 1, scales = "free")
