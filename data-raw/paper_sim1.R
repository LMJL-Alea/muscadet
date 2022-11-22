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

load("data-raw/sim3.RData")

plan(multisession, workers = 4)
df1 <- crossing(
  rmin_alpha = c(2, 32),
  rmin_alpha12 = c(2, 32),
  rmin_tau = c(2, 32),
  q = c(0.5, 1),
  divisor_marginal = c("r", "d"),
  divisor_cross = c("r", "d"),
  method = c("integration", "ip", "profiling", "direct")
) |>
  filter(!(method != "direct" & q == 0.5)) |>
  mutate(
    mod = furrr::future_pmap(list(rmin_alpha, rmin_alpha12, rmin_tau, q, divisor_marginal, divisor_cross, method), ~ {
      map_dfr(sim3[[1]], function(data) {
        estimate(
          X = data,
          model = "Gauss",
          rmin_alpha = ..1,
          rmin_alpha12 = ..2,
          rmin_tau = ..3,
          q = ..4,
          p = 2,
          divisor_marginal = ..5,
          divisor_cross = ..6,
          method = ..7)
      })
    })
  )
plan(sequential)

df1 |>
  unnest(mod) |>
  select(-c(alpha12, fmin)) |>
  group_by(rmin_alpha, rmin_alpha12, rmin_tau, q, divisor_marginal, divisor_cross, method) |>
  summarise(
    rho1 = yardstick::rmse_vec(rep(100, 100), rho1),
    rho2 = yardstick::rmse_vec(rep(100, 100), rho2),
    alpha1 = yardstick::rmse_vec(rep(0.03, 100), alpha1),
    alpha2 = yardstick::rmse_vec(rep(0.03, 100), alpha2),
    tau = yardstick::rmse_vec(rep(0.5, 100), tau),
    k12 = yardstick::rmse_vec(rep(0.5 * 100 * pi * 0.035^2, 100), k12)
  ) |>
  select(-rho1, -rho2) |>
  pivot_longer(alpha1:k12) |>
  filter(method == "direct") |>
  ggplot(aes(interaction(rmin_alpha, rmin_alpha12, rmin_tau), value, fill = as.factor(q))) +
  geom_col(position = position_dodge(), color = "black") +
  facet_grid(rows = vars(name), cols = vars(divisor_marginal, divisor_cross), scales = "free")

# useless to keep q=0.5 so filter out

df1 |>
  unnest(mod) |>
  select(-c(alpha12, fmin)) |>
  group_by(rmin_alpha, rmin_alpha12, rmin_tau, q, divisor_marginal, divisor_cross, method) |>
  summarise(
    rho1 = yardstick::rmse_vec(rep(100, 100), rho1),
    rho2 = yardstick::rmse_vec(rep(100, 100), rho2),
    alpha1 = yardstick::rmse_vec(rep(0.03, 100), alpha1),
    alpha2 = yardstick::rmse_vec(rep(0.03, 100), alpha2),
    tau = yardstick::rmse_vec(rep(0.5, 100), tau),
    k12 = yardstick::rmse_vec(rep(0.5 * 100 * pi * 0.035^2, 100), k12)
  ) |>
  select(-rho1, -rho2) |>
  pivot_longer(alpha1:k12) |>
  filter(q == 1, method %in% c("profiling", "direct")) |>
  ggplot(aes(interaction(rmin_alpha, rmin_alpha12, rmin_tau), value, fill = method)) +
  geom_col(color = "black", position = position_dodge()) +
  facet_grid(cols = vars(divisor_marginal, divisor_cross), rows = vars(name), scales = "free")

# --> profiling, divisor_marginal = d, rmin_tau = 2

df1 |>
  unnest(mod) |>
  select(-c(alpha12, fmin)) |>
  group_by(rmin_alpha, rmin_alpha12, rmin_tau, q, divisor_marginal, divisor_cross, method) |>
  summarise(
    rho1 = yardstick::rmse_vec(rep(100, 100), rho1),
    rho2 = yardstick::rmse_vec(rep(100, 100), rho2),
    alpha1 = yardstick::rmse_vec(rep(0.03, 100), alpha1),
    alpha2 = yardstick::rmse_vec(rep(0.03, 100), alpha2),
    tau = yardstick::rmse_vec(rep(0.5, 100), tau),
    k12 = yardstick::rmse_vec(rep(0.5 * 100 * pi * 0.035^2, 100), k12)
  ) |>
  select(-rho1, -rho2, -alpha1, -alpha2) |>
  pivot_longer(tau:k12) |>
  filter(q == 1, method == "profiling", divisor_marginal == "d", rmin_tau == 2) |>
  ggplot(aes(interaction(rmin_alpha, rmin_alpha12), value)) +
  geom_col(color = "black", position = position_dodge()) +
  facet_grid(cols = vars(divisor_cross), rows = vars(name), scales = "free")

# --> divisor_cross = d, rmin_alpha = 2, rmin_alpha12 = 2

df1 |>
  unnest(mod) |>
  select(-c(alpha12, fmin)) |>
  group_by(rmin_alpha, rmin_alpha12, rmin_tau, q, divisor_marginal, divisor_cross, method) |>
  summarise(
    rho1 = yardstick::rmse_vec(rep(100, 100), rho1),
    rho2 = yardstick::rmse_vec(rep(100, 100), rho2),
    alpha1 = yardstick::rmse_vec(rep(0.03, 100), alpha1),
    alpha2 = yardstick::rmse_vec(rep(0.03, 100), alpha2),
    tau = yardstick::rmse_vec(rep(0.5, 100), tau),
    k12 = yardstick::rmse_vec(rep(0.5 * 100 * pi * 0.035^2, 100), k12)
  ) |>
  select(-rho1, -rho2) |>
  pivot_longer(alpha1:k12) |>
  filter(q == 1, method == "profiling", divisor_marginal == "d") |>
  ggplot(aes(interaction(rmin_alpha, rmin_alpha12, rmin_tau), value)) +
  geom_col(color = "black", position = position_dodge()) +
  facet_grid(cols = vars(divisor_marginal, divisor_cross), rows = vars(name), scales = "free")
