library(mediator)

k1 <- 100 * pi * 0.03^2
k2 <- 100 * pi * 0.03^2
k12maxSq <- min(k1 * k2, (1 - k1) * (1 - k2))
beta_max <- 2 / (0.03^2 + 0.03^2)
(tau2max <- k12maxSq / 100 / 100 / pi^2 / 0.04^4)
sqrt(tau2max)

withr::with_seed(1234, {
  sim3 <- map(1:3, ~ {
    rbidpp(
      n = 100,
      seed = 1234,
      rho1 = 100,
      rho2 = 100,
      tau = 0.5,
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

system.time(
  marginalization1 <- map_dfr(sim3[[1]], estimate, method = "marginalization") %>%
    mutate(method = "marginalization")
)
system.time(
  profiling1 <- map_dfr(sim3[[1]], estimate, method = "profiling") %>%
    mutate(method = "profiling")
)
system.time(
  brute_force1 <- map_dfr(sim3[[1]], estimate, method = "direct") %>%
    mutate(method = "direct")
)
system.time(
  brute_force1alt <- map2_dfr(sim3[[1]], profiling1$fmin, ~ estimate(.x, external_fmin = .y, method = "direct")) %>%
    mutate(method = "direct_fmin")
)
(rmse_tau_mg <- yardstick::rmse_vec(rep(0.5, 100), marginalization1$tau))
(rmse_tau_pr <- yardstick::rmse_vec(rep(0.5, 100), profiling1$tau))
(rmse_tau_bf <- yardstick::rmse_vec(rep(0.5, 100), brute_force1$tau))
(rmse_k12_mg <- yardstick::rmse_vec(rep(0.5 * 100 * pi * 0.035^2, 100), marginalization1$k12))
(rmse_k12_pr <- yardstick::rmse_vec(rep(0.5 * 100 * pi * 0.035^2, 100), profiling1$k12))
(rmse_k12_bf <- yardstick::rmse_vec(rep(0.5 * 100 * pi * 0.035^2, 100), brute_force1$k12))

plot(marginalization1$fmin, profiling1$fmin)
abline(0, 1)
plot(marginalization1$fmin, brute_force1$fmin)
abline(0, 1)
plot(marginalization1$fmin, brute_force1alt$fmin)
abline(0, 1)
plot(profiling1$fmin, brute_force1$fmin)
abline(0, 1)
plot(profiling1$fmin, brute_force1alt$fmin)
abline(0, 1)
plot(brute_force1$fmin, brute_force1alt$fmin)
abline(0, 1)

bind_rows(marginalization1, profiling1, brute_force1, brute_force1alt) %>%
  select(method, fmin) %>%
  group_by(method) %>%
  mutate(id = 1:n()) %>%
  ungroup() %>%
  pivot_wider(names_from = method, values_from = fmin) %>%
  select(-1) %>%
  GGally::ggpairs()

bind_rows(marginalization2, profiling2, brute_force2) %>%
  select(method, fmin) %>%
  group_by(method) %>%
  mutate(id = 1:n()) %>%
  ungroup() %>%
  pivot_wider(names_from = method, values_from = fmin) %>%
  select(-1) %>%
  GGally::ggpairs()
plot(marginalization2$fmin, brute_force2$fmin)
abline(0, 1)
plot(marginalization2$fmin, profiling2$fmin)
abline(0, 1)

bind_rows(marginalization3, profiling3, brute_force3) %>%
  select(method, fmin) %>%
  group_by(method) %>%
  mutate(id = 1:n()) %>%
  ungroup() %>%
  pivot_wider(names_from = method, values_from = fmin) %>%
  select(-1) %>%
  GGally::ggpairs()
plot(marginalization3$fmin, brute_force3$fmin)
abline(0, 1)
plot(marginalization3$fmin, profiling3$fmin)
abline(0, 1)

system.time(
  marginalization2 <- map_dfr(sim3[[2]], estimate, method = "marginalization") %>%
    mutate(method = "marginalization")
)
system.time(
  profiling2 <- map_dfr(sim3[[2]], estimate, method = "profiling") %>%
    mutate(method = "profiling")
)
system.time(
  brute_force2 <- map_dfr(sim3[[2]], estimate, method = "direct") %>%
    mutate(method = "direct")
)
(rmse_tau_mg <- yardstick::rmse_vec(rep(0.5, 100), marginalization2$tau))
(rmse_tau_pr <- yardstick::rmse_vec(rep(0.5, 100), profiling2$tau))
(rmse_tau_bf <- yardstick::rmse_vec(rep(0.5, 100), brute_force2$tau))
(rmse_k12_mg <- yardstick::rmse_vec(rep(0.5 * 100 * pi * 0.035^2, 100), marginalization2$k12))
(rmse_k12_pr <- yardstick::rmse_vec(rep(0.5 * 100 * pi * 0.035^2, 100), profiling2$k12))
(rmse_k12_bf <- yardstick::rmse_vec(rep(0.5 * 100 * pi * 0.035^2, 100), brute_force2$k12))

system.time(
  marginalization3 <- map_dfr(sim3[[3]], estimate, method = "marginalization") %>%
    mutate(method = "marginalization")
)
system.time(
  profiling3 <- map_dfr(sim3[[3]], estimate, method = "profiling") %>%
    mutate(method = "profiling")
)
system.time(
  brute_force3 <- map_dfr(sim3[[3]], estimate, method = "direct") %>%
    mutate(method = "direct")
)
(rmse_tau_mg <- yardstick::rmse_vec(rep(0.5, 100), marginalization3$tau))
(rmse_tau_pr <- yardstick::rmse_vec(rep(0.5, 100), profiling3$tau))
(rmse_tau_bf <- yardstick::rmse_vec(rep(0.5, 100), brute_force3$tau))
(rmse_k12_mg <- yardstick::rmse_vec(rep(0.5 * 100 * pi * 0.035^2, 100), marginalization3$k12))
(rmse_k12_pr <- yardstick::rmse_vec(rep(0.5 * 100 * pi * 0.035^2, 100), profiling3$k12))
(rmse_k12_bf <- yardstick::rmse_vec(rep(0.5 * 100 * pi * 0.035^2, 100), brute_force3$k12))

save.image("data-raw/sim3.RData")

bind_rows(
  marginalization1 %>% mutate(L = 1),
  profiling1 %>% mutate(L = 1),
  brute_force1 %>% mutate(L = 1),
  marginalization2 %>% mutate(L = 2),
  profiling2 %>% mutate(L = 2),
  brute_force2 %>% mutate(L = 2),
  marginalization3 %>% mutate(L = 3),
  profiling3 %>% mutate(L = 3),
  brute_force3 %>% mutate(L = 3)
) %>%
  select(k12, tau, method, L) %>%
  group_by(method, L) %>%
  summarise(
    tau = yardstick::rmse_vec(rep(0.5, 100), tau),
    k12 = yardstick::rmse_vec(rep(0.5 * 100 * pi * 0.035^2, 100), k12)
  ) %>%
  ungroup() %>%
  mutate(time = c(
    104, 134, 285,
    100, 138, 290,
    96, 134, 289
  )) %>%
  pivot_longer(c(tau, k12, time)) %>%
  ggplot(aes(L, value, color = method)) +
  geom_line() +
  geom_point() +
  facet_wrap(vars(name), scales = "free")

# with high tau -----------------------------------------------------------

withr::with_seed(1234, {
  sim3bis <- map(1, ~ {
    rbidpp(
      n = 100,
      seed = 1234,
      rho1 = 100,
      rho2 = 100,
      tau = 0.7,
      alpha1 = 0.03,
      alpha2 = 0.03,
      alpha12 = 0.032,
      L = .x,
      Kspec = "Kspecgauss",
      testtau = "testtaugauss",
      progress = FALSE,
      ncores = 1
    )
  })
})

system.time(
  marginalization1bis <- map_dfr(sim3bis[[1]], estimate, method = "marginalization") %>%
    mutate(method = "marginalization")
)
system.time(
  profiling1bis <- map_dfr(sim3bis[[1]], estimate, method = "profiling") %>%
    mutate(method = "profiling")
)
system.time(
  brute_force1bis <- map_dfr(sim3bis[[1]], estimate, method = "direct") %>%
    mutate(method = "direct")
)

(rmse_tau_mg <- yardstick::rmse_vec(rep(0.7, 100), marginalization1bis$tau))
(rmse_tau_pr <- yardstick::rmse_vec(rep(0.7, 100), profiling1bis$tau))
(rmse_tau_bf <- yardstick::rmse_vec(rep(0.7, 100), brute_force1bis$tau))
(rmse_k12_mg <- yardstick::rmse_vec(rep(0.7 * 100 * pi * 0.032^2, 100), marginalization1bis$k12))
(rmse_k12_pr <- yardstick::rmse_vec(rep(0.7 * 100 * pi * 0.032^2, 100), profiling1bis$k12))
(rmse_k12_bf <- yardstick::rmse_vec(rep(0.7 * 100 * pi * 0.032^2, 100), brute_force1bis$k12))

bind_rows(marginalization1bis, profiling1bis, brute_force1bis) %>%
  select(method, fmin) %>%
  group_by(method) %>%
  mutate(id = 1:n()) %>%
  ungroup() %>%
  pivot_wider(names_from = method, values_from = fmin) %>%
  select(-1) %>%
  GGally::ggpairs()

plot(marginalization1bis$fmin, profiling1bis$fmin)
abline(0, 1)
plot(marginalization1bis$fmin, brute_force1bis$fmin)
abline(0, 1)
plot(profiling1bis$fmin, brute_force1bis$fmin)
abline(0, 1)

# with low tau -----------------------------------------------------------

withr::with_seed(1234, {
  sim3ter <- map(1, ~ {
    rbidpp(
      n = 100,
      seed = 1234,
      rho1 = 100,
      rho2 = 100,
      tau = 0.1,
      alpha1 = 0.03,
      alpha2 = 0.03,
      alpha12 = 0.04,
      L = .x,
      Kspec = "Kspecgauss",
      testtau = "testtaugauss",
      progress = FALSE,
      ncores = 1
    )
  })
})

system.time(
  marginalization1ter <- map_dfr(sim3ter[[1]], estimate, method = "marginalization") %>%
    mutate(method = "marginalization")
)
system.time(
  profiling1ter <- map_dfr(sim3ter[[1]], estimate, method = "profiling") %>%
    mutate(method = "profiling")
)
system.time(
  brute_force1ter <- map_dfr(sim3ter[[1]], estimate, method = "direct") %>%
    mutate(method = "direct")
)

(rmse_tau_mg <- yardstick::rmse_vec(rep(0.1, 100), marginalization1ter$tau))
(rmse_tau_pr <- yardstick::rmse_vec(rep(0.1, 100), profiling1ter$tau))
(rmse_tau_bf <- yardstick::rmse_vec(rep(0.1, 100), brute_force1ter$tau))
(rmse_k12_mg <- yardstick::rmse_vec(rep(0.1 * 100 * pi * 0.04^2, 100), marginalization1ter$k12))
(rmse_k12_pr <- yardstick::rmse_vec(rep(0.1 * 100 * pi * 0.04^2, 100), profiling1ter$k12))
(rmse_k12_bf <- yardstick::rmse_vec(rep(0.1 * 100 * pi * 0.04^2, 100), brute_force1ter$k12))

bind_rows(marginalization1ter, profiling1ter, brute_force1ter) %>%
  select(method, k12, tau) %>%
  pivot_longer(cols = c(k12, tau)) %>%
  ggplot(aes(name, value, fill = method)) +
  geom_boxplot()

plot(log(marginalization1ter$fmin), log(profiling1ter$fmin))
abline(0, 1)
plot(log(marginalization1ter$fmin), log(brute_force1ter$fmin))
abline(0, 1)
plot(profiling1ter$fmin, brute_force1ter$fmin)
abline(0, 1)

system.time(
  marginalization3ter <- map_dfr(sim3ter[[3]], estimate, method = "marginalization") %>%
    mutate(method = "marginalization")
)
system.time(
  profiling1ter <- map_dfr(sim3ter[[1]], estimate, method = "profiling") %>%
    mutate(method = "profiling")
)
system.time(
  brute_force1ter <- map_dfr(sim3ter[[1]], estimate, method = "direct") %>%
    mutate(method = "direct")
)

(rmse_tau_mg <- yardstick::rmse_vec(rep(0.1, 100), marginalization1ter$tau))
(rmse_tau_pr <- yardstick::rmse_vec(rep(0.1, 100), profiling1ter$tau))
(rmse_tau_bf <- yardstick::rmse_vec(rep(0.1, 100), brute_force1ter$tau))
(rmse_k12_mg <- yardstick::rmse_vec(rep(0.1 * 100 * pi * 0.04^2, 100), marginalization1ter$k12))
(rmse_k12_pr <- yardstick::rmse_vec(rep(0.1 * 100 * pi * 0.04^2, 100), profiling1ter$k12))
(rmse_k12_bf <- yardstick::rmse_vec(rep(0.1 * 100 * pi * 0.04^2, 100), brute_force1ter$k12))

bind_rows(marginalization1ter, profiling1ter, brute_force1ter) %>%
  select(method, k12, tau) %>%
  pivot_longer(cols = c(k12, tau)) %>%
  ggplot(aes(name, value, fill = method)) +
  geom_boxplot()

plot(log(marginalization1ter$fmin), log(profiling1ter$fmin))
abline(0, 1)
plot(log(marginalization1ter$fmin), log(brute_force1ter$fmin))
abline(0, 1)
plot(profiling1ter$fmin, brute_force1ter$fmin)
abline(0, 1)
