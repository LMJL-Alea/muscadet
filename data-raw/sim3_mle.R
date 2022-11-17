load("data-raw/sim3.RData")
sol1 <- map(sim3[[1]], estimate) # PCF
sol2 <- map(sim3[[1]], mle, num_threads = 4) # MLE-4-RandomInit
sol3 <- map2(sim3[[1]], sol1, mle, num_threads = 4) # MLE-4-PCFInit
sol4 <- map(sim3[[1]], mle, num_threads = 4, optimizer = "neldermead") # MLE-4-RandomInit
sol5 <- map2(sim3[[1]], sol1, mle, num_threads = 4, optimizer = "neldermead") # MLE-4-PCFInit

sol4 <- map(sim3[[1]], mle, num_threads = 4, fixed_marginal_parameters = TRUE) # MLE-2-RandomInit
sol5 <- map2(sim3[[1]], sol1, mle, num_threads = 4, fixed_marginal_parameters = TRUE) # MLE-2-PCFInit

# tau ---------------------------------------------------------------------

tibble(
  pcf = map_dbl(sol1, "tau"),
  mle4_irand = map(sol2, "par") %>% map_dbl(6),
  mle4_iPCF = map(sol3, "par") %>% map_dbl(6),
  mle4_irand_f = map(sol4, "par") %>% map_dbl(6),
  # mle4_iPCF_f = map(sol5, "par") %>% map_dbl(6)
) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(name, value)) +
  geom_boxplot()

yardstick::rmse_vec(rep(0.5, 100), map_dbl(sol1, "tau"))
yardstick::rmse_vec(rep(0.5, 100), map(sol2, "par") %>% map_dbl(6))
yardstick::rmse_vec(rep(0.5, 100), map(sol3, "par") %>% map_dbl(6))
yardstick::rmse_vec(rep(0.5, 100), map(sol4, "par") %>% map_dbl(6))

# alpha1 ------------------------------------------------------------------

tibble(
  pcf = map_dbl(sol1, "alpha1"),
  mle4_irand = map(sol2, "par") %>% map_dbl(3),
  mle4_iPCF = map(sol3, "par") %>% map_dbl(3)
) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(name, value)) +
  geom_boxplot()

yardstick::rmse_vec(rep(0.03, 100), map_dbl(sol1, "alpha1")) / 0.03 * 100
yardstick::rmse_vec(rep(0.03, 100), map(sol2, "par") %>% map_dbl(3)) / 0.03 * 100
yardstick::rmse_vec(rep(0.03, 100), map(sol3, "par") %>% map_dbl(3)) / 0.03 * 100

# alpha2 ------------------------------------------------------------------

tibble(
  pcf = map_dbl(sol1, "alpha2"),
  mle4_irand = map(sol2, "par") %>% map_dbl(4),
  mle4_iPCF = map(sol3, "par") %>% map_dbl(4),
  mle4_irand_f = map(sol4, "par") %>% map_dbl(4),
) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(name, value)) +
  geom_boxplot()

yardstick::rmse_vec(rep(0.03, 100), map_dbl(sol1, "alpha2")) / 0.03 * 100
yardstick::rmse_vec(rep(0.03, 100), map(sol2, "par") %>% map_dbl(4)) / 0.03 * 100
yardstick::rmse_vec(rep(0.03, 100), map(sol3, "par") %>% map_dbl(4)) / 0.03 * 100

