# Define fixed parameters
rho1 <- 100
rho2 <- 100
alpha1 <- 0.03
alpha2 <- 0.03
alpha12 <- 0.035
R <- 100

# Find correlation upper bound
k1 <- rho1 * pi * alpha1^2
k2 <- rho2 * pi * alpha2^2
max_k12 <- sqrt(max(min(k1 * k2, (1 - k1) * (1 - k2) - sqrt(.Machine$double.eps)), 0))
max_tau <- max_k12 / sqrt(rho1 * rho2) / (pi * alpha12^2)
tau <- seq(0, floor(10 * max_tau) / 10, by = 0.1)

# Run simulations
sims <- purrr::imap(tau, ~ {
  writeLines(paste("- Simulating bivariate Gaussian DPP with correlation", .x))
  sim <- mediator::simulate(
    n = R,
    seed = 1234,
    rho1 = rho1,
    rho2 = rho2,
    tau = .x,
    alpha1 = alpha1,
    alpha2 = alpha2,
    alpha12 = alpha12,
    progress = TRUE,
    Kspec = "Kspecgauss",
    testtau = "testtaugauss"
  )
  dataset <- paste0("sim_gauss", .y - 1)
  assign(dataset, sim)
  save(
    list = dataset,
    file = paste0("data/", dataset, ".rda"),
    compress = "xz"
  )
})

library(tidyverse)
library(mediator)
new_simulations <- tibble(
  rho1 = 100,
  rho2 = 100,
  alpha1 = c(0.05, 0.01),
  alpha2 = c(0.05, 0.01),
  alpha12 = c(0.05, 0.012),
  tau = c(0.2, 0.5)
) %>%
  mutate(
    sim = list(rho1, rho2, alpha1, alpha2, alpha12, tau) %>%
      pmap(~ {
        rbidpp(
          n = 100,
          seed = 1234,
          rho1 = ..1,
          rho2 = ..2,
          tau = ..6,
          alpha1 = ..3,
          alpha2 = ..4,
          alpha12 = ..5,
          progress = FALSE,
          Kspec = "Kspecgauss",
          testtau = "testtaugauss"
        )
      })
  )
