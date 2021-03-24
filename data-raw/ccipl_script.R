library(tidyverse)
library(mediator)

load("simulations.rda")

simulations <- mutate(simulations, pcf = map(sim, estimate))

save(simulations, file = "simulations_pcf.rda", version = 3, compress = "xz")

ncores <- 20

mle_grouped <- function(data_list,
                        initial_guess = NULL,
                        optimizer = "bobyqa",
                        fixed_marginal_parameters = FALSE,
                        num_threads = ncores) {
  if (is.null(initial_guess)) {
    fit <- data_list %>%
      imap(~ {
        print(paste0("Estimating parameters for data set #", .y))
        mle(
          data = .x,
          fixed_marginal_parameters = fixed_marginal_parameters,
          num_threads = num_threads, optimizer = optimizer
        )
      })
  } else {
    fit <- data_list %>%
      imap(~ {
        print(paste0("Estimating parameters for data set #", .y))
        mle(
          data = .x,
          fixed_marginal_parameters = fixed_marginal_parameters,
          num_threads = num_threads,
          initial_guess = initial_guess[[.y]], optimizer = optimizer
        )
      })
  }
  fit %>%
    map_df("par") %>%
    mutate(fmin = map_dbl(fit, "value"))
}

simulations_L1 <- simulations[1:3, ] %>%
  mutate(
    mle2_imid = map(sim, mle_grouped,
                    fixed_marginal_parameters = TRUE),
    mle2_ipcf = map(sim, mle_grouped,
                    fixed_marginal_parameters = TRUE,
                    initial_guess = pcf),
    mle4_imid = map(sim, mle_grouped,
                    fixed_marginal_parameters = FALSE),
    mle4_imid = map(sim, mle_grouped,
                    fixed_marginal_parameters = FALSE,
                    initial_guess = pcf)
  )

save(simulations_L1, file = "mle_L1.rda", version = 3, compress = "xz")

simulations_L2 <- simulations[4:6, ] %>%
  mutate(
    mle2_imid = map(sim, mle_grouped,
                    fixed_marginal_parameters = TRUE),
    mle2_ipcf = map(sim, mle_grouped,
                    fixed_marginal_parameters = TRUE,
                    initial_guess = pcf),
    mle4_imid = map(sim, mle_grouped,
                    fixed_marginal_parameters = FALSE),
    mle4_imid = map(sim, mle_grouped,
                    fixed_marginal_parameters = FALSE,
                    initial_guess = pcf)
  )

save(simulations_L2, file = "mle_L2.rda", version = 3, compress = "xz")

simulations_L3 <- simulations[7:9, ] %>%
  mutate(
    mle2_imid = map(sim, mle_grouped,
                    fixed_marginal_parameters = TRUE),
    mle2_ipcf = map(sim, mle_grouped,
                    fixed_marginal_parameters = TRUE,
                    initial_guess = pcf),
    mle4_imid = map(sim, mle_grouped,
                    fixed_marginal_parameters = FALSE),
    mle4_imid = map(sim, mle_grouped,
                    fixed_marginal_parameters = FALSE,
                    initial_guess = pcf)
  )

save(simulations_L3, file = "mle_L3.rda", version = 3, compress = "xz")

simulations <- bind_rows(simulations_L1, simulations_L2, simulations_L3)

save(simulations, file = "simulations_all.rda", version = 3, compress = "xz")
