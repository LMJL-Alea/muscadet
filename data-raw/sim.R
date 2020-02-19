# Some simple simulated data
rho1 <- 100
rho2 <- 100
alpha1 <- 0.03
alpha2 <- 0.03
alpha12 <- 0.035
R <- 100
# Gaussian tau_max
k1 <- rho1 * pi * alpha1^2
k2 <- rho2 * pi * alpha2^2
max_k12 <- sqrt(max(min(k1 * k2, (1 - k1) * (1 - k2) - sqrt(.Machine$double.eps)), 0))
max_tau_gauss <- max_k12 / sqrt(rho1 * rho2) / (pi * alpha12^2)
# Bessel tau_max
k1 <- rho1 * (2 * pi * alpha1^2)
k2 <- rho2 * (2 * pi * alpha2^2)
max_k12 <- sqrt(max(min(k1 * k2, (1 - k1) * (1 - k2) - sqrt(.Machine$double.eps)), 0))
max_tau_bessel <- max_k12 / sqrt(rho1 * rho2) / (2 * pi * alpha12^2)
# Produce admissible tau
tau <- 0.5 # min(max_tau_gauss, max_tau_bessel) / 2

# Gaussian simulations ----------------------------------------------------

# Run simulation
sim_gauss5 <- mediator::rbidppn(n = R, seed = 1234, rho1 = rho1, rho2 = rho2, tau = tau, alpha1 = alpha1, alpha2 = alpha2, alpha12 = alpha12, progress = TRUE, Kspec = "Kspecgauss", testtau = "testtaugauss")

usethis::use_data(sim_gauss5, overwrite = TRUE, compress = "xz")

# Bessel simulations ------------------------------------------------------

# Run simulation
sim_bessel <- mediator::rbidppn(n = R, seed = 1234, rho1 = rho1, rho2 = rho2, tau = tau, alpha1 = alpha1, alpha2 = alpha2, alpha12 = alpha12, progress = TRUE, Kspec = "Kspecbessel", testtau = "testtaubessel")

usethis::use_data(sim_bessel, overwrite = TRUE, compress = "xz")
