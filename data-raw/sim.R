# Some simple simulated data
rho1 <- 100
rho2 <- 100
alpha1 <- 0.03
alpha2 <- 0.03
alpha12 <- 0.05
tau <- 0.2

# Test parameter feasability
k1 <- rho1 * (2*pi*alpha1^2)
k2 <- rho2 * (2*pi*alpha2^2)
tau^2 * alpha12^4
alpha1^2 * alpha2^2 * min(1, (1/k1 - 1) * (1/k2 - 1))

# Run simulation
R <- 100
sim <- list()
for(i in 1:R){
  cat(i)
  sim[[i]] <- simdppxmatern(rho1, rho2, alpha1, alpha2, alpha12, tau, progress = 1, Kspec = "Kspecbessel", testtau = "testtaubessel")
}

usethis::use_data(sim)
