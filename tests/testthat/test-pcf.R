test_that("estimate() works", {
  set.seed(1234)
  df <- rbidpp(
    n = 1,
    rho1 = 100, rho2 = 100,
    alpha1 = 0.03, alpha2 = 0.03,
    alpha12 = 0.035, tau = 0.5,
    L = 1, model = "Gauss"
  )
  res <- estimate(df[[1]])
  expect_snapshot(res)
})
