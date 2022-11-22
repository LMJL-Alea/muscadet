test_that("get_feasible_repulsion_rate() works", {
  expect_snapshot(get_feasible_repulsion_rate(rho = 100))
})

test_that("get_feasible_space() works", {
  expect_snapshot(get_feasible_space(
    rho1 = 100, rho2 = 100,
    alpha1 = 0.03, alpha2 = 0.03
  ))
})
