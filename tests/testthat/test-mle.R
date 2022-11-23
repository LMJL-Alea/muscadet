test_that("fit_via_mle() works", {
  res <- fit_via_mle(sim_gauss0[[1]])
  expect_snapshot(res$value)
})
