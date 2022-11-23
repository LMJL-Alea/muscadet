test_that("fit_via_pcf() works", {
  res <- fit_via_pcf(sim_gauss0[[1]])
  expect_snapshot(res)
})
