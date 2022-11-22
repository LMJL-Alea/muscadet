generate_nd_grid <- function(n, d, lb = rep(0, d), ub = rep(1, d)) {
  l <- rlang::list2()
  for (i in 1:d) {
    name <- paste0("k", i)
    l <- c(l, rlang::list2(!!name := 0:n))
  }
  res <- tidyr::expand_grid(!!!l)
  res$sq_norm <- colSums((t(res[1:d]) / (ub - lb))^2)
  res$weight <- as.integer(2^rowSums(res[1:d] > 0))
  dplyr::arrange(res, .data$sq_norm)
}
