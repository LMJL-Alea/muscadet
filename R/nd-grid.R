generate_nd_grid <- function(n, d, lb = rep(0, d), ub = rep(1, d)) {
  l <- list2()
  for (i in 1:d) {
    name <- paste0("k", i)
    l <- c(l, list2(!!name := 0:n))
  }
  tidyr::expand_grid(!!!l) %>%
    dplyr::mutate(
      sq_norm = colSums((t(.[1:d]) / (ub - lb))^2),
      weight = as.integer(2^rowSums(.[1:d] > 0))
    ) %>%
    dplyr::arrange(sq_norm)
}
