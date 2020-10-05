generate_nd_grid <- function(n, d) {
  l <- list2()
  for (i in 1:d) {
    name <- paste0("k", i)
    l <- c(l, list2(!!name := 0:n))
  }
  tidyr::expand_grid(!!!l) %>%
    dplyr::mutate(
      sq_norm = rowSums(.[1:d]^2),
      weight = as.integer(2^rowSums(.[1:d] > 0))
    ) %>%
    dplyr::arrange(sq_norm)
}
