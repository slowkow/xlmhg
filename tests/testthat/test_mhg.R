context("mhg")

test_that("mhg works", {
  # Size of the population.
  N <- 5000L
  # Successes in the population.
  K <- 100L
  # Only consider enrichments in the first L observations.
  L <- N / 4L
  # Require at least X successes in the first L observations.
  X <- 5L
  
  set.seed(42)
  
  # This should be significant.
  v <- rep(0, N)
  v[sample(100, 5)] <- 1
  v[sample(200, 10)] <- 1
  
  res <- do_mHG_test(v, N, K, L, X)
  res$pvalue
  
  expect_equal(res$pvalue, 1.810658e-05, tolerance = 1e-6)
  
  # This is how you can plot the results.
#   plot_mhg(
#     fc = sort(rnorm(N, 0, 1)),
#     v = v,
#     res = res,
#     n = L,
#     main = "GO:123",
#     value = bquote("log"[2] ~ "FC"),
#     cex.lab = 2, cex.axis = 1.5, cex.main = 2, cex = 1.5
#   )
  
  # This should be non-significant.
  v <- rep(0, N)
  v[sample(N, K / 2)] <- 1

  res <- do_mHG_test(v, N, K, L, X)
  
  expect_equal(res$pvalue, 0.9902327, tolerance = 1e-6)
})
