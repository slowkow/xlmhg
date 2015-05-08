context("mhg")

test_that("mhg works", {
  # Size of the population.
  N <- 5000L
  # Successes in the population.
  K <- 100L
  # Only consider enrichments in the first L observations.
  L <- N / 8L
  # Require at least X successes in the first L observations.
  X <- 5L
  # ???
  mat <- matrix(1)
  
  set.seed(420)
  v <- rep(0, N)
  v[sample(100, 5)] <- 1
  v[sample(200, 5)] <- 1
  
  res <- do_mHG_test(v, N, K, L, X, mat)
  
  fc <- sort(rnorm(N, 0, 1))
  
  # expect_equal(grab(xs, regex = "\\/", offset = 2), as.character(6:10))
  
  # plot_gsea1(fc, res, n = 500)
  
  plot_gsea2(fc, res, n = 500, main = "GO:123", cex.axis = 2, cex.lab = 2, cex.main = 3, cex = 1.5)
})
