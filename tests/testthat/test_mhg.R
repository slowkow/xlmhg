context("mhg")

test_that("mhg works", {
  # Size of the population.
  N <- 5000L
  # Successes in the population.
  K <- 100L
  # Only consider enrichments in the first L observations.
  # L <- N / 4L
  L <- 400
  # L <- N
  # Require at least X successes in the first L observations.
  X <- 5L
  # ???
  mat <- matrix(1)
  
  set.seed(420)
  v <- rep(0, N)
  v[sample(100, 5)] <- 1
  v[sample(200, 15)] <- 1
#   v[sample(N, 50)] <- 1
  
  res <- do_mHG_test(v, N, K, L, X, mat)
  res$pvalue
  
  fc <- sort(rnorm(N, 0, 1))
  
  plot_mhg(fc, v, res, L, "GO:123",
           cex.lab = 2, cex.axis = 1.5, cex.main = 2, cex = 1.5, style = 'h')
  
  # # #
  
  # The hypergeometric distribution.
  dhyper(x = 1, m = K, n = N - K, k = 5)
  
  score <- function(k, K, n, N) {
    sum(sapply(k:min(n, K), function(i) {
      a <- choose(n, i) * choose(N - n, K - i)
      b <- choose(N, K)
      a / b
    }))
  }
  N <- 5000
  K <- 100
  sapply(c(1, 10, 100), function(n) round(-log10( score(k = n, K, n, N)) ))
  sapply(c(1, 10, 100), function(n) round(-log10( score(k = 0, K, n, N)) ))
  
  k <- 5
  n <- 25
  # Probability of k successes after n trials.
  p <- dhyper(x = k, m = K, n = N - K, k = n)
  
  florian_hg <- function(p, k, K, n, N) {
    pval <- p
    for (i in k:min(n, K)) {
      p <- p * (n - i) * (K - i) / ((i + 1) * (N - K - n + i + 1))
      pval <- pval + p
    }
    pval
  }
  score(k, K, n, N)
  florian_hg(p, k, K, n, N)

  # expect_equal(grab(xs, regex = "\\/", offset = 2), as.character(6:10))
})
