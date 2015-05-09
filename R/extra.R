.extra <- function() {
  # Size of the population.
  N <- 5000L
  # Successes in the population.
  K <- 100L
  # Only consider enrichments in the first L observations.
  L <- N / 4L
  # Require at least X successes in the first L observations.
  X <- 5L
  
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
  
  # -----------------------------------------------------------------------------
  # Update p at each step
  
  p = 1
  
  case0 <- function(n, N, k, K) {
    a <- (n + 1) * (N - K - n + k)  
    b <- (N - n) * (n - k + 1)
    a / b
  }
  
  case1 <- function(n, N, k, K) {
    a <- (n + 1) * (K - k)
    b <- (N - n) * (k + 1)
    a / b
  }
  
  # Case 0
  p * case0(n = 1, N, k = 0, K) # -0.65
  p * case0(n = 2, N, k = 0, K) * case0(n = 1, N, k = 0, K) # 0.96
  
  # Case 1
  p * case1(n = 1, N, k = 0, K) # 0.0063
  
  # n = 0 p = 0.02
  # n = 1 p = 0.0392078
  # n = 2 p = 0.0576468
  # n = 3 p = 0.0753396
  # n = 4 p = 0.0923084
  # n = 5 p = 0.108575
  # n = 6 p = 0.124159
  # n = 7 p = 0.139083
  # n = 8 p = 0.153365
# n = 9 p = 0.167026
}