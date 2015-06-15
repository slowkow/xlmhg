/**
 * Minimum-hypergeometric test for enrichment in ranked binary lists.
 * Copyright (C) 2015  Kamil Slowikowski <kslowikowski@fas.harvard.edu>
 *               2015  Florian Wagner <florian.wagner@duke.edu>
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <RcppArmadillo.h>
using namespace Rcpp;

// Return the minimum of two numbers.
// 
// @param a Numeric.
// @param b Numeric.
// @return The smaller of a or b.
long double min(long double a, long double b) {
  if (a < b) return a;
  return b;
}

// Return the maximum of two numbers.
// 
// @param a Numeric.
// @param b Numeric.
// @return The larger of a or b.
long double max(long double a, long double b) {
  if (a > b) return a;
  return b;
}

// Test for equality between numbers with 80 bit extended precision.
// 
// @param a Numeric.
// @param b Numeric.
// @param tol Numeric.
// @return Integer 0 or 1.
int is_equal(long double a, long double b, long double tol) {
  if (a == b) {
    return 1;
  } else if (fabsl(a - b) / max(fabsl(a), fabsl(b)) < tol) {
    return 1; 
  }
  return 0;
}

// Starting with the tail probability of finding some particular number of
// successes less than k, update it to the tail probability of finding k or
// more successes.
// 
// This is function is only useful in the context of this package, where we
// have split the probability computation across several functions.
// 
// @param p Probability of exactly k successes.
// @param k Number of successes after n draws.
// @param N Size of the population.
// @param K Number of successes in the population.
// @param n Number of draws, without replacement.
// @return Numeric between 0 and 1.
long double get_hypergeometric_pvalue(long double p, int k, int N, int K, int n) {
  long double pval = p;
  for (int i = k; i < min(K, n); i++) {
    p = p * (long double)( (n - i) * (K - i) ) /
        (long double)( (i + 1) * (N - K - n + i + 1) );
    pval = pval + p;
  }
  if (pval < 0) {
    pval = 0;
  }
  return pval;
}

// Compute the minimum hypergeometric score (mHG). It indicates the
// probability of the observed density of ones at the beginning of the vector
// under the assumption that all possible permutations of the list are equally
// probable.
// 
// The optimal size of the beginning is discovered by iteratively increasing
// the size and chooseing the one with the minimum hypergeometric score.
// 
// @param x Binary vector of ones and zeros.
// @param N Size of the population.
// @param K Number of successes in the population.
// @param L Only consider scores for the first L observations.
// @param X Require at least X ones to get a score less than 1.
// @param scores A vector of mHG scores.
// @param tol The tolerance for testing equality of two numbers.
int get_mHG(arma::vec x, int N, int K, int L, int X,
  arma::vec &scores, long double tol) {
  
  long double p = 1.0;
  long double pval;
  long double mHG = 1.1;
  int k = 0;
  int threshold = 0;
  
  if (K == 0 || K == N || K < X) {
    scores.at(0) = 1.0;
    return threshold;
  }
  
  for (int n = 0; n < L; n++) {
    // We see a zero in the presence vector.
    if (x[n] == 0) {
      // Compute P(k | N,K,n+1) from P(k | N,K,n)
      p = p * (long double)( (n + 1) * (N - K - n + k) ) /
          (long double)( (N - n) * (n - k + 1) );
    }
    // We see a one in the presence vector.
    else {
      // Compute P(k+1 | N,K,n+1) from P(k | N,K,n)
      p = p * (long double)( (n + 1) * (K - k) ) /
    			(long double)( (N - n) * (k + 1) );
			k = k + 1;
    }
    
    // Get the probability for k or more successes after n + 1 trials.
    pval = get_hypergeometric_pvalue(p, k, N, K, n + 1);
    scores.at(n) = pval;
    
    // Update the best score we have seen, and the corresponding threshold.
    if (x[n] != 0 && k >= X) {
      if (pval < mHG && !is_equal(pval, mHG, tol)) {
        mHG = pval;
        threshold = n + 1;
      }
    }
  }
  
  // We did not see enough ones in the first L elements of x.
  if (threshold == 0) {
    scores.at(0) = 1.0;
  }
  
  return threshold;
}

// Compute a minimum hypergeometric (mHG) p-value by dynamic programming.
// 
// @param N Size of the population.
// @param K Number of successes in the population.
// @param L Only consider scores for the first L observations.
// @param X Require at least X ones to get a score less than 1.
// @param mHG Find how many mHG scores are better than this observation.
// @param tol The tolerance for testing equality of two numbers.
long double get_mHG_pvalue(
  int N, int K, int L, int X, long double mHG, long double tol) {
    
  if (mHG > 1.0 || is_equal(mHG, 1.0, tol)) {
    return 1.0;
  } else if (mHG == 0) {
    return 0;
  } else if (K == 0 || K >= N || K < X) {
    return 0;
  } else if (L > N) {
    return 0;
  }
  
  // Number of failures in the population.
  int W = N - K;
  // Number of: trials, successes, failures.
  int n, k, w;
  
  long double p_start = 1.0;
  long double p;
  long double pval;
  
  arma::mat matrix = arma::mat(K + 1, N - K + 1).fill(0.0);
  matrix.at(0, 0) = 1.0;
  
  for (n = 1; n < N; n++) {
    // Number of sucesses in the population is at least number of trials.
    if (K >= n) {
      k = n;
      p_start = p_start * (long double)(K - n + 1) /
                (long double)(N - n + 1);
    } else {
      k = K;
      p_start = p_start * (long double)(n) /
                (long double)(n - K);
    }
    
    // We lack enough floating point precision to compute the p-value.
    if (p_start <= 0) {
      return NA_REAL;
    }
    
    p = p_start;
    pval = p_start;
    // Number of failures.
    w = n - k;
    
    // This trial is within the first L elements, and is at least X.
    if (n <= L && n >= X) {

      while (k >= X && w < W && (is_equal(pval, mHG, tol) || pval < mHG)) {
        matrix.at(k, w) = 0;
        p = p * (long double)( k * (N - K - n + k) ) /
            (long double)( (n - k + 1) * (K - k + 1) );
        pval = pval + p;
        w = w + 1;
        k = k - 1;
      }
    }
    
    while (k >= 0 && w <= W) {
      if (w > 0 && k > 0) {

        matrix.at(k, w) = matrix.at(k, w - 1) *
          (long double)(W - w + 1) / (long double)(N - n + 1) +
          matrix.at(k - 1, w) *
          (long double)(K - k + 1) / (long double)(N - n + 1);

      } else if (w > 0) {

        matrix.at(k, w) = matrix.at(k, w - 1) *
          (long double)(W - w + 1) / (long double)(N - n + 1);

      } else if (k > 0) {

        matrix.at(k, w) = matrix.at(k - 1, w) *
          (long double)(K - k + 1) / (long double)(N - n + 1);

      }
      w = w + 1;
      k = k - 1;
    }
  }
  
  return 1.0 - (matrix.at(K, W - 1) + matrix.at(K - 1, W));
}

//' Test for enrichment in a ranked binary list.
//' 
//' @description
//' Given a ranked binary list of ones and zeros, test if the ones are
//' enriched at the beginning of the list.
//'
//' @details
//' Suppose we have a set of \code{N = 5000} genes and \code{K = 100} of them
//' are annotated with a Gene Ontology (GO) term. Further, suppose that we find
//' some subset of these genes to be significantly differentially expressed
//' (DE) between two conditions. Within the DE genes, we notice that
//' \code{k = 15} of the DE genes are annotated with the Gene Ontology term. At
//' this point, we would like to know if the GO term is enriched for DE genes.
//'
//' 
//' We use the hypergeometric distribution to compute a probability that we
//' would observe a given number of DE genes annotated with a GO term. You
//' can find more details in the documentation for \code{\link{dhyper}}.
//' 
//' The method consists of three steps:
//' 
//' \itemize{
//' \item Compute a hypergeometric probability at each rank in the list.
//' \item Choose the minimum hypergeometric probability (mHG) as the test
//'       statistic.
//' \item Use dynamic programming to compute the exact permutation p-value
//'       for observing a test statistic at least as extreme by chance.
//' }
//' 
//' @param x Binary vector of ones and zeros.
//' @param N Size of the population.
//' @param K Number of successes in the population.
//' @param L Only consider scores for the first L observations.
//' @param X Require at least X ones to get a score less than 1.
//' @param upper_bound Instead of running a dynamic programming algorithm,
//'   return the upper bound for the p-value.
//' @param tol The tolerance for testing equality of two numbers.
//'
//' @return A list with items "threshold", "mHG", and "pvalue".
//'
//' @seealso \code{\link{plot_mhg}}
//'
//' @examples
//' # Size of the population.
//' N <- 5000L
//' # Successes in the population.
//' K <- 100L
//' # Only consider enrichments in the first L observations.
//' L <- N / 4L
//' # Require at least X successes in the first L observations.
//' X <- 5L
//'
//' set.seed(42)
//'
//' # Binary vector of successes and failures.
//' x <- rep(0, N)
//' x[sample(100, 5)] <- 1
//' x[sample(200, 10)] <- 1
//'
//' res <- mhg_test(x, N, K, L, X)
//'
//' abs(res$pvalue - 1.810658e-05) < 1e-6 # TRUE
//'
//' # Plot the result.
//' plot_mhg(sort(rnorm(N)), x, res, L)
//'
//' @references
//' Eden, E., Lipson, D., Yogev, S. & Yakhini, Z. Discovering motifs in ranked
//' lists of DNA sequences. PLoS Comput. Biol. 3, e39 (2007).
//' \url{http://dx.doi.org/10.1371/journal.pcbi.0030039}
//'
//' Wagner, F. GO-PCA: An Unsupervised Method to Explore Biological
//' Heterogeneity Based on Gene Expression and Prior Knowledge. bioRxiv (2015).
//' \url{http://dx.doi.org/10.1101/018705}
// [[Rcpp::export]]
Rcpp::List mhg_test(arma::vec x, int N, int K, int L, int X,
  bool upper_bound = false, long double tol = 0.0000000000000001) {
  
  if (N < 0) stop("Condition not met: N >= 0");
  if (K < 0 || K > N) stop("Condition not met: 0 <= K <= N");
  if (L < 0 || L > N) stop("Condition not met: 0 <= L <= N");
  if (X < 0 || X > K) stop("Condition not met: 0 <= X <= K");

  if (K == 0 || K == N) {
    return Rcpp::List::create(
      Rcpp::Named("threshold") = 0,
      Rcpp::Named("mhg") = 1.0,
      Rcpp::Named("pvalue") = 1.0
    );
  }
  
  int threshold;
  long double mHG, mHG_pvalue;
  double mHG_pvalue_double;
  
  // Get XL-mHG scores and the threshold for the best score.
  arma::vec scores = arma::vec(L).fill(0);
  threshold = get_mHG(x, N, K, L, X, scores, tol);
  mHG = scores.at(threshold);
  
  // There is no enrichment at all.
  if (is_equal(mHG, 1.0, tol)) {
    return Rcpp::List::create(
      Rcpp::Named("threshold") = 0,
      Rcpp::Named("mhg") = 1.0,
      Rcpp::Named("pvalue") = 1.0
    );
  }
  
  // Don't compute XL-mHG p-value, instead use upper bound.
  if (upper_bound) {
    mHG_pvalue_double = (double)min(1.0, mHG * (long double)(K));
  }
  // Compute an XL-mHG p-value.
  else {
    mHG_pvalue = get_mHG_pvalue(N, K, L, X, mHG, tol);
    mHG_pvalue_double = (double)(mHG_pvalue);
    
    // Check if floating point accuracy was insufficient for the p-value.
    if (mHG_pvalue_double <= 0) {
      // Use the upper bound instead.
      mHG_pvalue_double = (double)(mHG * (long double)(K));
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("threshold") = threshold,
    Rcpp::Named("mhg") = wrap(scores),
    Rcpp::Named("pvalue") = mHG_pvalue_double
  );
}
