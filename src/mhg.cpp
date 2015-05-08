#include <RcppArmadillo.h>
using namespace Rcpp;

//' Return the minimum of two numbers.
//' 
//' @param a Numeric.
//' @param b Numeric.
//' @return The smaller of a or b.
long double min(long double a, long double b) {
  if (a < b) return a;
  return b;
}

//' Return the maximum of two numbers.
//' 
//' @param a Numeric.
//' @param b Numeric.
//' @return The larger of a or b.
long double max(long double a, long double b) {
  if (a > b) return a;
  return b;
}

//' Test for equality between numbers with 80 bit extended precision.
//' 
//' @param a Numeric.
//' @param b Numeric.
//' @param tol Numeric.
//' @return Integer 0 or 1.
int is_equal(long double a, long double b, long double tol) {
  if (a == b) {
    return 1;
  } else if (fabsl(a - b) / max(fabsl(a), fabsl(b)) < tol) {
    return 1; 
  }
  return 0;
}

//' Starting with the tail probability of finding some particular number of
//' successes less than k, update it to the tail probability of finding k or
//' more successes.
//' 
//' This is function is only useful in the context of this package, where we
//' have split the probability computation across several functions.
//' 
//' @param p Probability of exactly k successes.
//' @param k Number of successes after n draws.
//' @param N Size of the population.
//' @param K Number of successes in the population.
//' @param n Number of draws, without replacement.
//' @return Numeric between 0 and 1.
// [[Rcpp::export]]
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

//' Add two numbers in log space.
//' @param log_a
//' @param log_b
//' @return The sum in log space: log(exp(log_a) + exp(log_b))
// [[Rcpp::export]]
double log_add(double log_a, double log_b) {
  if (log_a < log_b) {
    std::swap(log_a, log_b);
  }
  const double negdelta = log_b - log_a;
  if (arma::arma_isfinite(negdelta) == false) {
    return log_a;
  }
  return log_a + log1p(exp(negdelta));
}

double get_hypergeometric_pvalue2(double logp, int k, int N, int K, int n) {
  long double pval = logp;
  for (int i = k; i < min(K, n); i++) {
    logp = logp + log(
      (long double)( (n - i) * (K - i) ) /
      (long double)( (i + 1) * (N - K - n + i + 1) )
    );
    pval = log_add(pval, logp);
  }
  return pval;
}

//' Compute the minimum hypergeometric score (mHG). It indicates the
//' probability of the observed density of ones at the beginning of the vector
//' under the assumption that all possible permutations of the list are equally
//' probable.
//' 
//' The optimal size of the beginning is discovered by iteratively increasing
//' the size and chooseing the one with the minimum hypergeometric score.
//' 
//' @param v Binary vector of ones and zeros.
//' @param N Size of the population.
//' @param K Number of successes in the population.
//' @param L Only consider scores for the first L observations.
//' @param X Require at least X ones to get a score less than 1.
//' @param scores A vector of mHG scores.
//' @param tol The tolerance for testing equality of two numbers.
int get_mHG(arma::vec v, int N, int K, int L, int X,
  arma::vec &scores, long double tol) {
  
  if (K == 0 || K == N || K < X) {
    scores[0] = 1.0;
    return 0;
  }
  
  long double p = 1.0;
  long double pval;
  long double mHG = 1.1;
  int k = 0;
  int threshold = 0;
  
  for (int n = 0; n < L; n++) {
    // We see a zero in the presence vector.
    if (v[n] == 0) {
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
//      // Compute a p-value if we have seen enough elements.
//      if (k >= X) {
//        pval = get_hypergeometric_pvalue(p, k, N, K, n + 1);
//        
//        if (pval < mHG && !is_equal(pval, mHG, tol)) {
//          mHG = pval;
//          threshold = n + 1;
//        }
//        scores.at(n) = mHG;
//      }
    }
//    if (scores.at(n) == 0) {
//      scores.at(n) = p;
//    }
    pval = get_hypergeometric_pvalue(p, k, N, K, n + 1);
    if (k >= X) {
      if (pval < mHG && !is_equal(pval, mHG, tol)) {
        mHG = pval;
        threshold = n + 1;
      }
    }
    scores.at(n) = pval;
  }
  
  // We did not see enough ones in the first L elements of v.
  if (threshold == 0) {
    scores.at(0) = 1.0;
  }
//  else {
//    scores.at(0) = mHG;
//  }
  
  return threshold;
}

//' Compute a minimum hypergeometric (mHG) p-value by dynamic programming.
long double get_mHG_pvalue(int N, int K, int L, int X,
  long double mHG, arma::mat matrix, long double tol) {
    
  if (mHG > 1.0 || is_equal(mHG, 1.0, tol)) {
    return 1.0;
  } else if (mHG == 0) {
    return 0;
  } else if (K == 0 || K >= N || K < X) {
    return 0;
  } else if (L > N) {
    return 0;
  }
  
  int W = N - K;
  int n, k, w;
  
  long double p_start = 1.0;
  long double p;
  long double pval;
  matrix.at(0, 0) = 1.0;
  
  for (n = 1; n < N; n++) {
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
    w = n - k;
    
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

// [[Rcpp::export]]
Rcpp::List do_mHG_test(arma::vec v, int N, int K, int L, int X, arma::mat mat,
  bool use_upper_bound = false, bool verbose = false,
  long double tolerance = 0.0000000000000001) {
  
  if (N < 0) stop("Condition not met: N >= 0");
  if (K < 0 || K > N) stop("Condition not met: 0 <= K <= N");
  if (L < 0 || L > N) stop("Condition not met: 0 <= L <= N");
  if (X < 0 || X > K) stop("Condition not met: 0 <= X <= K");

  if (K == 0 || K == N) {
    return Rcpp::List::create(
      Rcpp::Named("threshold") = 0,
      Rcpp::Named("mHG") = 1.0,
      Rcpp::Named("pvalue") = 1.0
    );
  }
  
  arma::mat matrix = arma::mat(K + 1, N - K + 1).fill(0);
  
  if (mat.n_rows >= K + 1 && mat.n_cols >= N - K + 1) {
    matrix = mat;
  }
  
  int threshold;
  long double mHG, mHG_pvalue;
  double mHG_double, mHG_pvalue_double;
  
  // Get XL-mHG scores and the threshold for the best score.
  arma::vec scores = arma::vec(L).fill(0);
  threshold = get_mHG(v, N, K, L, X, scores, tolerance);
  // mHG = scores[0];
  mHG = scores.at(threshold);
  
  // There is no enrichment at all.
  if (is_equal(mHG, 1.0, tolerance)) {
    return Rcpp::List::create(
      Rcpp::Named("threshold") = 0,
      Rcpp::Named("mHG") = 1.0,
      Rcpp::Named("pvalue") = 1.0
    );
  }
  
  // Don't compute XL-mHG p-value, instead use upper bound.
  if (use_upper_bound) {
    mHG_pvalue_double = (double)min(1.0, mHG * (long double)(K));
  }
  // Compute an XL-mHG p-value.
  else {
    mHG_pvalue = get_mHG_pvalue(N, K, L, X, mHG, matrix, tolerance);
    mHG_pvalue_double = (double)(mHG_pvalue);
    
    // Check if floating point accuracy was insufficient for the p-value.
    if (mHG_pvalue_double <= 0) {
      // Use the upper bound instead.
      mHG_pvalue_double = (double)(mHG * (long double)(K));
    }
  }
  
  mHG_double = (double)(mHG);

  return Rcpp::List::create(
    Rcpp::Named("threshold") = threshold + 1,
    // Rcpp::Named("mHG") = mHG_double,
    Rcpp::Named("mHG") = wrap(scores),
    Rcpp::Named("pvalue") = mHG_pvalue_double
  );
}
