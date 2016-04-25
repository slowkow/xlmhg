#' Non-parametric rank enrichment test for binary data.
#'
#' This package implements a non-parametric enrichment test for
#' ranked binary lists. This is a general package suitable for many
#' applications. For example, it can be used for gene-set enrichment tests to
#' find enrichment for gene annotations.
#' 
#' Run the test with \code{\link{mhg_test}}
#'
#' Plot the results with \code{\link{plot_mhg}}
#' 
#' @references
#' Eden, E., Lipson, D., Yogev, S. & Yakhini, Z. Discovering motifs in ranked
#' lists of DNA sequences. PLoS Comput. Biol. 3, e39 (2007).
#' \url{http://dx.doi.org/10.1371/journal.pcbi.0030039}
#'
#' Wagner, F. GO-PCA: An Unsupervised Method to Explore Biological
#' Heterogeneity Based on Gene Expression and Prior Knowledge. bioRxiv (2015).
#' \url{http://dx.doi.org/10.1101/018705}
#'
#' @docType package
#' @name mhg
NULL

# Test a vector for enrichment of large values.
#' @param y A named vector of sorted values (e.g. F statistic)
#' @param lists A list with sets of names present in names(y).
#' @param L Consider the first L ranks.
#' @param X Consider enrichments with at least X items in the first L ranks.
#' @param workers Run the test on multiple cores.
do_mhg <- function(y, lists, L = NULL, X = 5, workers = 4) {
  N <- length(y)
  if (is.null(L)) {
    L <- N / 2
  }
  
  lengths <- lapply(lists, function(xs) sum(xs %in% names(y)))
  ylists <- lists[lengths >= X]
  
  # ylists <- lapply(lists, function(xs) xs[xs %in% names(y)])
  # ylists <- ylists[sapply(ylists, length) >= X]
  
  multicoreParam <- BiocParallel::MulticoreParam(workers = workers)
  
  retval <- bplapply(BPPARAM = multicoreParam, X = ylists, FUN = function(i) {
    K <- length(i)
    
    #     v_down <- names(sort(y, decreasing = FALSE)) %in% i
    # v_up <- names(sort(y, decreasing = TRUE)) %in% i
    v_up <- names(y) %in% i
    #     v_both <- names(sort(abs(y), decreasing = TRUE)) %in% i
    
    #     res_down <- mhg_test(v_down, N, K, L, X)
    res_up <- mhg_test(v_up, N, K, L, X)
    #     res_both <- mhg_test(v_both, N, K, L, X)
    
    #     if (res_down$threshold == L) res_down$pvalue <- 1
    if (res_up$threshold == L) res_up$pvalue <- 1
    #     if (res_both$threshold == L) res_both$pvalue <- 1
    
    list(
      "n" = length(i),
      #       "n.down" = sum(v_down[1:res_down$threshold]),
      #       "p.down" = res_down$pvalue,
      "n.up" = sum(v_up[1:res_up$threshold]),
      "p.up" = res_up$pvalue
      #       "n.abs" = sum(v_both[1:res_both$threshold]),
      #       "p.abs" = res_both$pvalue
    )
  })
  #   retval
  rnames <- names(retval)
  retval <- data.frame(matrix(
    unlist(retval), nrow = length(retval), byrow = TRUE
  ))
  rownames(retval) <- rnames
  #   colnames(retval) <- c(
  #     "n", "n.down", "p.down", "n.up", "p.up", "n.abs", "p.abs")
  colnames(retval) <- c("n", "n.up", "p.up")
  
  #   retval$p.down[retval$p.down == 0] <- 1
  retval$p.up[retval$p.up == 0] <- 1
  #   retval$p.abs[retval$p.abs == 0] <- 1
  #   retval$q.down <- qvalue(retval$p.down)$qvalues
  #   retval$q.up <- qvalue(retval$p.up)$qvalues
  #   retval$q.both <- qvalue(retval$p.both)$qvalues
  
  retval[order(retval[["p.up"]]), ]  
}