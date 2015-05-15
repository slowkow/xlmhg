#' Plot the results of the minimum hypergeometric test.
#' 
#' @param values The values corresponding to the ranked vector.
#' @param x A binary vector of successes and failures.
#' @param res A list returned by mhg_test.
#' @param n Plot the first n items in the list.
#' @param main The title of the plot.
#' @param value The name of the value used to rank the items.
#' @param cex The relative size of the legend text.
#' @param cex.lab The relative size of the axis labels.
#' @param cex.axis The relative size of the axis tickmarks.
#' @param cex.main The relative size of the title.
#'
#' @seealso \code{\link{mhg_test}}
#'
#' @examples
#' # Size of the population.
#' N <- 5000L
#' # Successes in the population.
#' K <- 100L
#' # Only consider enrichments in the first L observations.
#' L <- N / 4L
#' # Require at least X successes in the first L observations.
#' X <- 5L
#'   
#' set.seed(42)
#'   
#' # Binary vector of successes and failures.
#' x <- rep(0, N)
#' x[sample(100, 5)] <- 1
#' x[sample(200, 10)] <- 1
#'   
#' res <- mhg_test(x, N, K, L, X)
#'
#' abs(res$pvalue - 1.810658e-05) < 1e-6 # TRUE
#'
#' # Plot the result.
#' plot_mhg(sort(rnorm(N)), x, res, L)
plot_mhg <- function(
  values, x, res, n = 500, main = "", value = "Value",
  cex = 1.5, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2) {
  
  n <- min(length(res$mhg), n)
  
  layout(matrix(c(3, 2, 1)), heights = c(0.45, 0.1, 0.45))
  
  # c(bottom, left, top, right)
  par(mar = c(5, 7, 0, 2) + 0.1,
      mgp = c(3, 1, 0.2))
  xs <- 1:n
  y <- values[xs]
  barplot(
    height = y,
    col = "lightgrey",
    border = NA,
    space = 0,
    names.arg = NA,
    xlab = "Rank",
    ylab = "",
    xaxs = "i",
    las = 2,
    cex.lab = cex.lab * 1.4,
    cex.axis = cex.axis,
    cex = cex
  )
  mtext(side = 2, text = value, line = 4, cex = cex)
  axis(side = 1, at = c(1, xs[xs %% 100 == 0]), cex.axis = cex.axis)
  
  par(mar = c(0, 7, 0, 2) + 0.1)
  plot(
    x = xs,
    y = rep(0, length(xs)),
    yaxt = "n",
    xaxt = "n",
    pch = NA,
    ylab = "",
    xlab = "",
    xaxs = "i",
    bty = "n"
  )
  abline(v = which(x[xs] == 1), col = adjustcolor("blue", alpha.f = 0.6))
  
  par(mar = c(0, 7, 5, 2) + 0.1)
  # Remove NAs.
  y <- -log10( sapply(res$mhg, function(i) ifelse(is.na(i), 1, i)) )
  plot(
    x = xs,
    y = y[xs],
    xlab = "",
    xaxt = "n",
    xaxs = "i",
    ylab = "",
    ylim = c(0, max(y)),
    bty = "n",
    type = "l",
    col = "darkblue",
    main = main,
    las = 2,
    cex.axis = cex.axis,
    cex.main = cex.main
  )
  mtext(side = 2, text = "Enrichment", line = 4, cex = cex)
  abline(v = res$threshold, lty = 2, lwd = 2, col = "red")
  legend(
    x = res$threshold,
    y = 0.25 * (max(y) - min(y)), #quantile(y, probs = 0.99),
    legend = sprintf("P = %.2g", res$pvalue),
    bg = adjustcolor("white", alpha.f = 0.8),
    box.lty = 0,
    cex = cex * 1.3
  )
#   text(
#     x = res$threshold,
#     y = 0.25 * (max(y) - min(y)), #quantile(y, probs = 0.99),
#     adj = -0.2,
#     labels = sprintf("P = %.2g", res$pvalue),
#     ...
#   )
}
