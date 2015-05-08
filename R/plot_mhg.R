#' Plot the results of the minimum hypergeometric test.
#' @param fc
#' @param v
#' @param res
#' @param n
#' @param main
#' @param cex.main
#' @param style
#' @param ...
plot_mhg <- function(fc, v, res, n = 500, main = "", value = "Value", ...) {
  
  n <- min(length(res$mHG), n)
  
  layout(matrix(c(3, 2, 1)), heights = c(0.45, 0.1, 0.45))
  
  # c(bottom, left, top, right)
  par(mar = c(5, 7, 0, 2) + 0.1,
      mgp = c(3, 1, 0.2))
  x <- 1:n
  y <- fc[x]
  barplot(
    height = y,
    col = "lightgrey",
    border = NA,
    space = 0,
    xlab = "Rank",
    ylab = "",
    xaxs = "i",
    las = 2,
    ...
  )
  mtext(side = 2, text = value, line = 4, ...)
  axis(side = 1, at = c(1, x[x %% 100 == 0]), ...)
  
  par(mar = c(0, 7, 0, 2) + 0.1)
  plot(
    x = x,
    y = rep(0, length(x)),
    yaxt = "n",
    xaxt = "n",
    pch = NA,
    ylab = "",
    xlab = "",
    xaxs = "i",
    bty = "n"
  )
  abline(v = which(v[x] == 1), col = adjustcolor("blue", alpha.f = 0.6))
  
  par(mar = c(0, 7, 5, 2) + 0.1)
  y <- -log10(res$mHG[x])
  plot(
    x = x,
    y = y,
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
    ...
  )
  mtext(side = 2, text = "Enrichment score", line = 4, ...)
  abline(v = res$threshold, lty = 2, lwd = 2, col = "red")
  text(
    x = res$threshold,
    y = 0.25 * (max(y) - min(y)), #quantile(y, probs = 0.99),
    adj = -0.2,
    labels = sprintf("P = %.2g", res$pvalue),
    ...
  )
}
