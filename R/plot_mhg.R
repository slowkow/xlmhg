#' Plot the results of the minimum hypergeometric test.
#' @param fc
#' @param v
#' @param res
#' @param n
#' @param main
#' @param cex.main
#' @param style
#' @param ...
plot_mhg <- function(fc, v, res, n = 500, main = "", cex.main = 1, style = "h", ...) {
  if (style == "h") {
    .plot_mhg_horizontal(
      fc, v, res, n, main, cex.main = cex.main, ...)
  } else if (style == "v") {
    .plot_mhg_vertical(
      fc, v, res, n, main, cex.main = cex.main, ...)
  } else {
    stop(sprintf("style must be 'v' or 'h'"))
  }
}

.plot_mhg_vertical <- function(fc, v, res, n = 500, main = "", cex.main = 1,
                               ...) {
  
  n <- min(length(res$mHG), n)  
  
  layout(matrix(c(1, 2, 3), ncol = 3), widths = c(0.45, 0.1, 0.45))
  
  # c(bottom, left, top, right) 
  par(mar = c(5, 5, 3, 0) + 0.1)
  x <- fc[1:n]
  y <- 1:n
  plot(
    x = x,
    y = y,
    ylim = rev(range(y)),
    yaxt = "n",
    ylab = "Rank",
    xlab = "Value",
    type = "l",
    ...
  )
  axis(side = 2, at = c(1, rev(y[y %% 100 == 0])), las = 2, ...)
  
  par(mar = c(5, 0, 3, 0) + 0.1)
  plot(
    x = rep(0, length(y)),
    y = y,
    yaxt = "n",
    xaxt = "n",
    ylim = rev(range(y)),
    pch = NA,
    ylab = "",
    xlab = ""
  )
  abline(h = which(v[y] == 1), col = "blue")
  
  par(mar = c(5, 0, 3, 2) + 0.1)
  x <- -log10(res$mHG[y])
  plot(
    x = x,
    y = y,
    ylim = rev(range(y)),
    yaxt = "n",
    ylab = "",
    xlab = "Enrichment score",
    type = "l",
    ...
  )
  abline(h = res$threshold, lty = 2)
  text(
    x = quantile(x, probs = 0.05),
    y = res$threshold,
    adj = c(0, 2),
    labels = sprintf("P = %.2g", res$pvalue),
    ...
  )
  
  title(main, outer = TRUE, line = -2, cex.main = cex.main)
}

.plot_mhg_horizontal <- function(fc, v, res, n = 500, main = "", ...) {
  
  n <- min(length(res$mHG), n)
  
  layout(matrix(c(3, 2, 1)), heights = c(0.45, 0.1, 0.45))
  
  # c(bottom, left, top, right)
  par(mar = c(5, 6, 0, 2) + 0.1,
      mgp = c(3, 1, 0.2))
  x <- 1:n
  y <- fc[x]
#   plot(
#     x = x,
#     y = y,
#     xlab = "Rank",
#     ylab = "Value",
#     type = "l",
#     ...
#   )
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
  mtext(side = 2, text = "Value", line = 4, ...)
  axis(side = 1, at = c(1, x[x %% 100 == 0]), ...)
  
  par(mar = c(0, 6, 0, 2) + 0.1)
  plot(
    x = x,
    y = rep(0, length(x)),
    yaxt = "n",
    xaxt = "n",
    pch = NA,
    ylab = "",
    xlab = "",
    xaxs = "i"
  )
  abline(v = which(v[x] == 1), col = adjustcolor("blue", alpha.f = 0.6))
  
  par(mar = c(0, 6, 5, 2) + 0.1)
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
    y = quantile(y, probs = 0.05),
    adj = -0.2,
    labels = sprintf("P = %.2g", res$pvalue),
    ...
  )
}
