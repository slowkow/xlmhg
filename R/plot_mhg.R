
plot_gsea1 <- function(fc, res, n = 500) {
  
  n <- min(length(res$mHG), n)  
  
  layout(matrix(c(1, 2), ncol = 2))
  
  # c(bottom, left, top, right) 
  par(mar = c(4, 4, 1, 0) + 0.1)
  x <- fc[1:n]
  y <- 1:n
  plot(
    x = x,
    y = y,
    ylim = rev(range(y)),
    yaxt = "n",
    ylab = "Rank",
    xlab = "Value",
    type = "l"
  )
  axis(side = 2, at = c(1, rev(y[y %% 100 == 0])), las = 2)
  abline(h = which(v[1:n] == 1), col = "blue")
  
  par(mar = c(4, 1, 1, 2) + 0.1)
  x <- -log10(res$mHG[1:n])
  y <- 1:n
  plot(
    x = x,
    y = y,
    ylim = rev(range(y)),
    yaxt = "n",
    ylab = "",
    xlab = "Enrichment score",
    type = "l"
  )
  abline(h = res$threshold, lty = 2, col = "red", lwd = 1)
}

plot_gsea2 <- function(fc, res, n = 500, main = "", ...) {
  
  n <- min(length(res$mHG), n)
  
  layout(matrix(c(3, 2, 1)), heights = c(0.45, 0.1, 0.45))
  
  # c(bottom, left, top, right)
  par(mar = c(5, 5, 0, 2) + 0.1)
  x <- 1:n
  y <- fc[x]
  plot(
    x = x,
    y = y,
    xlab = "Rank",
    ylab = "Value",
    type = "l",
    ...
  )
  
  par(mar = c(0, 5, 0, 2) + 0.1)
  plot(
    x = x,
    y = rep(0, length(x)),
    yaxt = "n",
    xaxt = "n",
    pch = NA,
    ylab = "",
    xlab = ""
  )
  abline(v = which(v[x] == 1), col = "blue")
  
  par(mar = c(0, 5, 5, 2) + 0.1)
  y <- -log10(res$mHG[x])
  plot(
    x = x,
    y = y,
    xlab = "",
    xaxt = "n",
    ylab = "Enrichment score",
    type = "l",
    main = main,
    ...
  )
  abline(v = res$threshold, lty = 2)
  text(
    x = res$threshold,
    y = quantile(y, probs = 0.05),
    adj = -0.2,
    labels = sprintf("P = %.2g", res$pvalue),
    ...
  )
}
