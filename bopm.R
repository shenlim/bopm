# Binomial Options Pricing Model - Oct 2018. Author: Shen Lim.

bopm <- function(s, k, r, t, n, sigma, type, div) {
  # Computes option price, option delta and risk-neutral probability based on the
  # binomial tree approach (Cox, Ross, Rubinstein 1979).
  #
  # Args:
  #   s    : Current underlying asset price.-----------------[numeric]
  #   k    : Option strike price.----------------------------[numeric]
  #   r    : Risk-free rate of interest.---------------------[numeric]*
  #   t    : Time period, in years.--------------------------[numeric]
  #   n    : Number of steps.--------------------------------[integer]**
  #   sigma: Underlying asset price volatility, in decimals.-[numeric]
  #   type : "ca" (American call)--------------------------|
  #        : "pa" (American put)---------------------------|
  #        : "ce" (European call)--------------------------|
  #        : "pe" (European put)---------------------------[character]
  #   div  : Dividend yield (optional)-----------------------[numeric]*
  #
  #   *    : If r or div is > 1, it will be assummed to be in percentages.
  #   **   : n must be an integer. If not, it will be rounded down to the nearest integer.
  #
  # Returns and Stores:
  #   1. Risk-neutral probability (p).
  #   2. Option delta (delta).
  #   3. Option price (price).
  #
  # Stores:
  #   4. Stock tree (tree_stock).
  #   5. Exercise tree (tree_exercise).
  #   6. Option tree (tree_option).
  
  if (missing(div)) {div <- 0}
  if (div > 1) {div <- div / 100}
  if (r > 1) {r <- r / 100}
  n <- as.integer(n)
  
  u <- exp(sigma * sqrt(t / n))
  d <- 1 / u
  p <<- (exp((r - div) * t / n) - d) / (u - d)
  
  tree_stock <- matrix(NA, nrow = n + 1, ncol = n + 1)
  for (j in 1:(n + 1)) {
    for (i in 1:j) {
      tree_stock[i, j] <- s * u^(j - i) * d^(i - 1)
    }
  }
  tree_stock <<- tree_stock
  
  tree_exercise <- matrix(NA, nrow = n + 1, ncol = n + 1)
  for (j in 1:(n + 1)) {
    for (i in 1:j) {
      tree_exercise[i, j] <- if (type == "ca" || type == "ce") {
        max(tree_stock[i, j] - k, 0)
      } else {
        max(k - tree_stock[i, j], 0)
      }
    }
  }
  tree_exercise <<- tree_exercise
  
  tree_option <- matrix(NA, nrow = n + 1, ncol = n + 1)
  tree_option[, n + 1] <- if (type == "ca" || type == "ce") {
    pmax(tree_stock[, n + 1] - k, 0)
  } else {
    pmax(k - tree_stock[, n + 1], 0)
  }
  for (j in n:1) {
    for (i in 1:j) {
      if (type == "ca" || type == "pa") {
        exercise_payoff <- if (type == "ca" || type == "ce") {
          max(tree_stock[i, j] - k, 0)
        } else {
          max(k - tree_stock[i, j], 0)
        }
        hold_payoff <- (p * tree_option[i, j + 1] + (1 - p) * tree_option[i + 1, j + 1]) / exp(r * t / n)
        tree_option[i, j] <- max(exercise_payoff, hold_payoff)
      } else {
        tree_option[i, j] <- (p * tree_option[i, j + 1] + (1 - p) * tree_option[i + 1, j + 1]) / exp(r * t / n)
      }
    }
  }
  tree_option <<- tree_option
  
  delta <<- (tree_option[1, 2] - tree_option[2, 2]) / (tree_stock[1, 2] - tree_stock[2, 2])
  price <<- tree_option[1, 1]
  return(list(p = p, delta = delta, price = price))
}

stocktree <- function(cex, col, pch) {
  if (missing(cex)) {cex = 1}
  if (missing(col)) {col = 2}
  if (missing(pch)) {pch = 1}
  n <- ncol(tree_stock) - 1
  plot(x = c(0, n), y = c(-n, n), type = "n", main = "Stock Tree", xlab = "Steps (n)", ylab = "", axes = FALSE)
  axis(side = 1)
  points(x = 0, y = 0, pch = pch)
  text(0, 0.3, format(round(as.numeric(deparse(tree_stock[1, 1])), digits = 2), nsmall = 2), cex = cex)
  for (i in 1:n) {
    y = seq(from = -i, by = 2, length = i + 1)
    x = rep(i, times = length(y))
    points(x, y, pch = pch)
    for(j in 1:length(x))
      text(x[j], y[j] + 0.3, format(round(as.numeric(deparse(tree_stock[i + 2 - j, i + 1])),
                                          digits = 2), nsmall = 2), cex = cex)
    y = (-i):i
    x = rep(c(i, i - 1), times = 2 * i)[1:length(y)]
    lines(x, y, col = col)
  }
}

exercisetree <- function(cex, col, pch) {
  if (missing(cex)) {cex = 1}
  if (missing(col)) {col = 2}
  if (missing(pch)) {pch = 1}
  n <- ncol(tree_exercise) - 1
  plot(x = c(0, n), y = c(-n, n), type = "n", main = "Exercise Tree", xlab = "Steps (n)", ylab = "", axes = FALSE)
  axis(side = 1)
  points(x = 0, y = 0, pch = pch)
  text(0, 0.3, format(round(as.numeric(deparse(tree_exercise[1, 1])), digits = 2), nsmall = 2), cex = cex)
  for (i in 1:n) {
    y = seq(from = -i, by = 2, length = i + 1)
    x = rep(i, times = length(y))
    points(x, y, pch = pch)
    for(j in 1:length(x))
      text(x[j], y[j] + 0.3, format(round(as.numeric(deparse(tree_exercise[i + 2 - j, i + 1])),
                                          digits = 2), nsmall = 2), cex = cex)
    y = (-i):i
    x = rep(c(i, i - 1), times = 2 * i)[1:length(y)]
    lines(x, y, col = col)
  }
}

optiontree <- function(cex, col, pch) {
  if (missing(cex)) {cex = 1}
  if (missing(col)) {col = 2}
  if (missing(pch)) {pch = 1}
  n <- ncol(tree_option) - 1
  plot(x = c(0, n), y = c(-n, n), type = "n", main = "Option Tree", xlab = "Steps (n)", ylab = "", axes = FALSE)
  axis(side = 1)
  points(x = 0, y = 0, pch = pch)
  text(0, 0.3, format(round(as.numeric(deparse(tree_option[1, 1])), digits = 2), nsmall = 2), cex = cex)
  for (i in 1:n) {
    y = seq(from = -i, by = 2, length = i + 1)
    x = rep(i, times = length(y))
    points(x, y, pch = pch)
    for(j in 1:length(x))
      text(x[j], y[j] + 0.3, format(round(as.numeric(deparse(tree_option[i + 2 - j, i + 1])),
                                          digits = 2), nsmall = 2), cex = cex)
    y = (-i):i
    x = rep(c(i, i - 1), times = 2 * i)[1:length(y)]
    lines(x, y, col = col)
  }
}