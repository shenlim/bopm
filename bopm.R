# Binomial Options Pricing Model - Oct 2018. Author: Shen Lim

bopm <- function(s, k, r, t, n, sigma, type, style, div) {
  # Computes option price based on binomial tree (discrete-time) approach (Cox, Ross, Rubinstein 1979).
  #
  # Args:
  #   s:     Current underlying asset price.---------------------[numeric]
  #   k:     Option strike price.--------------------------------[numeric]
  #   r:     Risk-free rate of interest expressed in decimals.---[numeric]
  #   t:     Time period expressed in years.---------------------[numeric]
  #   n:     Number of steps.------------------------------------[integer]
  #   sigma: Underlying asset price volatility.------------------[numeric]
  #   type:  "call" or "put".------------------------------------[character]
  #   style: "e" for European options, "a" for American options.-[character]
  #   div:   Dividend yield expressed in decimals; "0" if none.--[numeric]
  #
  # Returns and stores:
  #   (1) Risk neutral probability, (2) Stock tree, (3) Exercise payoff tree,
  #   (4) Option tree, (5) Option delta, (6) Option price at t = 0.
  
  u <- exp(sigma * sqrt(t / n))
  d <- 1 / u
  p <<- (exp((r - div) * t / n) - d) / (u - d)
  
  tree_stock <- matrix(NA, nrow = n + 1, ncol = n + 1)
  for (i in 1:(n + 1)) {
    for (j in 1:i) {
      tree_stock[i, j] <- s * u^(j - 1) * d^(i - j)
    }
  }
  tree_stock <<- tree_stock
  
  tree_exercise <- matrix(NA, nrow = n + 1, ncol = n + 1)
  for (i in 1:(n + 1)) {
    for (j in 1:i) {
      tree_exercise[i, j] <- if (type == "call") max(tree_stock[i, j] - k, 0) else max(k - tree_stock[i, j], 0)
    }
  }
  tree_exercise <<- tree_exercise
  
  tree_option <- matrix(NA, nrow = n + 1, ncol = n + 1)
  tree_option[n + 1, ] <- if (type == "call") pmax(tree_stock[n + 1, ] - k, 0) else pmax(k - tree_stock[n + 1, ], 0)
  for (i in n:1) {
    for (j in 1:i) {
      if (style == "a") {
        exercise_payoff <- if (type == "call") max(tree_stock[i, j] - k, 0) else max(k - tree_stock[i, j], 0)
        hold_payoff <- (p * tree_option[i + 1, j + 1] + (1 - p) * tree_option[i + 1, j]) / exp(r * t / n)
        tree_option[i, j] <- max(exercise_payoff, hold_payoff)
      } else {
        tree_option[i, j] <- (p * tree_option[i + 1, j + 1] + (1 - p) * tree_option[i + 1, j]) / exp(r * t / n)
      }
    }
  }
  tree_option <<- tree_option
  
  delta <<- (tree_option[2, 2] - tree_option[2, 1]) / (tree_stock[2, 2] - tree_stock[2, 1])
  price <<- tree_option[1, 1]
  return(list(p = p, stock_tree = tree_stock, exercise_tree = tree_exercise, option_tree = tree_option, delta = delta, price = price))
}