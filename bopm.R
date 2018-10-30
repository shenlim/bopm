# Binomial Options Pricing Model - Oct 2018. Author: Shen Lim
# Copyright 2018, Shen Lim, All Rights Reserved.

p_prob <- function(r, delta_t, sigma, div) {
  u = exp(sigma*sqrt(delta_t))
  d = (1/u)
  return((exp((r-div)*delta_t)-d)/(u-d))
}

tree_stock <- function(S, N, delta_t, sigma) {
  tree = matrix(NA, nrow=N+1, ncol=N+1)
  u = exp(sigma*sqrt(delta_t))
  d = (1/u)

  for (i in 1:(N+1)) {
    for (j in 1:i) {
      tree[i,j] = S*u^(j-1)*d^(i-j)
    }
  }
  return(tree)
}

exercise <- function(S, K, N, delta_t, sigma, type) {
  ex = matrix(NA, nrow=N+1, ncol=N+1)
  u = exp(sigma*sqrt(delta_t))
  d = (1/u)

  if(type == "call") {
    for (i in 1:(N+1)) {
      for (j in 1:i) {
        ex[i,j] = max(((S*u^(j-1)*d^(i-j))-K),0)
      }
    }
    return(ex)
  }
  else {
    for (i in 1:(N+1)) {
      for (j in 1:i) {
        ex[i,j] = max((K-(S*u^(j-1)*d^(i-j))),0)
      }
    }
    return(ex)
  }
}

value_binomial_option <- function(K, r, delta_t, sigma, type, american, tree, ex, div) {
  p = p_prob(r, delta_t, sigma, div)
  option_tree = matrix(NA, nrow=nrow(tree), ncol=ncol(tree))
  if(type == "call") {
    option_tree[nrow(option_tree),] = pmax(tree[nrow(tree),]-K,0)
  }
  else {
    option_tree[nrow(option_tree),] = pmax(K-tree[nrow(tree),],0)
  }
  if (american == TRUE) {
    for (i in (nrow(tree)-1):1) {
      for(j in 1:i) {
        exercise.payoff <- if (type == "call") max(tree[i,j] - K, 0) else max(K - tree[i,j], 0)
        hold.payoff <- (p*option_tree[i+1,j+1] + (1-p)*option_tree[i+1,j])/exp(r*delta_t)
        option_tree[i,j] <- max(exercise.payoff, hold.payoff)
      }
    }
    return(option_tree)
  }
  else {
    for (i in (nrow(tree)-1):1) {
      for(j in 1:i) {
        option_tree[i,j] = (p*option_tree[i+1,j+1] + (1-p)*option_tree[i+1,j])/exp(r*delta_t)
      }
    }
    return(option_tree)
  }
}

binomial_option <- function(S, K, r, T, N, sigma, type, american, div) {
  p = p_prob(r= r, delta_t= T/N, sigma= sigma, div= div)
  tree = tree_stock(S= S, N= N, delta_t= T/N, sigma= sigma)
  ex = exercise(S= S, K= K, N= N, delta_t= T/N, sigma= sigma, type= type)
  option = value_binomial_option(K= K, r= r, delta_t= T/N, sigma= sigma, type= type, american, tree, ex, div)
  delta = (option[2,2]-option[2,1])/(tree[2,2]-tree[2,1])
  return(list(p= p, delta= delta, stock= tree, exercise= ex, option= option, price= option[1,1]))
}