p_prob <- function(r, sigma, delta_t) {
  u = exp(sigma*sqrt(delta_t))
  d = (1/u)
  return((exp(r*delta_t)-d)/(u-d))
}

tree_stock <- function(S, sigma, delta_t, N) {
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

value_binomial_option <- function(tree, sigma, delta_t, r, K, type) {
  p = p_prob(r, delta_t, sigma)
  option_tree = matrix(NA, nrow=nrow(tree), ncol=ncol(tree))
  if(type == 'put') {
    option_tree[nrow(option_tree),] = pmax(K-tree[nrow(tree),],0)
  }
  else {
    option_tree[nrow(option_tree),] = pmax(tree[nrow(tree),]-K,0)
  }
  for (i in (nrow(tree)-1):1) {
    for(j in 1:i) {
      option_tree[i,j] = (p*option_tree[i+1,j+1] + (1-p)*option_tree[i+1,j])/exp(r*delta_t)
    }
  }
  return(option_tree)
}

binomial_option <- function(S, K, r, T, N, sigma, type) {
  p = p_prob(r=r, delta_t=T/N, sigma=sigma)
  tree = tree_stock(S=S, sigma=sigma, delta_t=T/N, N=N)
  option = value_binomial_option(tree, sigma=sigma, delta_t=T/N, r=r, K=K, type=type)
  delta = (option[2,2]-option[2,1])/(tree[2,2]-tree[2,1])
  return(list(p=p, delta=delta, stock=tree, option=option, price=option[1,1]))
}