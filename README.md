# Binomial Options Pricing Model
Building generalized pricing models for options in discrete-time. The model was first proposed by [Cox, Ross, and Rubinstein](https://www.sciencedirect.com/science/article/pii/0304405X79900151?via%3Dihub) in 1979. This `README.md` will briefly explain the functions in the R-Script `bopm.R` and how they match with BOPM methods. Please refer to `bopm.pdf` for a more theoretical explanation of the model.

## R-Script
First, under the risk neutrality assumption, define `p`:

<a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{120}&space;p=\frac{e^{rT}-d}{u-d}" target="_blank"><img src="https://latex.codecogs.com/png.latex?\dpi{120}&space;p=\frac{e^{rT}-d}{u-d}" title="p=\frac{e^{rT}-d}{u-d}" /></a>

where

<a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{120}&space;u=e^{\sigma&space;\sqrt{\Delta{T}}}&space;\text{&space;;&space;}&space;d=e^{-\sigma&space;\sqrt{\Delta{T}}}" target="_blank"><img src="https://latex.codecogs.com/png.latex?\dpi{120}&space;u=e^{\sigma&space;\sqrt{\Delta{T}}}&space;\text{&space;;&space;}&space;d=e^{-\sigma&space;\sqrt{\Delta{T}}}" title="u=e^{\sigma \sqrt{\Delta{T}}} \text{ ; } d=e^{-\sigma \sqrt{\Delta{T}}}" /></a>

so that:
```r
p_prob <- function(r, delta_t, sigma, div) {
  u = exp(sigma*sqrt(delta_t))
  d = (1/u)
  return((exp((r-div)*delta_t)-d)/(u-d))
}
```
With `u` and `d`, construct the binomial price tree:
```r
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
```
The value of the stock (underlying asset) at each node can be calculated as:

<a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{120}&space;S_n=S_0\cdot{u^{N_u}}\cdot{d^{N_d}}&space;\text{&space;;&space;}&space;N_u=\text{number&space;of&space;up&space;mov.&space;;&space;}&space;N_d&space;=&space;\text{number&space;of&space;down&space;mov.}" target="_blank"><img src="https://latex.codecogs.com/png.latex?\dpi{120}&space;S_n=S_0\cdot{u^{N_u}}\cdot{d^{N_d}}&space;\text{&space;;&space;}&space;N_u=\text{number&space;of&space;up&space;mov.&space;;&space;}&space;N_d&space;=&space;\text{number&space;of&space;down&space;mov.}" title="S_n=S_0\cdot{u^{N_u}}\cdot{d^{N_d}} \text{ ; } N_u=\text{number of up mov. ; } N_d = \text{number of down mov.}" /></a>

Generalizing this to the matrix (tree) results in `tree[i,j] = S*u^(j-1)*d^(i-j)`.

Next, build the option payoff (value) tree. Note that the option value is `max[(Sn-K),0]` for a call option and `max[(K-Sn),0]` for a put option; where `K` is the strike price and `Sn` is the spot price of the stock at the *n*th period.
```r
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
```
The option price can be expressed as:

<a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{120}&space;f=[f_{u}(p)&plus;f_{d}(1-p)]e^{-rT}" target="_blank"><img src="https://latex.codecogs.com/png.latex?\dpi{120}&space;f=[f_{u}(p)&plus;f_{d}(1-p)]e^{-rT}" title="f=[f_{u}(p)+f_{d}(1-p)]e^{-rT}" /></a>

resulting in `option_tree[i,j] = (p*option_tree[i+1,j+1] + (1-p)*option_tree[i+1,j])/exp(r*delta_t)`. Note that for an American option (`american = TRUE` or `american = 1`), we must compare the exercise payoff and the option value at each node and choose whichever is higher.

Putting it together:
```r
binomial_option <- function(S, K, r, T, N, sigma, type, american, div) {
  p = p_prob(r= r, delta_t= T/N, sigma= sigma, div= div)
  tree = tree_stock(S= S, N= N, delta_t= T/N, sigma= sigma)
  ex = exercise(S= S, K= K, N= N, delta_t= T/N, sigma= sigma, type= type)
  option = value_binomial_option(K= K, r= r, delta_t= T/N, sigma= sigma, type= type, american, tree, ex, div)
  delta = (option[2,2]-option[2,1])/(tree[2,2]-tree[2,1])
  return(list(p= p, delta= delta, stock= tree, exercise= ex, option= option, price= option[1,1]))
}
```
Recall that

<a href="https://www.codecogs.com/eqnedit.php?latex=\dpi{120}&space;\Delta=\frac{f_{u}-f_{d}}{S_{0}u-S_{0}d}" target="_blank"><img src="https://latex.codecogs.com/png.latex?\dpi{120}&space;\Delta=\frac{f_{u}-f_{d}}{S_{0}u-S_{0}d}" title="\Delta=\frac{f_{u}-f_{d}}{S_{0}u-S_{0}d}" /></a>

resulting in `delta = (option[2,2]-option[2,1])/(tree[2,2]-tree[2,1])`.

## Generating the Output
The option function is given by:
```r
binomial_option <- function(S, K, r, T, N, sigma, type, american, div)
```
Inputs:

* `S`: Current underlying asset price.
* `K`: Option strike price.
* `r`: Risk-free rate of interest; convert percantages to decimals.
* `T`: Time period; expressed in years.
* `N`: Number of steps.
* `sigma`: Underlying asset price volatility.
* `type`: `call` for call option, `put` for put option.
* `american`: `0` or `FALSE` for European option, `1` or `TRUE` for American option.
* `div`: Dividend yield, if any; convert percentages to decimals.

## Built With
* [R](https://www.r-project.org/) - Software environment for statistical computing.
* [RStudio](https://www.rstudio.com/) - Free and open-source IDE for R.
* [TeX](https://www.latex-project.org/get/)
* [CodeCogs](https://www.codecogs.com/latex/eqneditor.php) - Online equation editor.

## License
This project is licensed under the MIT License - see the `LICENSE.md` file for details.
