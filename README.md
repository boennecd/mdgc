
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mdgc

TODO: The goal of mdgc is to …

## Example

``` r
library(mdgc)

# simulates a data set and puts into a format to pass to C++
#
# Args: 
#   n: number of observations. 
#   p: number of variables. 
# 
# Returns: 
#   Simluated mask data set in the format needed to pass to C++.
sim_dat <- function(n, p = 4){
  # get the covariance matrix
  Sb <- diag(p)
  Sb[lower.tri(Sb)] <- Sb[upper.tri(Sb)] <- .5
  Sb <- Sb / p / 5
  Sig <- cov2cor(drop(rWishart(1L, 5L * p, Sb)))
  
  # draw the observations
  Z <- truth <- crossprod(chol(Sig), matrix(rnorm(n * p), p))
  
  # mask 
  is_mask <- matrix(runif(n * p) < .33, p)
  is_int <- ceiling(p / 3):p
  is_mask[is_int, ] <- is_mask[is_int, ] & Z[is_int, ] < 0
  
  Z[ is_int, ][is_mask[ is_int, ]] <- 0
  Z[-is_int, ][is_mask[-is_int, ]] <- NA_real_
  
  # create matrix in the Z format to pass to c++
  lower <- matrix(-Inf, p, n)
  upper <- Z
  # codes are: 
  #  0: latent Z is observed (upper is the observed point).
  #  1: latent Z can be anything.. 
  #  2: latent Z is in an interval. 
  code <- matrix(0L, p, n)
  code[-is_int, ][is_mask[-is_int, ]] <- 1L 
  code[ is_int, ][is_mask[ is_int, ]] <- 2L 
  
  list(lower = lower, upper = upper, code = code, Sigma = Sig, 
       truth = truth)
}

set.seed(2)
p <- 4L
dat <- sim_dat(2000L, p = p)
dat$lower[, 1:10]
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,] -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf  -Inf
#> [2,] -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf  -Inf
#> [3,] -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf  -Inf
#> [4,] -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf  -Inf
dat$upper[, 1:10]
#>        [,1]   [,2]  [,3]   [,4]  [,5]   [,6]   [,7]   [,8]    [,9]  [,10]
#> [1,] -1.034 -0.976 -1.03 -0.388 0.885     NA     NA  0.183 -0.4573     NA
#> [2,] -0.312  0.809  0.00  0.739 2.442  0.000 -0.999  0.479  0.0000 -1.506
#> [3,]  0.000  0.133  0.00 -2.121 0.962 -0.283  0.453 -0.153 -0.4099 -0.488
#> [4,]  0.000 -0.782  0.47  0.936 1.259  0.000  1.150  0.991  0.0437  0.000
dat$code [, 1:10]
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,]    0    0    0    0    0    1    1    0    0     1
#> [2,]    0    0    2    0    0    2    0    0    2     0
#> [3,]    2    0    2    0    0    0    0    0    0     0
#> [4,]    2    0    0    0    0    2    0    0    0     2
dat$truth[, 1:10]
#>        [,1]   [,2]   [,3]   [,4]  [,5]   [,6]   [,7]   [,8]    [,9]  [,10]
#> [1,] -1.034 -0.976 -1.034 -0.388 0.885 -0.642 -0.495  0.183 -0.4573 -0.853
#> [2,] -0.312  0.809 -0.797  0.739 2.442 -0.328 -0.999  0.479 -0.5222 -1.506
#> [3,] -0.638  0.133 -1.269 -2.121 0.962 -0.283  0.453 -0.153 -0.4099 -0.488
#> [4,] -1.341 -0.782  0.470  0.936 1.259 -1.586  1.150  0.991  0.0437 -0.273

# Get pointers to objects in C++
ptr <- mdgc:::get_log_lm_terms(lower = dat$lower, upper = dat$upper, 
                               code = dat$code)

# define log marginal likelihood function
log_ml <- function(vcov, releps = 1e-3, n_threads = 1L)
  mdgc:::eval_log_lm_terms(
    ptr = ptr, vcov = vcov, indices = 0:(NCOL(dat$lower) - 1L), 
    maxpts = 100000L, abseps = -1, releps = releps, n_threads = n_threads)

# print the marginal likelihood at the true parameters
set.seed(1)
print(log_ml(dat$Sigma), digits = 7)
#> [1] -9609.76
print(log_ml(dat$Sigma, n_threads = 4L), digits = 7)
#> [1] -9609.76
sd(replicate(20, log_ml(dat$Sigma, n_threads = 4L)))
#> [1] 0.00143

# check computation time
library(microbenchmark)
microbenchmark(
  `1 thread`     = log_ml(dat$Sigma), 
  `2 threads`    = log_ml(dat$Sigma, n_threads = 2L),
  `4 threads`    = log_ml(dat$Sigma, n_threads = 4L),
  `lower releps` = log_ml(dat$Sigma, releps = 1e-4), times = 10)
#> Unit: milliseconds
#>          expr   min    lq  mean median    uq   max neval
#>      1 thread  76.8  77.3  78.2   77.7  79.0  81.1    10
#>     2 threads  39.8  39.8  40.2   40.3  40.5  40.9    10
#>     4 threads  21.2  21.5  21.7   21.8  21.8  22.0    10
#>  lower releps 185.0 186.8 188.3  188.1 190.1 192.4    10

# estimate the parameters
par_fn <- function(par){
  set.seed(1)
  # use log-cholesky parametrization
  log_sds <- par[1:p]
  L <- diag(exp(log_sds), p)
  L[lower.tri(L)] <- par[-(1:p)]
  Arg <- tcrossprod(L)
  
  -log_ml(cov2cor(Arg), n_threads = 4L)
}

start_val <- c(rep(log(1), p), numeric(p * (p - 1) / 2))
system.time(opt_out <- optim(start_val, par_fn))
#>    user  system elapsed 
#>   34.69    0.00    8.67

# check the result
local({
  par <- opt_out$par
  log_sds <- par[1:p]
  L <- diag(exp(log_sds), p)
  L[lower.tri(L)] <- par[-(1:p)]
  cov2cor(tcrossprod(L))
})
#>       [,1]    [,2]  [,3]    [,4]
#> [1,] 1.000  0.2512 0.272  0.2122
#> [2,] 0.251  1.0000 0.167 -0.0083
#> [3,] 0.272  0.1665 1.000  0.1763
#> [4,] 0.212 -0.0083 0.176  1.0000
dat$Sigma # the truth 
#>       [,1]    [,2]  [,3]    [,4]
#> [1,] 1.000  0.2753 0.276  0.2045
#> [2,] 0.275  1.0000 0.155 -0.0168
#> [3,] 0.276  0.1546 1.000  0.1989
#> [4,] 0.205 -0.0168 0.199  1.0000
```
