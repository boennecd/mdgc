
# mdgc

This package contains a marginal likelihood approach to estimating the
model discussed by D. Hoff (2007) and Zhao and Udell (2019). We have
modified the Fortran code by Genz and Bretz (2002) to supply an
approximation gradient for the log marginal likelihood and to use an
approximation of the marginal likelihood similar to the CDF
approximation in Genz and Bretz (2002).

However, the package can be useful for a lot of other models. For
instance, the methods are directly applicable to other Gaussian copula
models and some mixed effect models. All methods are implemented in C++,
support computation in parallel, and should easily be able to be ported
to other languages.

## Example

Below, we provide an example similar to Zhao and Udell (2019 Section
7.1). The authors use a data set with a random correlation matrix, 5
continuous variables, 5 binary variables, and 5 ordinal variables with 5
levels. There is a total of 2000 observations and 30% of the variables
are missing completely at random.

To summarize Zhao and Udell (2019) results, they show that their
approximate EM algorithm converges in what seems to be 20-25 seconds
while it takes more than 150 seconds for the MCMC algorithm used by D.
Hoff (2007). In essence, these figures should be kept in mind when
looking at the results below. Importantly, Zhao and Udell (2019) use an
approximation in the E-step of the EM algorithm which is fast but might
be crude is some settings. Using a potentially arbitrarily precise
approximation of the marginal likelihood is useful if this can be done
quickly enough.

To summarize, we do the following below:

1.  simulate the data set we will use.
2.  show how to use the C++ functions and that these provide an
    approximation of the log marginal likelihood and its gradient.
    Moreover, we show that the methods scales well in the number of
    threads.
3.  define functions to perform maximum likelihood estimation.
4.  estimate the parameters using a simple gradient descent algorithm
    and stochastic gradient descent using ADAM.
5.  show how to improve 4. by using better starting values which are
    quick to compute. As of this writing, this reduces the estimation
    time to about 4 seconds using four threads and about 12 seconds
    using one thread.

Presumably/likely, computing the marginals of each variable should be
extremely fast, as should the imputation once the model parameters are
estimated. Given that this is true, then the main concern should be the
time it takes to estimate the model parameters. As we show, this can be
done quickly using our code. We start of with the example using the high
end API. Then we turn to a more detailed explanation.

### Quick Example

We first simulate a data set and provide an example which shows how to
use the package.

``` r
# load the packages we need
library(bench)
library(mdgc)
```

``` r
# simulates a data set and mask some of the data.
#
# Args: 
#   n: number of observations. 
#   p: number of variables. 
#   n_lvls: number of levels for the ordinal variables. 
# 
# Returns: 
#   Simluated masked data and true covariance matrix. 
sim_dat <- function(n, p = 4, n_lvls = 5L){
  # get the covariance matrix
  Sb <- diag(p)
  Sb[lower.tri(Sb)] <- Sb[upper.tri(Sb)] <- .5
  Sb <- Sb / p / 5
  Sig <- cov2cor(drop(rWishart(1L, 5L * p, Sb)))
  
  # draw the observations
  truth <- matrix(rnorm(n * p), n) %*% chol(Sig)
  
  # determine the type
  type <- rep(1:3, each = floor((p + 3 - 1) / 3))[1:p]
  is_con <- type == 1L
  is_bin <- type == 2L
  is_ord <- type == 3L
  
  # sample which are masked data 
  is_mask <- matrix(runif(n * p) < .3, n)
  
  # create observed data
  truth_obs <- data.frame(truth)
  truth_obs[, is_con] <- qexp(pnorm(as.matrix(truth_obs[, is_con])))
  
  bs_bin <- c(-Inf, 0., Inf)
  truth_obs[, is_bin] <- truth_obs[, is_bin] > bs_bin[2]
  
  bs_ord <- qnorm(seq(0, 1, length.out = n_lvls + 1L))
  truth_obs[, is_ord] <- as.integer(cut(truth[, is_ord], breaks = bs_ord))
  for(i in which(is_ord)){
    truth_obs[, i] <- ordered(truth_obs[, i])
    levels(truth_obs[, i]) <- 
      LETTERS[seq_len(length(unique(truth_obs[, i])))]
  }

  # mask the data
  seen_obs <- truth_obs
  seen_obs[is_mask] <- NA
  
  list(truth = truth, truth_obs = truth_obs, seen_obs = seen_obs, 
       Sigma = Sig)
}

# simulate and show the data
set.seed(1)
p <- 15L
dat <- sim_dat(2000L, p = p)

# how an observed data set could look
head(dat$seen_obs)
#>      X1   X2    X3    X4    X5    X6    X7    X8    X9   X10  X11  X12  X13 X14
#> 1 0.560   NA    NA 0.812 0.800  TRUE    NA  TRUE    NA    NA    C    D    C   C
#> 2    NA 1.85 0.132 0.215    NA    NA  TRUE    NA    NA FALSE <NA>    A <NA>   A
#> 3 1.435   NA    NA 0.575 0.891  TRUE  TRUE  TRUE    NA  TRUE    A <NA>    B   C
#> 4 0.636   NA 0.455 0.227 1.727  TRUE FALSE  TRUE FALSE  TRUE    C    B <NA>   B
#> 5 0.664   NA 1.334    NA    NA  TRUE  TRUE FALSE    NA    NA    E <NA>    D   D
#> 6 0.285   NA 0.309 0.178 0.156 FALSE    NA FALSE    NA    NA <NA>    A    B   A
#>    X15
#> 1    D
#> 2    B
#> 3 <NA>
#> 4 <NA>
#> 5 <NA>
#> 6    B

# assign objects needed for model estimation
mdgc_obj <- get_mgdc(dat$seen_obs)
log_ml_ptr <- get_mgdc_log_ml(mdgc_obj)
start_val <- mgdc_start_value(mdgc_obj)

# this is very fast so we can neglect this when we consider the computation 
# time
mark(`Setup time` = {
  mdgc_obj <- get_mgdc(dat$seen_obs)
  log_ml_ptr <- get_mgdc_log_ml(mdgc_obj)
  start_val <- mgdc_start_value(mdgc_obj)
}, min_iterations = 10)
#> # A tibble: 1 x 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 Setup time   10.7ms     11ms      87.8    7.33MB     17.1

# fit the model using two different methods
set.seed(60941821)
system.time(
  fit_adam <- mgdc_fit(
    ptr = log_ml_ptr, vcov = start_val, n_threads = 4L, 
    lr = 1e-2, maxit = 5L, batch_size = 100L, method = "adam"))
#>    user  system elapsed 
#>  15.105   0.003   4.052
system.time(
  fit_svrg <- mgdc_fit(
    ptr = log_ml_ptr, vcov = start_val, n_threads = 4L, 
    lr = 1e-3, maxit = 5L, batch_size = 100L, method = "svrg"))
#>    user  system elapsed 
#>  31.619   0.003   8.564

# compare the log marginal likelihood 
mgdc_log_ml(vcov = fit_adam$result, ptr = log_ml_ptr, releps = 1e-3)
#> [1] -21703
mgdc_log_ml(vcov = fit_svrg$result, ptr = log_ml_ptr, releps = 1e-3)
#> [1] -21701

# compare the estimated correlation matrix with the true value
do_plot <- function(est, truth, main){
  par_old <- par(mfcol = c(1, 3), mar  = c(1, 1, 4, 1))
  on.exit(par(par_old))
  sc <- colorRampPalette(c("Red", "White", "Blue"))(51)
  
  f <- function(x, main)
    image(x[, NCOL(x):1], main = main, col = sc, zlim = c(-1, 1), 
          xaxt = "n", yaxt = "n", bty = "n")
  f(est, main)
  f(truth, "Truth")
  f(est - truth, "Difference")
}

do_plot(fit_adam$result, dat$Sigma, "Estimates (ADAM)")
```

<img src="man/figures/README-sim_dat-1.png" width="100%" />

``` r
do_plot(fit_svrg$result, dat$Sigma, "Estimates (SVRG)")
```

<img src="man/figures/README-sim_dat-2.png" width="100%" />

### Detailed Example

In this section, we show a more detailed example where we use some of
the non-exported functions package. This section is mainly included to
give an idea of what is going on under the hood.

``` r
# assign function to evalute the log marginal likelihood
log_ml <- function(...)
  mgdc_log_ml(ptr = log_ml_ptr, ...)

# print the approximate log marginal likelihood at the true parameters
set.seed(1)
print(log_ml(dat$Sigma, n_threads = 1L), digits = 7)
#> [1] -21748.55

# check standard error
sd(replicate(20, log_ml(dat$Sigma, n_threads = 4L)))
#> [1] 0.0244

# without reordering
print(log_ml(dat$Sigma, n_threads = 4L, do_reorder = FALSE), digits = 7)
#> [1] -21748.53

# check standard error
sd(replicate(20, log_ml(dat$Sigma, n_threads = 4L, do_reorder = FALSE)))
#> [1] 0.0644

# check computation time
mark(
  `1 thread                 ` = 
    log_ml(dat$Sigma                    , n_threads = 1L), 
  `1 thread  (w/o rordering)` = 
    log_ml(dat$Sigma, do_reorder = FALSE, n_threads = 1L), 
  `2 threads                ` = 
    log_ml(dat$Sigma                    , n_threads = 2L),
  `2 threads (w/o rordering)` = 
    log_ml(dat$Sigma, do_reorder = FALSE, n_threads = 2L),
  `4 threads                ` = 
    log_ml(dat$Sigma                    , n_threads = 4L), 
  `4 threads (w/o rordering)` = 
    log_ml(dat$Sigma, do_reorder = FALSE, n_threads = 4L), 
  min_iterations = 5, check = FALSE)
#> # A tibble: 6 x 6
#>   expression                     min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 1 thread                     494ms    495ms      2.01    33.9KB        0
#> 2 1 thread  (w/o rordering)    760ms    771ms      1.29    33.9KB        0
#> 3 2 threads                    260ms    265ms      3.79    33.9KB        0
#> 4 2 threads (w/o rordering)    398ms    421ms      2.40    33.9KB        0
#> 5 4 threads                    135ms    138ms      7.16    33.9KB        0
#> 6 4 threads (w/o rordering)    225ms    230ms      4.37    33.9KB        0

#####
# we can also get an approximation of the gradient
t1 <- log_ml(dat$Sigma, comp_derivs = TRUE, n_threads = 1L, releps = 1e-3)
t2 <- log_ml(dat$Sigma, comp_derivs = TRUE, n_threads = 4L, releps = 1e-3)
all.equal(t1, t2, tolerance = 1e-2)
#> [1] TRUE
  
mark(
  `1 thread                 ` = 
    log_ml(dat$Sigma, comp_derivs = TRUE                    , n_threads = 1L), 
  `1 thread  (w/o rordering)` = 
    log_ml(dat$Sigma, comp_derivs = TRUE, do_reorder = FALSE, n_threads = 1L), 
  `2 threads                ` = 
    log_ml(dat$Sigma, comp_derivs = TRUE                    , n_threads = 2L),
  `2 threads (w/o rordering)` = 
    log_ml(dat$Sigma, comp_derivs = TRUE, do_reorder = FALSE, n_threads = 2L),
  `4 threads                ` = 
    log_ml(dat$Sigma, comp_derivs = TRUE                    , n_threads = 4L), 
  `4 threads (w/o rordering)` = 
    log_ml(dat$Sigma, comp_derivs = TRUE, do_reorder = FALSE, n_threads = 4L), 
  min_iterations = 5, check = FALSE)
#> # A tibble: 6 x 6
#>   expression                     min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 1 thread                     2.35s    2.38s     0.419    35.7KB        0
#> 2 1 thread  (w/o rordering)    4.13s    4.28s     0.234    35.7KB        0
#> 3 2 threads                    1.22s    1.26s     0.784    35.7KB        0
#> 4 2 threads (w/o rordering)    2.13s    2.17s     0.456    35.7KB        0
#> 5 4 threads                 664.58ms 685.15ms     1.45     35.7KB        0
#> 6 4 threads (w/o rordering)    1.16s    1.23s     0.822    35.7KB        0

#####
# the main code in the packages provides an approximation to the CDF similar 
# to the one in the mvtnorm package. We provide an example below to 
# illustrate this. Feel free to skip this part of the README.
library(mvtnorm)
set.seed(1)
p_ex <- 5L
S_ex <- diag(p_ex)
S_ex[lower.tri(S_ex)] <- S_ex[upper.tri(S_ex)] <- .25
m_ex <- seq(-2, 2, length.out = p_ex)
lower_ex <- m_ex + drop(rnorm(p_ex) %*% chol(S_ex)) - 1
upper_ex <- lower_ex + 1

use_mvtnorm <- function()
  pmvnorm(
    lower = lower_ex, upper = upper_ex, sigma = S_ex, mean = m_ex, 
    algorithm = GenzBretz(maxpts = 100000L, abseps = -1, releps = 1e-5))
use_this_pkg <- function(derivs = FALSE)
  mdgc:::pmvnorm(lower = lower_ex, upper = upper_ex, mu = m_ex, 
                 Sigma = S_ex, maxvls = 100000L, abseps = -1, releps = 1e-5, 
                 derivs = derivs)
use_mvtnorm()
#> [1] 0.00136
#> attr(,"error")
#> [1] 1.09e-08
#> attr(,"msg")
#> [1] "Normal Completion"
use_this_pkg()
#>         [,1]
#> [1,] 0.00136
#> attr(,"minvls")
#> [1] 6992
#> attr(,"inform")
#> [1] 0
#> attr(,"abserr")
#> [1] 1.12e-08
all.equal(c(use_mvtnorm()), c(use_this_pkg()), tolerance = 1e-5)
#> [1] TRUE
mark(mvtnorm = use_mvtnorm(), mdgc = use_this_pkg(), 
     min_iterations = 25, check = FALSE)
#> # A tibble: 2 x 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 mvtnorm      1.07ms   3.61ms      258.    4.43KB     2.01
#> 2 mdgc         3.02ms   7.59ms      134.    2.49KB     0

sd(replicate(25, use_mvtnorm()))
#> [1] 3.29e-09
sd(replicate(25, use_this_pkg()))
#> [1] 2.96e-09

# the latter function can also provide gradients with respect to the mean 
# and covariance matrix
library(numDeriv)
gr_hat <- jacobian(function(a){
  m <- a[1:p_ex]
  S <- matrix(nr = p_ex, nc = p_ex)
  S[upper.tri(S, TRUE)] <- a[-(1:p_ex)]
  S[lower.tri(S)] <- t(S)[lower.tri(S)]
  
  set.seed(1L)
  res <- pmvnorm(
    lower = lower_ex, upper = upper_ex, sigma = S, mean = m, 
    algorithm = GenzBretz(maxpts = 10000L, abseps = -1, releps = 1e-6))
  c(res)
}, c(m_ex, S_ex[upper.tri(S_ex, TRUE)]))
gr <- use_this_pkg(TRUE)

# the off diagonal elements of the covariance matrix are not scaled by 2
gr_hat / gr[-1]
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
#> [1,]    1    1    1    1    1    1    2    1    2     2     1     2     2     2
#>      [,15] [,16] [,17] [,18] [,19] [,20]
#> [1,]     1     2     2     2     2     1

# creates a matrix from a log-Cholesky decomposition.
# 
# Args:
#   par: p (p + 1) / 2 elements in the log-Cholesky decomposition.
get_lchol_inv <- function(par){
  # use log-cholesky parametrization
  p <- (sqrt(8 * length(par) + 1) - 1) / 2
  L <- matrix(0, p, p)
  L[lower.tri(L, TRUE)] <- par
  diag(L) <- exp(diag(L))
  tcrossprod(L)
}

# creates the log-Cholesky decomposition. 
# 
# Args: 
#   par: positive definite matrix to decompose
get_lchol <- function(par){
  lSig <- t(chol(par))
  diag(lSig) <- log(diag(lSig))
  lSig[lower.tri(lSig, TRUE)]
}

# indices used to apply a matrix product with a get_commutation matrix
com_vec <- mdgc:::get_commutation_vec(p, p, FALSE)

# computes the approximate log marginal likelihood. 
#
# Args:
#   par: log-Cholesky decomposition.
#   seed: seed to use.  
#   comp_derivs: logical for whether to approximate the gradient. 
#   n_threads: number of threads. 
#   releps: relative error for each term.
#   indices: integer vector with which terms to include. 
par_fn <- function(par, seed = NULL, comp_derivs = FALSE, 
                   n_threads = 1L, releps = 1e-2, 
                   indices = 0:(NROW(dat$seen_obs) - 1L)){
  if(!is.null(seed))
    set.seed(seed)
  Arg <- get_lchol_inv(par)
  
  res <- log_ml(Arg, comp_derivs = comp_derivs, indices = indices,
                n_threads = n_threads, releps = releps)
  log_ml <- c(res)
  if(comp_derivs){
    gr <- attr(res, "grad")
    tmp <- matrix(0, p, p)
    tmp[lower.tri(tmp, TRUE)] <- par
    diag(tmp) <- exp(diag(tmp))
    gr <- gr[com_vec] + c(gr)
    gr <- mdgc:::x_dot_X_kron_I(x = gr, X = tmp, l = p)
    gr <- gr[, lower.tri(tmp, TRUE)]
    idx_diag <- c(1L, 1L + cumsum(NCOL(tmp):2)) 
    gr[idx_diag] <- gr[idx_diag] * diag(tmp)
      
    attr(log_ml, "grad") <- gr
    
  }
  
  log_ml
}

# check that the function gives the correct log marginal likelihood
# approximation and gradient approximation.
lSig <- get_lchol(dat$Sigma)
r1 <- par_fn(lSig, comp_derivs = TRUE, n_threads = 4L, releps = 1e-3, 
             indices = 1:100)
r2 <- numDeriv::jacobian(par_fn, lSig, seed = 1L, n_threads = 6L, 
                         releps = 1e-3, indices = 1:100)
all.equal(attr(r1, "grad"), drop(r2), tolerance = 1e-2)
#> [1] TRUE

#####
# performs gradient descent. 
# 
# Args: 
#   val: starting value. 
#   step_start: starting value for the step length. 
#   n_threads: number of threads to use. 
#   maxit: maximum number of iteration. 
#   eps: convergence threshold to use. 
#   seed: seed to use.
#   c1, c2: parameters for Wolfe condition.
naiv_gradient_descent <- function(val, step_start, n_threads = 4L, 
                                  maxit = 10L, eps = 1e-3, seed = 1L, 
                                  c1 = 1e-3, c2 = .1){
  fun_vals <- step_sizes <- rep(NA_real_, maxit)
  
  gr_new <- par_fn(val, comp_derivs = TRUE, n_threads = n_threads, 
                   seed = seed)
  for(i in 1:maxit){
    gr <- gr_new
    fun_vals[i] <- prev_fun <- c(gr)
    dir <- attr(gr, "grad")
    step <- step_start
    m <- drop(dir %*% attr(gr, "grad"))
    
    max_j <- 11L
    for(j in 1:max_j){
      if(j == max_j)
        warning("Failed to find a decrease")
      new_val <- val + step * dir
      new_val <- get_lchol(cov2cor(get_lchol_inv(new_val)))
      new_fun <- par_fn(new_val, comp_derivs = FALSE, n_threads = n_threads, 
                        seed = seed)
      
      # strong Wolfe conditions
      if(-new_fun <= -prev_fun + c1 * step * m){
        gr_new <- par_fn(new_val, comp_derivs = TRUE, n_threads = n_threads, 
                         seed = seed)
      
        if(abs(drop(attr(gr_new, "grad") %*% dir)) >= c2 * m){
          val <- new_val
          break
        }
      }
      
      step <- step / 2
    }
    
    step_sizes[i] <- step
  }
  
  list(result = get_lchol_inv(val), logml = prev_fun, 
       nit = i, step_sizes = step_sizes, fun_vals = fun_vals)
}

# estimate model parameters
start_val <- numeric(p * (p + 1) / 2)
system.time(res <- naiv_gradient_descent(val = start_val, step_start = .001, 
                                         maxit = 20L, eps = 1e-2))
#>    user  system elapsed 
#>  85.122   0.004  21.862

# compare estimates with truth
norm(res$result - dat$Sigma)
#> [1] 0.615
res$result
#>        [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11] [,12]
#>  [1,] 1.000 0.568 0.416 0.539 0.499 0.588 0.503 0.561 0.575 0.496 0.479 0.542
#>  [2,] 0.568 1.000 0.384 0.504 0.330 0.428 0.521 0.535 0.500 0.484 0.531 0.512
#>  [3,] 0.416 0.384 1.000 0.562 0.510 0.475 0.560 0.514 0.486 0.513 0.567 0.456
#>  [4,] 0.539 0.504 0.562 1.000 0.507 0.594 0.526 0.479 0.548 0.492 0.472 0.539
#>  [5,] 0.499 0.330 0.510 0.507 1.000 0.508 0.437 0.490 0.512 0.506 0.477 0.534
#>  [6,] 0.588 0.428 0.475 0.594 0.508 1.000 0.526 0.509 0.550 0.460 0.443 0.617
#>  [7,] 0.503 0.521 0.560 0.526 0.437 0.526 1.000 0.566 0.438 0.555 0.565 0.564
#>  [8,] 0.561 0.535 0.514 0.479 0.490 0.509 0.566 1.000 0.534 0.622 0.450 0.527
#>  [9,] 0.575 0.500 0.486 0.548 0.512 0.550 0.438 0.534 1.000 0.501 0.537 0.497
#> [10,] 0.496 0.484 0.513 0.492 0.506 0.460 0.555 0.622 0.501 1.000 0.522 0.452
#> [11,] 0.479 0.531 0.567 0.472 0.477 0.443 0.565 0.450 0.537 0.522 1.000 0.467
#> [12,] 0.542 0.512 0.456 0.539 0.534 0.617 0.564 0.527 0.497 0.452 0.467 1.000
#> [13,] 0.578 0.528 0.483 0.464 0.511 0.495 0.557 0.600 0.598 0.539 0.527 0.491
#> [14,] 0.526 0.435 0.568 0.519 0.596 0.576 0.512 0.537 0.542 0.455 0.528 0.457
#> [15,] 0.643 0.484 0.500 0.639 0.507 0.564 0.547 0.566 0.573 0.517 0.497 0.513
#>       [,13] [,14] [,15]
#>  [1,] 0.578 0.526 0.643
#>  [2,] 0.528 0.435 0.484
#>  [3,] 0.483 0.568 0.500
#>  [4,] 0.464 0.519 0.639
#>  [5,] 0.511 0.596 0.507
#>  [6,] 0.495 0.576 0.564
#>  [7,] 0.557 0.512 0.547
#>  [8,] 0.600 0.537 0.566
#>  [9,] 0.598 0.542 0.573
#> [10,] 0.539 0.455 0.517
#> [11,] 0.527 0.528 0.497
#> [12,] 0.491 0.457 0.513
#> [13,] 1.000 0.640 0.640
#> [14,] 0.640 1.000 0.580
#> [15,] 0.640 0.580 1.000
dat$Sigma
#>        [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11] [,12]
#>  [1,] 1.000 0.534 0.394 0.521 0.449 0.536 0.469 0.508 0.528 0.457 0.446 0.504
#>  [2,] 0.534 1.000 0.417 0.500 0.319 0.414 0.486 0.524 0.473 0.504 0.521 0.517
#>  [3,] 0.394 0.417 1.000 0.570 0.509 0.461 0.551 0.516 0.467 0.502 0.528 0.461
#>  [4,] 0.521 0.500 0.570 1.000 0.489 0.540 0.500 0.458 0.550 0.474 0.471 0.504
#>  [5,] 0.449 0.319 0.509 0.489 1.000 0.484 0.383 0.426 0.491 0.456 0.451 0.496
#>  [6,] 0.536 0.414 0.461 0.540 0.484 1.000 0.475 0.459 0.461 0.393 0.409 0.555
#>  [7,] 0.469 0.486 0.551 0.500 0.383 0.475 1.000 0.500 0.438 0.556 0.523 0.516
#>  [8,] 0.508 0.524 0.516 0.458 0.426 0.459 0.500 1.000 0.502 0.612 0.416 0.507
#>  [9,] 0.528 0.473 0.467 0.550 0.491 0.461 0.438 0.502 1.000 0.475 0.506 0.452
#> [10,] 0.457 0.504 0.502 0.474 0.456 0.393 0.556 0.612 0.475 1.000 0.514 0.420
#> [11,] 0.446 0.521 0.528 0.471 0.451 0.409 0.523 0.416 0.506 0.514 1.000 0.420
#> [12,] 0.504 0.517 0.461 0.504 0.496 0.555 0.516 0.507 0.452 0.420 0.420 1.000
#> [13,] 0.517 0.469 0.445 0.465 0.498 0.454 0.545 0.531 0.569 0.469 0.480 0.450
#> [14,] 0.476 0.445 0.565 0.533 0.548 0.544 0.474 0.473 0.513 0.421 0.515 0.433
#> [15,] 0.627 0.506 0.495 0.638 0.481 0.535 0.514 0.511 0.539 0.537 0.493 0.488
#>       [,13] [,14] [,15]
#>  [1,] 0.517 0.476 0.627
#>  [2,] 0.469 0.445 0.506
#>  [3,] 0.445 0.565 0.495
#>  [4,] 0.465 0.533 0.638
#>  [5,] 0.498 0.548 0.481
#>  [6,] 0.454 0.544 0.535
#>  [7,] 0.545 0.474 0.514
#>  [8,] 0.531 0.473 0.511
#>  [9,] 0.569 0.513 0.539
#> [10,] 0.469 0.421 0.537
#> [11,] 0.480 0.515 0.493
#> [12,] 0.450 0.433 0.488
#> [13,] 1.000 0.598 0.615
#> [14,] 0.598 1.000 0.548
#> [15,] 0.615 0.548 1.000

# or plot both of them and compare
do_plot(res$result, dat$Sigma, "Estimates")
```

<img src="man/figures/README-detailed-1.png" width="100%" />

``` r

res$fun_vals # log marginal likelihood estimates at each iteration
#>  [1] -25889 -22498 -22152 -21914 -21796 -21762 -21759 -21722 -21721 -21706
#> [11] -21705 -21703 -21703 -21702 -21702 -21701 -21701 -21701 -21701 -21701

#####
# performs stochastic gradient descent instead (using ADAM).
# 
# Args: 
#   val: starting value. 
#   batch_size: number of observations in each batch. 
#   n_threads: number of threads to use. 
#   maxit: maximum number of iteration. 
#   seed: seed to use.
#   epsilon, alpha, beta_1, beta_2: ADAM parameters.
adam <- function(val, batch_size, n_threads = 4L, maxit = 10L, 
                 seed = 1L, epsilon = 1e-8, alpha = .001, beta_1 = .9, 
                 beta_2 = .999){
  indices <- sample(0:(NROW(dat$seen_obs) - 1L), replace = FALSE)
  blocks <- tapply(indices, (seq_along(indices) - 1L) %/% batch_size, 
                   identity, simplify = FALSE)
  
  n_blocks <- length(blocks)
  n_par <- length(val)
  m <- v <- numeric(n_par)
  fun_vals <- numeric(maxit)
  estimates <- matrix(NA_real_, n_par, maxit)
  i <- -1L
  
  for(k in 1:maxit){
    for(ii in 1:n_blocks){
      i <- i + 1L
      idx_b <- (i %% n_blocks) + 1L
      m_old <- m
      v_old <- v
      res <- par_fn(val, comp_derivs = TRUE, n_threads = n_threads, 
                    seed = seed, indices = blocks[[idx_b]])
      fun_vals[(i %/% n_blocks) + 1L] <- 
        fun_vals[(i %/% n_blocks) + 1L] + c(res)
      
      gr <- attr(res, "grad")
      
      m <- beta_1 * m_old + (1 - beta_1) * gr
      v <- beta_2 * v_old + (1 - beta_2) * gr^2
      
      m_hat <- m / (1 - beta_1^(i + 1))
      v_hat <- v / (1 - beta_2^(i + 1))
      
      val <- val + alpha * m_hat / (sqrt(v_hat) + epsilon)
      val <- get_lchol(cov2cor(get_lchol_inv(val)))
    }
    
    estimates[, k] <- val
  }
  
  list(result = get_lchol_inv(val), fun_vals = fun_vals, 
       estimates = estimates)
}

# estimate the model parameters
set.seed(1)
system.time(res_adam  <- adam(
  val = start_val, alpha = 1e-2, maxit = 10L, batch_size = 100L))
#>    user  system elapsed 
#>  34.358   0.004   9.392

# compare estimates with the truth
norm(res_adam$result - dat$Sigma)
#> [1] 0.599
res_adam$result
#>        [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11] [,12]
#>  [1,] 1.000 0.563 0.405 0.529 0.490 0.582 0.494 0.555 0.567 0.485 0.470 0.530
#>  [2,] 0.563 1.000 0.384 0.504 0.329 0.431 0.522 0.540 0.504 0.486 0.532 0.511
#>  [3,] 0.405 0.384 1.000 0.559 0.502 0.476 0.566 0.506 0.480 0.508 0.565 0.447
#>  [4,] 0.529 0.504 0.559 1.000 0.500 0.595 0.528 0.478 0.547 0.486 0.465 0.532
#>  [5,] 0.490 0.329 0.502 0.500 1.000 0.508 0.432 0.485 0.508 0.503 0.467 0.524
#>  [6,] 0.582 0.431 0.476 0.595 0.508 1.000 0.530 0.509 0.549 0.459 0.437 0.613
#>  [7,] 0.494 0.522 0.566 0.528 0.432 0.530 1.000 0.564 0.440 0.555 0.559 0.560
#>  [8,] 0.555 0.540 0.506 0.478 0.485 0.509 0.564 1.000 0.530 0.619 0.445 0.520
#>  [9,] 0.567 0.504 0.480 0.547 0.508 0.549 0.440 0.530 1.000 0.493 0.534 0.491
#> [10,] 0.485 0.486 0.508 0.486 0.503 0.459 0.555 0.619 0.493 1.000 0.518 0.441
#> [11,] 0.470 0.532 0.565 0.465 0.467 0.437 0.559 0.445 0.534 0.518 1.000 0.459
#> [12,] 0.530 0.511 0.447 0.532 0.524 0.613 0.560 0.520 0.491 0.441 0.459 1.000
#> [13,] 0.568 0.530 0.483 0.456 0.509 0.488 0.554 0.599 0.592 0.541 0.523 0.482
#> [14,] 0.512 0.431 0.569 0.512 0.586 0.576 0.515 0.532 0.537 0.451 0.526 0.444
#> [15,] 0.638 0.482 0.495 0.639 0.497 0.568 0.542 0.560 0.574 0.510 0.488 0.500
#>       [,13] [,14] [,15]
#>  [1,] 0.568 0.512 0.638
#>  [2,] 0.530 0.431 0.482
#>  [3,] 0.483 0.569 0.495
#>  [4,] 0.456 0.512 0.639
#>  [5,] 0.509 0.586 0.497
#>  [6,] 0.488 0.576 0.568
#>  [7,] 0.554 0.515 0.542
#>  [8,] 0.599 0.532 0.560
#>  [9,] 0.592 0.537 0.574
#> [10,] 0.541 0.451 0.510
#> [11,] 0.523 0.526 0.488
#> [12,] 0.482 0.444 0.500
#> [13,] 1.000 0.635 0.632
#> [14,] 0.635 1.000 0.575
#> [15,] 0.632 0.575 1.000
dat$Sigma
#>        [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11] [,12]
#>  [1,] 1.000 0.534 0.394 0.521 0.449 0.536 0.469 0.508 0.528 0.457 0.446 0.504
#>  [2,] 0.534 1.000 0.417 0.500 0.319 0.414 0.486 0.524 0.473 0.504 0.521 0.517
#>  [3,] 0.394 0.417 1.000 0.570 0.509 0.461 0.551 0.516 0.467 0.502 0.528 0.461
#>  [4,] 0.521 0.500 0.570 1.000 0.489 0.540 0.500 0.458 0.550 0.474 0.471 0.504
#>  [5,] 0.449 0.319 0.509 0.489 1.000 0.484 0.383 0.426 0.491 0.456 0.451 0.496
#>  [6,] 0.536 0.414 0.461 0.540 0.484 1.000 0.475 0.459 0.461 0.393 0.409 0.555
#>  [7,] 0.469 0.486 0.551 0.500 0.383 0.475 1.000 0.500 0.438 0.556 0.523 0.516
#>  [8,] 0.508 0.524 0.516 0.458 0.426 0.459 0.500 1.000 0.502 0.612 0.416 0.507
#>  [9,] 0.528 0.473 0.467 0.550 0.491 0.461 0.438 0.502 1.000 0.475 0.506 0.452
#> [10,] 0.457 0.504 0.502 0.474 0.456 0.393 0.556 0.612 0.475 1.000 0.514 0.420
#> [11,] 0.446 0.521 0.528 0.471 0.451 0.409 0.523 0.416 0.506 0.514 1.000 0.420
#> [12,] 0.504 0.517 0.461 0.504 0.496 0.555 0.516 0.507 0.452 0.420 0.420 1.000
#> [13,] 0.517 0.469 0.445 0.465 0.498 0.454 0.545 0.531 0.569 0.469 0.480 0.450
#> [14,] 0.476 0.445 0.565 0.533 0.548 0.544 0.474 0.473 0.513 0.421 0.515 0.433
#> [15,] 0.627 0.506 0.495 0.638 0.481 0.535 0.514 0.511 0.539 0.537 0.493 0.488
#>       [,13] [,14] [,15]
#>  [1,] 0.517 0.476 0.627
#>  [2,] 0.469 0.445 0.506
#>  [3,] 0.445 0.565 0.495
#>  [4,] 0.465 0.533 0.638
#>  [5,] 0.498 0.548 0.481
#>  [6,] 0.454 0.544 0.535
#>  [7,] 0.545 0.474 0.514
#>  [8,] 0.531 0.473 0.511
#>  [9,] 0.569 0.513 0.539
#> [10,] 0.469 0.421 0.537
#> [11,] 0.480 0.515 0.493
#> [12,] 0.450 0.433 0.488
#> [13,] 1.000 0.598 0.615
#> [14,] 0.598 1.000 0.548
#> [15,] 0.615 0.548 1.000

# use plot instead
do_plot(res_adam$result, dat$Sigma, "Estimates (ADAM)")
```

<img src="man/figures/README-detailed-2.png" width="100%" />

``` r

# look at the maximum log marginal likelihood both at the end and after 
# each iteration
log_ml(res_adam$result)
#> [1] -21703
funvals_adam_org <- 
  apply(res_adam$estimates, 2L, function(x) log_ml(get_lchol_inv(x)))
funvals_adam_org
#>  [1] -22679 -22055 -21846 -21775 -21740 -21720 -21710 -21706 -21704 -21703
res_adam$fun_vals # likely lower bounds on the log-marginal likelihood
#>  [1] -23903 -22346 -21953 -21829 -21780 -21755 -21742 -21735 -21732 -21731
res_adam_org <- res_adam

#####
# performs stochastic gradient descent instead (using SVRG).
# 
# Args: 
#   val: starting value. 
#   batch_size: number of observations in each batch. 
#   n_threads: number of threads to use. 
#   maxit: maximum number of iteration. 
#   seed: seed to use.
#   lr: learning rate. 
svrg <- function(val, batch_size, n_threads = 4L, maxit = 10L, 
                 seed = 1L, lr){
  all_indices <- 0:(NROW(dat$seen_obs) - 1L)
  indices <- sample(all_indices, replace = FALSE)
  blocks <- tapply(indices, (seq_along(indices) - 1L) %/% batch_size, 
                   identity, simplify = FALSE)
  
  n_blocks <- length(blocks)
  n_par <- length(val)
  estimates <- matrix(NA_real_, n_par, maxit + 1L)
  fun_vals <- numeric(maxit + 1L)
  estimates[, 1L] <- val
  
  for(k in 1:maxit + 1L){
    old_val <- estimates[, k - 1L]
    old_grs <- sapply(1:n_blocks - 1L, function(ii){
      idx_b <- (ii %% n_blocks) + 1L
      res_old <- par_fn(old_val, comp_derivs = TRUE, n_threads = n_threads, 
                        seed = seed, indices = blocks[[idx_b]])
      c(res_old, attr(res_old, "grad"))
    })
    
    fun_vals[k - 1L] <- sum(old_grs[1, ])
    old_grs <- old_grs[-1L, , drop = FALSE ]
    old_gr <- rowSums(old_grs) / n_blocks
    
    for(ii in 1:n_blocks - 1L){
      idx_b <- (ii %% n_blocks) + 1L
      res <- par_fn(val, comp_derivs = TRUE, n_threads = n_threads, 
                    seed = seed, indices = blocks[[idx_b]])
      fun_vals[k] <- fun_vals[k] + c(res)
      dir <- attr(res, "grad") - old_grs[, ii + 1L] + old_gr
      
      val <- val + lr * dir
      val <- get_lchol(cov2cor(get_lchol_inv(val)))
    }
    
    estimates[, k] <- val
  }
  
  list(result = get_lchol_inv(val), fun_vals = fun_vals[-1L], 
       estimates = estimates[, -1L, drop = FALSE])
}

# estimate the model parameters
set.seed(1)
system.time(res_svrg  <- svrg(
  val = start_val, lr = 1e-3, maxit = 10L, batch_size = 100L))
#>    user  system elapsed 
#>    68.4     0.0    18.7

# compare estimates with the truth
norm(res_svrg$result - dat$Sigma)
#> [1] 0.596
res_svrg$result
#>        [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11] [,12]
#>  [1,] 1.000 0.564 0.410 0.532 0.493 0.581 0.493 0.553 0.568 0.487 0.472 0.536
#>  [2,] 0.564 1.000 0.386 0.504 0.330 0.429 0.522 0.536 0.501 0.485 0.531 0.513
#>  [3,] 0.410 0.386 1.000 0.561 0.509 0.474 0.560 0.513 0.484 0.511 0.567 0.456
#>  [4,] 0.532 0.504 0.561 1.000 0.506 0.594 0.526 0.479 0.547 0.491 0.472 0.539
#>  [5,] 0.493 0.330 0.509 0.506 1.000 0.508 0.436 0.489 0.511 0.505 0.476 0.533
#>  [6,] 0.581 0.429 0.474 0.594 0.508 1.000 0.525 0.507 0.549 0.460 0.442 0.616
#>  [7,] 0.493 0.522 0.560 0.526 0.436 0.525 1.000 0.565 0.441 0.554 0.564 0.563
#>  [8,] 0.553 0.536 0.513 0.479 0.489 0.507 0.565 1.000 0.533 0.619 0.450 0.524
#>  [9,] 0.568 0.501 0.484 0.547 0.511 0.549 0.441 0.533 1.000 0.498 0.535 0.495
#> [10,] 0.487 0.485 0.511 0.491 0.505 0.460 0.554 0.619 0.498 1.000 0.522 0.451
#> [11,] 0.472 0.531 0.567 0.472 0.476 0.442 0.564 0.450 0.535 0.522 1.000 0.465
#> [12,] 0.536 0.513 0.456 0.539 0.533 0.616 0.563 0.524 0.495 0.451 0.465 1.000
#> [13,] 0.571 0.529 0.483 0.463 0.510 0.495 0.556 0.598 0.596 0.538 0.526 0.490
#> [14,] 0.519 0.436 0.568 0.519 0.595 0.575 0.511 0.536 0.541 0.454 0.527 0.456
#> [15,] 0.639 0.484 0.499 0.638 0.505 0.563 0.545 0.564 0.571 0.516 0.495 0.512
#>       [,13] [,14] [,15]
#>  [1,] 0.571 0.519 0.639
#>  [2,] 0.529 0.436 0.484
#>  [3,] 0.483 0.568 0.499
#>  [4,] 0.463 0.519 0.638
#>  [5,] 0.510 0.595 0.505
#>  [6,] 0.495 0.575 0.563
#>  [7,] 0.556 0.511 0.545
#>  [8,] 0.598 0.536 0.564
#>  [9,] 0.596 0.541 0.571
#> [10,] 0.538 0.454 0.516
#> [11,] 0.526 0.527 0.495
#> [12,] 0.490 0.456 0.512
#> [13,] 1.000 0.639 0.638
#> [14,] 0.639 1.000 0.579
#> [15,] 0.638 0.579 1.000
dat$Sigma
#>        [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11] [,12]
#>  [1,] 1.000 0.534 0.394 0.521 0.449 0.536 0.469 0.508 0.528 0.457 0.446 0.504
#>  [2,] 0.534 1.000 0.417 0.500 0.319 0.414 0.486 0.524 0.473 0.504 0.521 0.517
#>  [3,] 0.394 0.417 1.000 0.570 0.509 0.461 0.551 0.516 0.467 0.502 0.528 0.461
#>  [4,] 0.521 0.500 0.570 1.000 0.489 0.540 0.500 0.458 0.550 0.474 0.471 0.504
#>  [5,] 0.449 0.319 0.509 0.489 1.000 0.484 0.383 0.426 0.491 0.456 0.451 0.496
#>  [6,] 0.536 0.414 0.461 0.540 0.484 1.000 0.475 0.459 0.461 0.393 0.409 0.555
#>  [7,] 0.469 0.486 0.551 0.500 0.383 0.475 1.000 0.500 0.438 0.556 0.523 0.516
#>  [8,] 0.508 0.524 0.516 0.458 0.426 0.459 0.500 1.000 0.502 0.612 0.416 0.507
#>  [9,] 0.528 0.473 0.467 0.550 0.491 0.461 0.438 0.502 1.000 0.475 0.506 0.452
#> [10,] 0.457 0.504 0.502 0.474 0.456 0.393 0.556 0.612 0.475 1.000 0.514 0.420
#> [11,] 0.446 0.521 0.528 0.471 0.451 0.409 0.523 0.416 0.506 0.514 1.000 0.420
#> [12,] 0.504 0.517 0.461 0.504 0.496 0.555 0.516 0.507 0.452 0.420 0.420 1.000
#> [13,] 0.517 0.469 0.445 0.465 0.498 0.454 0.545 0.531 0.569 0.469 0.480 0.450
#> [14,] 0.476 0.445 0.565 0.533 0.548 0.544 0.474 0.473 0.513 0.421 0.515 0.433
#> [15,] 0.627 0.506 0.495 0.638 0.481 0.535 0.514 0.511 0.539 0.537 0.493 0.488
#>       [,13] [,14] [,15]
#>  [1,] 0.517 0.476 0.627
#>  [2,] 0.469 0.445 0.506
#>  [3,] 0.445 0.565 0.495
#>  [4,] 0.465 0.533 0.638
#>  [5,] 0.498 0.548 0.481
#>  [6,] 0.454 0.544 0.535
#>  [7,] 0.545 0.474 0.514
#>  [8,] 0.531 0.473 0.511
#>  [9,] 0.569 0.513 0.539
#> [10,] 0.469 0.421 0.537
#> [11,] 0.480 0.515 0.493
#> [12,] 0.450 0.433 0.488
#> [13,] 1.000 0.598 0.615
#> [14,] 0.598 1.000 0.548
#> [15,] 0.615 0.548 1.000

# use plot instead
do_plot(res_svrg$result, dat$Sigma, "Estimates (SVRG)")
```

<img src="man/figures/README-detailed-3.png" width="100%" />

``` r

# look at the maximum log marginal likelihood both at the end and after 
# each iteration
funvals_svrg_org <- res_svrg$fun_vals
funvals_svrg_org[length(funvals_svrg_org)] <- log_ml(res_svrg$result)
funvals_svrg_org
#>  [1] -22413 -21929 -21788 -21738 -21718 -21709 -21705 -21703 -21702 -21702

#####
# we can use better starting values. E.g. something heuristic like: 
#   - transform back into the [0, 1] scale. 
#   - take the middle of the interval and map back. 
#   - compute the partial correlations. 
get_z_hat <- function(lower, upper, code){
  out <- mapply(function(l, u, co){
    if(co <= 1)
      return(u)
    
    a <- if(is.infinite(l)) 0 else pnorm(l)
    b <- if(is.infinite(u)) 1 else pnorm(u)
    qnorm((a + b) / 2)
  }, l = lower, u = upper, c = code)
  dim(out) <- dim(lower)
  out
}
tmp <- get_z_hat(mdgc_obj$lower, mdgc_obj$upper, mdgc_obj$code)

# we also have a C++ function to do this which is faster
all.equal(tmp, mdgc:::get_z_hat(
  mdgc_obj$lower, mdgc_obj$upper, mdgc_obj$code, n_threads = 4L))
#> [1] TRUE

# the latter is faster but both are relatively fast
mark(
  `R version  ` = get_z_hat(mdgc_obj$lower, mdgc_obj$upper, mdgc_obj$code), 
  `C++ verison` = mdgc:::get_z_hat(
  mdgc_obj$lower, mdgc_obj$upper, mdgc_obj$code, n_threads = 4L), 
  min_iterations = 10)
#> Warning: Some expressions had a GC in every iteration; so filtering is disabled.
#> # A tibble: 2 x 6
#>   expression       min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>  <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 R version     76.8ms   82.5ms      12.1    1019KB    32.6 
#> 2 C++ verison  698.6µs  767.6µs    1246.      234KB     6.00

# then we can compute an approximation of the covariance matrix as follows
system.time(chat <- cov2cor(cov(t(tmp), use = "pairwise.complete.obs")))
#>    user  system elapsed 
#>   0.003   0.000   0.003

# the starting value is already quite close
norm(chat - dat$Sigma)
#> [1] 1.48
do_plot(chat, dat$Sigma, "Starting value")
```

<img src="man/figures/README-detailed-4.png" width="100%" />

``` r

# run ADAM again 
start_val <- get_lchol(chat)
set.seed(1)
system.time(res_adam  <- adam(
  val = start_val, alpha = 1e-2, maxit = 5L, batch_size = 100L))
#>    user  system elapsed 
#>  16.618   0.004   4.503

# for comparisons, we also run the code using one thread
set.seed(1)
system.time(res_adam_ser  <- adam(
  val = start_val, alpha = 1e-2, maxit = 5L, batch_size = 100L, 
  n_threads = 1L))
#>    user  system elapsed 
#>    12.6     0.0    12.6

# we get (roughly) the same
norm(res_adam$result - res_adam_ser$result)
#> [1] 0.00453

# plot estimate
norm(res_adam$result - dat$Sigma)
#> [1] 0.611
do_plot(res_adam$result, dat$Sigma, "Estimates (ADAM)")
```

<img src="man/figures/README-detailed-5.png" width="100%" />

``` r

# check log marginal likelihood like before
log_ml(res_adam$result)
#> [1] -21702
funvals_adam <- 
  apply(res_adam$estimates, 2L, function(x) log_ml(get_lchol_inv(x)))
funvals_adam
#> [1] -21715 -21705 -21702 -21702 -21702
res_adam$fun_vals # likely lower bounds on the log-marginal likelihood
#> [1] -21814 -21726 -21734 -21732 -21731

# do the same with SVRG
set.seed(1)
system.time(res_svrg  <- svrg(
  val = start_val, lr = 1e-3, maxit = 5L, batch_size = 100L))
#>    user  system elapsed 
#>  32.797   0.004   8.970

# compare estimates with the truth
norm(res_svrg$result - dat$Sigma)
#> [1] 0.575
res_svrg$result
#>        [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11] [,12]
#>  [1,] 1.000 0.570 0.416 0.537 0.498 0.587 0.501 0.560 0.574 0.494 0.478 0.542
#>  [2,] 0.570 1.000 0.378 0.500 0.324 0.423 0.514 0.529 0.495 0.477 0.526 0.509
#>  [3,] 0.416 0.378 1.000 0.562 0.510 0.474 0.559 0.512 0.484 0.510 0.566 0.455
#>  [4,] 0.537 0.500 0.562 1.000 0.505 0.593 0.524 0.477 0.545 0.488 0.471 0.537
#>  [5,] 0.498 0.324 0.510 0.505 1.000 0.507 0.435 0.487 0.510 0.502 0.475 0.532
#>  [6,] 0.587 0.423 0.474 0.593 0.507 1.000 0.522 0.506 0.545 0.457 0.441 0.613
#>  [7,] 0.501 0.514 0.559 0.524 0.435 0.522 1.000 0.560 0.441 0.548 0.562 0.561
#>  [8,] 0.560 0.529 0.512 0.477 0.487 0.506 0.560 1.000 0.530 0.611 0.448 0.523
#>  [9,] 0.574 0.495 0.484 0.545 0.510 0.545 0.441 0.530 1.000 0.496 0.533 0.495
#> [10,] 0.494 0.477 0.510 0.488 0.502 0.457 0.548 0.611 0.496 1.000 0.519 0.450
#> [11,] 0.478 0.526 0.566 0.471 0.475 0.441 0.562 0.448 0.533 0.519 1.000 0.464
#> [12,] 0.542 0.509 0.455 0.537 0.532 0.613 0.561 0.523 0.495 0.450 0.464 1.000
#> [13,] 0.577 0.524 0.482 0.462 0.509 0.492 0.553 0.596 0.595 0.538 0.525 0.489
#> [14,] 0.525 0.430 0.568 0.517 0.594 0.574 0.509 0.534 0.540 0.454 0.526 0.456
#> [15,] 0.643 0.480 0.500 0.638 0.505 0.562 0.544 0.563 0.570 0.516 0.495 0.512
#>       [,13] [,14] [,15]
#>  [1,] 0.577 0.525 0.643
#>  [2,] 0.524 0.430 0.480
#>  [3,] 0.482 0.568 0.500
#>  [4,] 0.462 0.517 0.638
#>  [5,] 0.509 0.594 0.505
#>  [6,] 0.492 0.574 0.562
#>  [7,] 0.553 0.509 0.544
#>  [8,] 0.596 0.534 0.563
#>  [9,] 0.595 0.540 0.570
#> [10,] 0.538 0.454 0.516
#> [11,] 0.525 0.526 0.495
#> [12,] 0.489 0.456 0.512
#> [13,] 1.000 0.639 0.638
#> [14,] 0.639 1.000 0.579
#> [15,] 0.638 0.579 1.000
dat$Sigma
#>        [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11] [,12]
#>  [1,] 1.000 0.534 0.394 0.521 0.449 0.536 0.469 0.508 0.528 0.457 0.446 0.504
#>  [2,] 0.534 1.000 0.417 0.500 0.319 0.414 0.486 0.524 0.473 0.504 0.521 0.517
#>  [3,] 0.394 0.417 1.000 0.570 0.509 0.461 0.551 0.516 0.467 0.502 0.528 0.461
#>  [4,] 0.521 0.500 0.570 1.000 0.489 0.540 0.500 0.458 0.550 0.474 0.471 0.504
#>  [5,] 0.449 0.319 0.509 0.489 1.000 0.484 0.383 0.426 0.491 0.456 0.451 0.496
#>  [6,] 0.536 0.414 0.461 0.540 0.484 1.000 0.475 0.459 0.461 0.393 0.409 0.555
#>  [7,] 0.469 0.486 0.551 0.500 0.383 0.475 1.000 0.500 0.438 0.556 0.523 0.516
#>  [8,] 0.508 0.524 0.516 0.458 0.426 0.459 0.500 1.000 0.502 0.612 0.416 0.507
#>  [9,] 0.528 0.473 0.467 0.550 0.491 0.461 0.438 0.502 1.000 0.475 0.506 0.452
#> [10,] 0.457 0.504 0.502 0.474 0.456 0.393 0.556 0.612 0.475 1.000 0.514 0.420
#> [11,] 0.446 0.521 0.528 0.471 0.451 0.409 0.523 0.416 0.506 0.514 1.000 0.420
#> [12,] 0.504 0.517 0.461 0.504 0.496 0.555 0.516 0.507 0.452 0.420 0.420 1.000
#> [13,] 0.517 0.469 0.445 0.465 0.498 0.454 0.545 0.531 0.569 0.469 0.480 0.450
#> [14,] 0.476 0.445 0.565 0.533 0.548 0.544 0.474 0.473 0.513 0.421 0.515 0.433
#> [15,] 0.627 0.506 0.495 0.638 0.481 0.535 0.514 0.511 0.539 0.537 0.493 0.488
#>       [,13] [,14] [,15]
#>  [1,] 0.517 0.476 0.627
#>  [2,] 0.469 0.445 0.506
#>  [3,] 0.445 0.565 0.495
#>  [4,] 0.465 0.533 0.638
#>  [5,] 0.498 0.548 0.481
#>  [6,] 0.454 0.544 0.535
#>  [7,] 0.545 0.474 0.514
#>  [8,] 0.531 0.473 0.511
#>  [9,] 0.569 0.513 0.539
#> [10,] 0.469 0.421 0.537
#> [11,] 0.480 0.515 0.493
#> [12,] 0.450 0.433 0.488
#> [13,] 1.000 0.598 0.615
#> [14,] 0.598 1.000 0.548
#> [15,] 0.615 0.548 1.000

# use plot instead
do_plot(res_svrg$result, dat$Sigma, "Estimates (SVRG)")
```

<img src="man/figures/README-detailed-6.png" width="100%" />

``` r

# look at the maximum log marginal likelihood both at the end and after 
# each iteration
funvals_svrg <- res_svrg$fun_vals
funvals_svrg[length(funvals_svrg)] <- log_ml(res_svrg$result)
funvals_svrg
#> [1] -21725 -21707 -21703 -21702 -21701

#####
# compare convergence of the different methods 
#  line: gradient descent. 
#  dashed: ADAM with poor starting values. 
#  dotted: ADAM with better starting values
#  blue: same as ADAM but using SVRG.
lls <- matrix(NA_real_, max(length(res$fun_vals), length(funvals_adam_org), 
                            length(funvals_adam), length(funvals_svrg_org), 
                            length(funvals_svrg)), 5L)
lls[seq_along(res$fun_vals)    , 1] <- res$fun_vals
lls[seq_along(funvals_adam_org), 2] <- funvals_adam_org
lls[seq_along(funvals_adam)    , 3] <- funvals_adam
lls[seq_along(funvals_svrg_org), 4] <- funvals_svrg_org
lls[seq_along(funvals_svrg)    , 5] <- funvals_svrg

par(mfcol = c(1, 1), mar = c(5, 5, 1, 1))
matplot(
  lls, lty = c(1:3, 2:3), col = c(rep("black", 3), rep("darkblue", 2)), 
  bty = "l", type = "l", xlab = "Iteration", 
  ylab = "Log marginal likelihood")
```

<img src="man/figures/README-detailed-7.png" width="100%" />

``` r

# skipping the first n steps
n_skip <- 5L
matplot(
  lls, lty = c(1:3, 2:3), col = c(rep("black", 3), rep("darkblue", 2)), 
  ylim = range(lls[-(1:n_skip), ], na.rm = TRUE), bty = "l", 
  type = "l", xlab = "Iteration", ylab = "Log marginal likelihood")
```

<img src="man/figures/README-detailed-8.png" width="100%" />

## References

<div id="refs" class="references">

<div id="ref-hoff07">

D. Hoff, Peter. 2007. “Extending the Rank Likelihood for Semiparametric
Copula Estimation.” *Ann. Appl. Stat.* 1 (1). The Institute of
Mathematical Statistics: 265–83. <https://doi.org/10.1214/07-AOAS107>.

</div>

<div id="ref-Genz02">

Genz, Alan, and Frank Bretz. 2002. “Comparison of Methods for the
Computation of Multivariate T Probabilities.” *Journal of Computational
and Graphical Statistics* 11 (4). Taylor & Francis: 950–71.
<https://doi.org/10.1198/106186002394>.

</div>

<div id="ref-zhao19">

Zhao, Yuxuan, and Madeleine Udell. 2019. “Missing Value Imputation for
Mixed Data via Gaussian Copula.”

</div>

</div>
