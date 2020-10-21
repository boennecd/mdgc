
# mdgc

[![Build Status on
Travis](https://travis-ci.org/boennecd/mdgc.svg?branch=master,osx)](https://travis-ci.org/boennecd/mdgc)

This package contains a marginal likelihood approach to estimating the
model discussed by Hoff (2007) and Zhao and Udell (2019). That is, a
missing data approach where one uses Gaussian copulas. We have modified
the Fortran code by Genz and Bretz (2002) to supply an approximation of
the gradient for the log marginal likelihood and to use an approximation
of the marginal likelihood similar to the CDF approximation in Genz and
Bretz (2002). We have also used the same Fortran code to perform the
imputation conditional on a correlation matrix. Slides from a
presentation on the packages is provided at
[rpubs.com/boennecd/Gaussian-copula-KTH](https://rpubs.com/boennecd/Gaussian-copula-KTH)
and paper has not been published yet.

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
(this is with a pure R implementation to be fair) while it takes more
than 150 seconds for the MCMC algorithm used by Hoff (2007). These
figures should be kept in mind when looking at the results below.
Importantly, Zhao and Udell (2019) use an approximation in the E-step of
an EM algorithm which is fast but might be crude is some settings. Using
a potentially arbitrarily precise approximation of the log marginal
likelihood is useful if this can be done quickly enough.

We will provide a [quick-example](#quick-example) and [an even shorter
example](#an-even-shorter-example) where we show how to use the methods
in the package to estimate the correlation matrix and to perform the
imputation. We then show a [simulation study](#simulation-study) where
we compare with the method suggested by Zhao and Udell (2019).

We end by providing a [detailed example](#detailed-example) where we:

1.  show how to use the C++ functions and that these provide an
    approximation of the log marginal likelihood and its gradient.
    Moreover, we show that the methods scales well in the number of
    threads.
2.  define functions to perform maximum likelihood estimation.
3.  estimate the parameters using a simple gradient descent algorithm,
    and stochastic gradient descent methods. This serves as an example
    to show how to implement other gradient based methods to estimate
    the model.
4.  show how to improve 4. by using better starting values which are
    quick to compute. As of this writing, this reduces the estimation
    time to about 4 seconds using four threads and about 12 seconds
    using one thread.

The last section is added to give an idea about what is going on under
the hood and can likely be skipped.

## Installation

The packages can be installed from Github by calling:

``` r
remotes::install_github("boennecd/mdgc")
```

### Quick Example

We first simulate a data set and provide an example which shows how to
use the package. The [An Even Shorter Example](#an-even-shorter-example)
section shows a shorter example then what is shown here. You may want to
see this first if you just want to perform some quick imputation.

``` r
# load the packages we need
library(bench)
library(mdgc)
library(missForest, quietly = TRUE)
#> randomForest 4.6-14
#> Type rfNews() to see new features/changes/bug fixes.
# remotes::install_github("udellgroup/mixedgcImp", ref = "5ad6d523d")
library(mixedgcImp)
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
#   Simluated masked data, the true data, and true covariance matrix. 
sim_dat <- function(n, p = 3L, n_lvls = 5L){
  # get the covariance matrix
  Sig <- cov2cor(drop(rWishart(1L, p, diag(p))))
    
  # draw the observations
  truth <- matrix(rnorm(n * p), n) %*% chol(Sig)
  
  # determine the type
  type <- rep(1:3, each = floor((p + 3 - 1) / 3))[1:p]
  is_con <- type == 1L
  is_bin <- type == 2L
  is_ord <- type == 3L
  
  # sample which are masked data 
  is_mask <- matrix(runif(n * p) < .3, n)
  
  # make sure we have no rows with all missing data
  while(any(all_nans <- rowSums(is_mask) == NCOL(is_mask)))
    is_mask[all_nans, ] <- runif(sum(all_nans) * p) < .3
  
  # create observed data
  truth_obs <- data.frame(truth)
  truth_obs[, is_con] <- qexp(pnorm(as.matrix(truth_obs[, is_con])))
  
  bs_border <- 0
  truth_obs[, is_bin] <- 
    truth_obs[, is_bin] > rep(bs_border, each = NROW(truth_obs))
  
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
#>      X1    X2    X3     X4       X5    X6    X7    X8    X9   X10  X11  X12
#> 1 0.237 0.693 0.798 0.0666       NA FALSE FALSE FALSE FALSE  TRUE    E    C
#> 2 0.142    NA    NA 0.0927 0.000152 FALSE    NA  TRUE    NA    NA    E    B
#> 3    NA 0.748 0.629 0.4280       NA    NA  TRUE    NA    NA  TRUE <NA>    A
#> 4 2.702    NA    NA 2.1776 1.700870 FALSE  TRUE  TRUE    NA  TRUE    A <NA>
#> 5 0.925    NA 0.205 0.6046 0.171311  TRUE  TRUE FALSE FALSE FALSE    E    B
#> 6 0.115    NA 1.341     NA       NA FALSE  TRUE  TRUE    NA    NA    E <NA>
#>    X13 X14  X15
#> 1    B   B <NA>
#> 2    A   A    C
#> 3 <NA>   C    E
#> 4    D   B <NA>
#> 5 <NA>   D <NA>
#> 6    A   B <NA>

# assign objects needed for model estimation
mdgc_obj <- get_mdgc(dat$seen_obs)
log_ml_ptr <- get_mdgc_log_ml(mdgc_obj)
start_val <- mdgc_start_value(mdgc_obj)

# this is very fast so we can neglect this when we consider the computation 
# time
mark(`Setup time` = {
  mdgc_obj <- get_mdgc(dat$seen_obs)
  log_ml_ptr <- get_mdgc_log_ml(mdgc_obj)
  start_val <- mdgc_start_value(mdgc_obj)
}, min_iterations = 10)
#> # A tibble: 1 x 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 Setup time   17.2ms   18.2ms      54.0    8.86MB     12.3

# fit the model using two different methods
set.seed(60941821)
system.time(
  fit_adam <- mdgc_fit(
    ptr = log_ml_ptr, vcov = start_val, n_threads = 4L, 
    lr = 1e-3, maxit = 25L, batch_size = 100L, method = "adam", 
     rel_eps = 1e-3, maxpts = 5000L))
#>    user  system elapsed 
#>   36.96    0.00    9.25
set.seed(fit_seed <- 19570958L)
system.time(
  fit_svrg <- mdgc_fit(
    ptr = log_ml_ptr, vcov = start_val, n_threads = 4L, 
    lr = 1e-3, maxit = 25L, batch_size = 100L, method = "svrg", 
    verbose = TRUE, rel_eps = 1e-3, maxpts = 5000L))
#> End if iteration    1 with learning rate 0.00100000
#> Log marginal likelihood approximation is    -23440.91
#> Previous approximate gradient norm was        3399.64
#> 
#> End if iteration    2 with learning rate 0.00098000
#> Log marginal likelihood approximation is    -23389.86
#> Previous approximate gradient norm was        1701.56
#> 
#> End if iteration    3 with learning rate 0.00096040
#> Log marginal likelihood approximation is    -23367.41
#> Previous approximate gradient norm was        1146.69
#> 
#> End if iteration    4 with learning rate 0.00094119
#> Log marginal likelihood approximation is    -23355.64
#> Previous approximate gradient norm was         844.93
#> 
#> End if iteration    5 with learning rate 0.00092237
#> Log marginal likelihood approximation is    -23348.67
#> Previous approximate gradient norm was         667.57
#> 
#> End if iteration    6 with learning rate 0.00090392
#> Log marginal likelihood approximation is    -23344.30
#> Previous approximate gradient norm was         554.26
#> 
#> End if iteration    7 with learning rate 0.00088584
#> Log marginal likelihood approximation is    -23341.28
#> Previous approximate gradient norm was         474.18
#> 
#> End if iteration    8 with learning rate 0.00086813
#> Log marginal likelihood approximation is    -23339.24
#> Previous approximate gradient norm was         419.90
#> 
#> End if iteration    9 with learning rate 0.00085076
#> Log marginal likelihood approximation is    -23337.69
#> Previous approximate gradient norm was         377.09
#> 
#> End if iteration   10 with learning rate 0.00083375
#> Log marginal likelihood approximation is    -23336.49
#> Previous approximate gradient norm was         343.38
#> 
#> End if iteration   11 with learning rate 0.00081707
#> Log marginal likelihood approximation is    -23335.62
#> Previous approximate gradient norm was         320.76
#> 
#> End if iteration   12 with learning rate 0.00080073
#> Log marginal likelihood approximation is    -23334.89
#> Previous approximate gradient norm was         302.69
#> 
#> End if iteration   13 with learning rate 0.00078472
#> Log marginal likelihood approximation is    -23334.32
#> Previous approximate gradient norm was         288.15
#> 
#> End if iteration   14 with learning rate 0.00076902
#> Log marginal likelihood approximation is    -23333.86
#> Previous approximate gradient norm was         276.06
#> 
#> End if iteration   15 with learning rate 0.00075364
#> Log marginal likelihood approximation is    -23333.45
#> Previous approximate gradient norm was         263.95
#> 
#> End if iteration   16 with learning rate 0.00073857
#> Log marginal likelihood approximation is    -23333.13
#> Previous approximate gradient norm was         253.99
#> 
#> End if iteration   17 with learning rate 0.00072380
#> Log marginal likelihood approximation is    -23332.87
#> Previous approximate gradient norm was         272.55
#> 
#> End if iteration   18 with learning rate 0.00070932
#> Log marginal likelihood approximation is    -23332.66
#> Previous approximate gradient norm was         248.81
#> 
#> End if iteration   19 with learning rate 0.00069514
#> Log marginal likelihood approximation is    -23332.46
#> Previous approximate gradient norm was         242.78
#> 
#> End if iteration   20 with learning rate 0.00068123
#> Log marginal likelihood approximation is    -23332.29
#> Previous approximate gradient norm was         238.61
#> 
#> End if iteration   21 with learning rate 0.00066761
#> Log marginal likelihood approximation is    -23332.28
#> Previous approximate gradient norm was         243.06
#> 
#> End if iteration   22 with learning rate 0.00065426
#> Log marginal likelihood approximation is    -23332.10
#> Previous approximate gradient norm was         232.73
#> 
#> End if iteration   23 with learning rate 0.00064117
#> Log marginal likelihood approximation is    -23331.98
#> Previous approximate gradient norm was         228.22
#> 
#> End if iteration   24 with learning rate 0.00062835
#> Log marginal likelihood approximation is    -23331.73
#> Previous approximate gradient norm was         239.21
#> 
#> End if iteration   25 with learning rate 0.00061578
#> Log marginal likelihood approximation is    -23331.65
#> Previous approximate gradient norm was         230.12
#>    user  system elapsed 
#>  77.585   0.004  19.413

# compare the log marginal likelihood 
mdgc_log_ml(vcov = fit_adam$result, ptr = log_ml_ptr, rel_eps = 1e-3)
#> [1] -23351
mdgc_log_ml(vcov = fit_svrg$result, ptr = log_ml_ptr, rel_eps = 1e-3)
#> [1] -23332

# we can use an approximation in the method
set.seed(fit_seed)
system.time(
  fit_svrg_aprx <- mdgc_fit(
    ptr = log_ml_ptr, vcov = start_val, n_threads = 4L, 
    lr = 1e-3, maxit = 25L, batch_size = 100L, method = "svrg", 
    verbose = TRUE, rel_eps = 1e-3, maxpts = 5000L, 
    use_aprx = TRUE))
#> End if iteration    1 with learning rate 0.00100000
#> Log marginal likelihood approximation is    -23440.91
#> Previous approximate gradient norm was        3399.64
#> 
#> End if iteration    2 with learning rate 0.00098000
#> Log marginal likelihood approximation is    -23389.86
#> Previous approximate gradient norm was        1701.56
#> 
#> End if iteration    3 with learning rate 0.00096040
#> Log marginal likelihood approximation is    -23367.41
#> Previous approximate gradient norm was        1146.69
#> 
#> End if iteration    4 with learning rate 0.00094119
#> Log marginal likelihood approximation is    -23355.64
#> Previous approximate gradient norm was         844.93
#> 
#> End if iteration    5 with learning rate 0.00092237
#> Log marginal likelihood approximation is    -23348.67
#> Previous approximate gradient norm was         667.57
#> 
#> End if iteration    6 with learning rate 0.00090392
#> Log marginal likelihood approximation is    -23344.30
#> Previous approximate gradient norm was         554.26
#> 
#> End if iteration    7 with learning rate 0.00088584
#> Log marginal likelihood approximation is    -23341.28
#> Previous approximate gradient norm was         474.18
#> 
#> End if iteration    8 with learning rate 0.00086813
#> Log marginal likelihood approximation is    -23339.24
#> Previous approximate gradient norm was         419.90
#> 
#> End if iteration    9 with learning rate 0.00085076
#> Log marginal likelihood approximation is    -23337.69
#> Previous approximate gradient norm was         377.09
#> 
#> End if iteration   10 with learning rate 0.00083375
#> Log marginal likelihood approximation is    -23336.49
#> Previous approximate gradient norm was         343.38
#> 
#> End if iteration   11 with learning rate 0.00081707
#> Log marginal likelihood approximation is    -23335.62
#> Previous approximate gradient norm was         320.76
#> 
#> End if iteration   12 with learning rate 0.00080073
#> Log marginal likelihood approximation is    -23334.89
#> Previous approximate gradient norm was         302.69
#> 
#> End if iteration   13 with learning rate 0.00078472
#> Log marginal likelihood approximation is    -23334.32
#> Previous approximate gradient norm was         288.15
#> 
#> End if iteration   14 with learning rate 0.00076902
#> Log marginal likelihood approximation is    -23333.86
#> Previous approximate gradient norm was         276.06
#> 
#> End if iteration   15 with learning rate 0.00075364
#> Log marginal likelihood approximation is    -23333.45
#> Previous approximate gradient norm was         263.95
#> 
#> End if iteration   16 with learning rate 0.00073857
#> Log marginal likelihood approximation is    -23333.13
#> Previous approximate gradient norm was         253.99
#> 
#> End if iteration   17 with learning rate 0.00072380
#> Log marginal likelihood approximation is    -23332.87
#> Previous approximate gradient norm was         272.55
#> 
#> End if iteration   18 with learning rate 0.00070932
#> Log marginal likelihood approximation is    -23332.66
#> Previous approximate gradient norm was         248.81
#> 
#> End if iteration   19 with learning rate 0.00069514
#> Log marginal likelihood approximation is    -23332.46
#> Previous approximate gradient norm was         242.78
#> 
#> End if iteration   20 with learning rate 0.00068123
#> Log marginal likelihood approximation is    -23332.29
#> Previous approximate gradient norm was         238.61
#> 
#> End if iteration   21 with learning rate 0.00066761
#> Log marginal likelihood approximation is    -23332.28
#> Previous approximate gradient norm was         243.06
#> 
#> End if iteration   22 with learning rate 0.00065426
#> Log marginal likelihood approximation is    -23332.10
#> Previous approximate gradient norm was         232.73
#> 
#> End if iteration   23 with learning rate 0.00064117
#> Log marginal likelihood approximation is    -23331.98
#> Previous approximate gradient norm was         228.22
#> 
#> End if iteration   24 with learning rate 0.00062835
#> Log marginal likelihood approximation is    -23331.73
#> Previous approximate gradient norm was         239.21
#> 
#> End if iteration   25 with learning rate 0.00061578
#> Log marginal likelihood approximation is    -23331.65
#> Previous approximate gradient norm was         230.12
#>    user  system elapsed 
#>  55.589   0.003  13.900
norm(fit_svrg_aprx$result - fit_svrg$result, "F") # essentially the same
#> [1] 4.45e-08

# compare the estimated correlation matrix with the true value
do_plot <- function(est, truth, main){
  par_old <- par(mfcol = c(1, 3), mar  = c(1, 1, 4, 1))
  on.exit(par(par_old))
  sc <- colorRampPalette(c("Red", "White", "Blue"))(201)
  
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

``` r
# perform the imputation
system.time(
  imp_res <- mdgc_impute(mdgc_obj, fit_svrg$result, rel_eps = 1e-3,
                         maxit = 10000L, n_threads = 4L))
#>    user  system elapsed 
#>   11.07    0.00    3.23

# look at the result for one of the observations
imp_res[2L]
#> [[1]]
#> [[1]]$X1
#> [1] 0.142
#> 
#> [[1]]$X2
#> [1] 2.09
#> 
#> [[1]]$X3
#> [1] 0.249
#> 
#> [[1]]$X4
#> [1] 0.0927
#> 
#> [[1]]$X5
#> [1] 0.000152
#> 
#> [[1]]$X6
#> FALSE  TRUE 
#>     1     0 
#> 
#> [[1]]$X7
#> FALSE  TRUE 
#> 0.196 0.804 
#> 
#> [[1]]$X8
#> FALSE  TRUE 
#>     0     1 
#> 
#> [[1]]$X9
#> FALSE  TRUE 
#> 0.811 0.189 
#> 
#> [[1]]$X10
#> FALSE  TRUE 
#>  0.25  0.75 
#> 
#> [[1]]$X11
#> A B C D E 
#> 0 0 0 0 1 
#> 
#> [[1]]$X12
#> A B C D E 
#> 0 1 0 0 0 
#> 
#> [[1]]$X13
#> A B C D E 
#> 1 0 0 0 0 
#> 
#> [[1]]$X14
#> A B C D E 
#> 1 0 0 0 0 
#> 
#> [[1]]$X15
#> A B C D E 
#> 0 0 1 0 0

# compare with the observed and true data
rbind(truth = dat$truth_obs[2L, ], observed = dat$seen_obs[2L, ])
#>             X1   X2    X3     X4       X5    X6   X7   X8    X9  X10 X11 X12
#> truth    0.142 2.63 0.338 0.0927 0.000152 FALSE TRUE TRUE FALSE TRUE   E   B
#> observed 0.142   NA    NA 0.0927 0.000152 FALSE   NA TRUE    NA   NA   E   B
#>          X13 X14 X15
#> truth      A   A   C
#> observed   A   A   C

# we can threshold the data like this
threshold <- function(org_data, imputed){
  # checks
  stopifnot(NROW(org_data) == length(imputed), 
            is.list(imputed), is.data.frame(org_data))
  
  # threshold
  is_cont <- which(sapply(org_data, is.numeric))
  is_bin  <- which(sapply(org_data, is.logical)) 
  is_ord  <- which(sapply(org_data, is.ordered))
  stopifnot(
    length(is_cont) + length(is_bin) + length(is_ord) == NCOL(org_data))
  is_cat <- c(is_bin, is_ord)
  
  trans_to_df <- function(x){
    if(is.matrix(x))
      as.data.frame(t(x))
    else
      as.data.frame(  x )
  }
  
  out_cont <- trans_to_df(sapply(imputed, function(x) unlist(x[is_cont])))
  out_cat <- trans_to_df(sapply(imputed, function(x) 
    sapply(x[is_cat], which.max)))
  out <- cbind(out_cont, out_cat)
  
  # set factor levels etc. 
  out <- out[, order(c(is_cont, is_bin, is_ord))]
  if(length(is_bin) > 0)
    out[, is_bin] <- out[, is_bin] > 1L
  if(length(is_ord) > 0)
    for(i in is_ord)
      out[[i]] <- ordered(
        unlist(out[[i]]), labels = levels(org_data[, i]))
  
  colnames(out) <- colnames(org_data)
  out
}
thresh_dat <- threshold(dat$seen_obs, imp_res)

# compare thresholded data with observed and true data
head(thresh_dat)
#>      X1    X2    X3     X4       X5    X6    X7    X8    X9   X10 X11 X12 X13
#> 1 0.237 0.693 0.798 0.0666 1.210228 FALSE FALSE FALSE FALSE  TRUE   E   C   B
#> 2 0.142 2.085 0.249 0.0927 0.000152 FALSE  TRUE  TRUE FALSE  TRUE   E   B   A
#> 3 1.301 0.748 0.629 0.4280 0.572526 FALSE  TRUE  TRUE  TRUE  TRUE   E   A   E
#> 4 2.702 0.780 1.136 2.1776 1.700870 FALSE  TRUE  TRUE  TRUE  TRUE   A   B   D
#> 5 0.925 0.913 0.205 0.6046 0.171311  TRUE  TRUE FALSE FALSE FALSE   E   B   D
#> 6 0.115 1.067 1.341 0.2369 0.233466 FALSE  TRUE  TRUE FALSE  TRUE   E   A   A
#>   X14 X15
#> 1   B   D
#> 2   A   C
#> 3   C   E
#> 4   B   E
#> 5   D   E
#> 6   B   B
head(dat$seen_obs)  # observed data
#>      X1    X2    X3     X4       X5    X6    X7    X8    X9   X10  X11  X12
#> 1 0.237 0.693 0.798 0.0666       NA FALSE FALSE FALSE FALSE  TRUE    E    C
#> 2 0.142    NA    NA 0.0927 0.000152 FALSE    NA  TRUE    NA    NA    E    B
#> 3    NA 0.748 0.629 0.4280       NA    NA  TRUE    NA    NA  TRUE <NA>    A
#> 4 2.702    NA    NA 2.1776 1.700870 FALSE  TRUE  TRUE    NA  TRUE    A <NA>
#> 5 0.925    NA 0.205 0.6046 0.171311  TRUE  TRUE FALSE FALSE FALSE    E    B
#> 6 0.115    NA 1.341     NA       NA FALSE  TRUE  TRUE    NA    NA    E <NA>
#>    X13 X14  X15
#> 1    B   B <NA>
#> 2    A   A    C
#> 3 <NA>   C    E
#> 4    D   B <NA>
#> 5 <NA>   D <NA>
#> 6    A   B <NA>
head(dat$truth_obs) # true data
#>      X1    X2    X3     X4       X5    X6    X7    X8    X9   X10 X11 X12 X13
#> 1 0.237 0.693 0.798 0.0666 0.950476 FALSE FALSE FALSE FALSE  TRUE   E   C   B
#> 2 0.142 2.630 0.338 0.0927 0.000152 FALSE  TRUE  TRUE FALSE  TRUE   E   B   A
#> 3 2.864 0.748 0.629 0.4280 1.341650 FALSE  TRUE  TRUE FALSE  TRUE   C   A   D
#> 4 2.702 1.153 0.459 2.1776 1.700870 FALSE  TRUE  TRUE  TRUE  TRUE   A   C   D
#> 5 0.925 0.365 0.205 0.6046 0.171311  TRUE  TRUE FALSE FALSE FALSE   E   B   B
#> 6 0.115 0.563 1.341 0.7184 0.306274 FALSE  TRUE  TRUE FALSE  TRUE   E   A   A
#>   X14 X15
#> 1   B   E
#> 2   A   C
#> 3   C   E
#> 4   B   E
#> 5   D   E
#> 6   B   B

# compare correct categories
get_classif_error <- function(impu_dat, truth = dat$truth_obs, 
                              observed = dat$seen_obs){
  is_cat <- sapply(truth, function(x)
    is.logical(x) || is.ordered(x))
  is_match <- impu_dat[, is_cat] == truth[, is_cat]
  is_match[!is.na(observed[, is_cat])] <- NA_integer_
  1 - colMeans(is_match, na.rm = TRUE)
}
get_classif_error(thresh_dat)
#>    X6    X7    X8    X9   X10   X11   X12   X13   X14   X15 
#> 0.276 0.294 0.225 0.318 0.286 0.574 0.638 0.616 0.592 0.550

# compute RMSE
get_rmse <- function(impu_dat, truth = dat$truth_obs,
                     observed = dat$seen_obs){
  is_con <- sapply(truth, is.numeric)
  err <- as.matrix(impu_dat[, is_con] - truth[, is_con])
  err[!is.na(observed[, is_con])] <- NA_real_
  sqrt(colMeans(err^2, na.rm = TRUE))
}
get_rmse(thresh_dat)
#>    X1    X2    X3    X4    X5 
#> 0.644 0.783 0.651 0.796 0.746

# we can compare this with missForest
miss_forest_arg <- dat$seen_obs
is_log <- sapply(miss_forest_arg, is.logical)
miss_forest_arg[, is_log] <- lapply(miss_forest_arg[, is_log], as.factor)
set.seed(1)
system.time(miss_res <- missForest(miss_forest_arg))
#>   missForest iteration 1 in progress...done!
#>   missForest iteration 2 in progress...done!
#>   missForest iteration 3 in progress...done!
#>   missForest iteration 4 in progress...done!
#>   missForest iteration 5 in progress...done!
#>   missForest iteration 6 in progress...done!
#>   missForest iteration 7 in progress...done!
#>   missForest iteration 8 in progress...done!
#>   missForest iteration 9 in progress...done!
#>    user  system elapsed 
#>  45.841   0.044  45.886

# turn binary variables back to logicals
miss_res$ximp[, is_log] <- lapply(
  miss_res$ximp[, is_log], function(x) as.integer(x) > 1L)

rbind(mdgc       = get_classif_error(thresh_dat),
      missForest = get_classif_error(miss_res$ximp))
#>               X6    X7    X8    X9   X10   X11   X12   X13   X14   X15
#> mdgc       0.276 0.294 0.225 0.318 0.286 0.574 0.638 0.616 0.592 0.550
#> missForest 0.315 0.340 0.304 0.371 0.319 0.651 0.726 0.680 0.673 0.612
rbind(mdgc       = get_rmse(thresh_dat),
      missForest = get_rmse(miss_res$ximp))
#>               X1    X2    X3    X4    X5
#> mdgc       0.644 0.783 0.651 0.796 0.746
#> missForest 0.806 0.848 0.755 0.845 0.842
```

#### An Even Shorter Example

Here is an example where we use the `mdgc` function to do the model
estimation and the imputation:

``` r
# have a data set with missing continous, binary, and ordinal variables
head(dat$seen_obs)
#>      X1    X2    X3     X4       X5    X6    X7    X8    X9   X10  X11  X12
#> 1 0.237 0.693 0.798 0.0666       NA FALSE FALSE FALSE FALSE  TRUE    E    C
#> 2 0.142    NA    NA 0.0927 0.000152 FALSE    NA  TRUE    NA    NA    E    B
#> 3    NA 0.748 0.629 0.4280       NA    NA  TRUE    NA    NA  TRUE <NA>    A
#> 4 2.702    NA    NA 2.1776 1.700870 FALSE  TRUE  TRUE    NA  TRUE    A <NA>
#> 5 0.925    NA 0.205 0.6046 0.171311  TRUE  TRUE FALSE FALSE FALSE    E    B
#> 6 0.115    NA 1.341     NA       NA FALSE  TRUE  TRUE    NA    NA    E <NA>
#>    X13 X14  X15
#> 1    B   B <NA>
#> 2    A   A    C
#> 3 <NA>   C    E
#> 4    D   B <NA>
#> 5 <NA>   D <NA>
#> 6    A   B <NA>
  
# perform the estimation and imputation
set.seed(1)
system.time(res <- mdgc(dat$seen_obs, verbose = TRUE, maxpts = 5000L, 
                        n_threads = 4L, maxit = 25L, use_aprx = TRUE))
#> Estimating the model...
#> End if iteration    1 with learning rate 0.00100000
#> Log marginal likelihood approximation is    -23442.63
#> Previous approximate gradient norm was        3393.96
#> 
#> End if iteration    2 with learning rate 0.00098000
#> Log marginal likelihood approximation is    -23390.91
#> Previous approximate gradient norm was        1694.62
#> 
#> End if iteration    3 with learning rate 0.00096040
#> Log marginal likelihood approximation is    -23368.15
#> Previous approximate gradient norm was        1138.76
#> 
#> End if iteration    4 with learning rate 0.00094119
#> Log marginal likelihood approximation is    -23356.07
#> Previous approximate gradient norm was         841.78
#> 
#> End if iteration    5 with learning rate 0.00092237
#> Log marginal likelihood approximation is    -23348.94
#> Previous approximate gradient norm was         658.80
#> 
#> End if iteration    6 with learning rate 0.00090392
#> Log marginal likelihood approximation is    -23344.40
#> Previous approximate gradient norm was         543.07
#> 
#> End if iteration    7 with learning rate 0.00088584
#> Log marginal likelihood approximation is    -23341.44
#> Previous approximate gradient norm was         466.07
#> 
#> End if iteration    8 with learning rate 0.00086813
#> Log marginal likelihood approximation is    -23339.27
#> Previous approximate gradient norm was         411.44
#> 
#> End if iteration    9 with learning rate 0.00085076
#> Log marginal likelihood approximation is    -23337.75
#> Previous approximate gradient norm was         370.94
#> 
#> End if iteration   10 with learning rate 0.00083375
#> Log marginal likelihood approximation is    -23336.53
#> Previous approximate gradient norm was         338.69
#> 
#> End if iteration   11 with learning rate 0.00081707
#> Log marginal likelihood approximation is    -23335.66
#> Previous approximate gradient norm was         316.54
#> 
#> Performing imputation...
#>    user  system elapsed 
#>   14.74    0.00    4.15

# compare the estimated correlation matrix with the truth
norm(dat$Sigma - res$vcov, "F") / norm(dat$Sigma, "F")
#> [1] 0.0956

# compute the classifcation error and RMSE
get_classif_error(res$ximp)
#>    X6    X7    X8    X9   X10   X11   X12   X13   X14   X15 
#> 0.274 0.290 0.225 0.316 0.286 0.566 0.642 0.616 0.599 0.541
get_rmse(res$ximp)
#>    X1    X2    X3    X4    X5 
#> 0.644 0.783 0.652 0.796 0.747
```

We can compare this with the `mixedgcImp` which uses the method
described in Zhao and Udell (2019):

``` r
# turn the data to a format that can be bassed
dat_pass <- dat$seen_obs
is_cat <- sapply(dat_pass, function(x) is.logical(x) | is.ordered(x))
dat_pass[, is_cat] <- lapply(dat_pass[, is_cat], as.integer)

system.time(imp_apr_em <- impute_mixedgc(dat_pass, eps = 1e-4))
#>    user  system elapsed 
#>    19.7     0.0    19.7

# compare the estimated correlation matrix with the truth
get_rel_err <- function(est, keep = seq_len(NROW(truth)), truth = dat$Sigma)
  norm(truth[keep, keep] - est[keep, keep], "F") / 
  norm(truth, "F")

c(mdgc               = get_rel_err(res$vcov), 
  mixedgcImp         = get_rel_err(imp_apr_em$R), 
  `mdgc bin/ordered` = get_rel_err(res$vcov    , is_cat),
  `mdgc bin/ordered` = get_rel_err(imp_apr_em$R, is_cat),
  `mdgc continuous`  = get_rel_err(res$vcov    , !is_cat),
  `mdgc continuous`  = get_rel_err(imp_apr_em$R, !is_cat))
#>             mdgc       mixedgcImp mdgc bin/ordered mdgc bin/ordered 
#>           0.0956           0.1243           0.0742           0.1083 
#>  mdgc continuous  mdgc continuous 
#>           0.0241           0.0257

# compare the classifcation error and RMSE
imp_apr_res <- as.data.frame(imp_apr_em$Ximp)
is_bin <- sapply(dat$seen_obs, is.logical)
imp_apr_res[, is_bin] <- lapply(imp_apr_res[, is_bin], `>`, e2 = 0)
is_ord <- sapply(dat$seen_obs, is.ordered)
imp_apr_res[, is_ord] <- mapply(function(x, idx)
  ordered(x, labels = levels(dat$seen_obs[[idx]])), 
  x = imp_apr_res[, is_ord], i = which(is_ord), SIMPLIFY = FALSE)

rbind(mdgc       = get_classif_error(res$ximp),
      mixedgcImp = get_classif_error(imp_apr_res))
#>               X6    X7    X8    X9   X10   X11   X12   X13   X14   X15
#> mdgc       0.274 0.290 0.225 0.316 0.286 0.566 0.642 0.616 0.599 0.541
#> mixedgcImp 0.281 0.328 0.232 0.320 0.288 0.626 0.694 0.688 0.609 0.556
rbind(mdgc       = get_rmse(res$ximp),
      mixedgcImp = get_rmse(imp_apr_res))
#>               X1    X2    X3    X4    X5
#> mdgc       0.644 0.783 0.652 0.796 0.747
#> mixedgcImp 0.645 0.789 0.655 0.810 0.755
```

### Simulation Study

We will perform a simulation study in this section to compare different
methods in terms of their computation time and performance. We first
perform the simulation.

``` r
# the seeds we will use
seeds <- c(293498804L, 311878062L, 370718465L, 577520465L, 336825034L, 661670751L, 750947065L, 255824398L, 281823005L, 721186455L, 251974931L, 643568686L, 273097364L, 328663824L, 490259480L, 517126385L, 651705963L, 43381670L, 503505882L, 609251792L, 643001919L, 244401725L, 983414550L, 850590318L, 714971434L, 469416928L, 237089923L, 131313381L, 689679752L, 344399119L, 330840537L, 6287534L, 735760574L, 477353355L, 579527946L, 83409653L, 710142087L, 830103443L, 94094987L, 422058348L, 117889526L, 259750108L, 180244429L, 762680168L, 112163383L, 10802048L, 440434442L, 747282444L, 736836365L, 837303896L, 50697895L, 231661028L, 872653438L, 297024405L, 719108161L, 201103881L, 485890767L, 852715172L, 542126886L, 155221223L, 18987375L, 203133067L, 460377933L, 949381283L, 589083178L, 820719063L, 543339683L, 154667703L, 480316186L, 310795921L, 287317945L, 30587393L, 381290126L, 178269809L, 939854883L, 660119506L, 825302990L, 764135140L, 433746745L, 173637986L, 100446967L, 333304121L, 225525537L, 443031789L, 587486506L, 245392609L, 469144801L, 44073812L, 462948652L, 226692940L, 165285895L, 546908869L, 550076645L, 872290900L, 452044364L, 620131127L, 600097817L, 787537854L, 15915195L, 64220696L)

# gather or compute the results (you may skip this)
res <- lapply(seeds, function(s){
  file_name <- file.path("sim-res", sprintf("seed-%d.RDS", s))
  
  if(file.exists(file_name)){
    message(sprintf("Reading '%s'", file_name))
    out <- readRDS(file_name)
  } else { 
    message(sprintf("Running '%s'", file_name))
    
    # simulate the data
    set.seed(s)
    dat <- sim_dat(2000L, p = 15L)
    
    # fit models and impute
    mdgc_time <- system.time(
      mdgc_res <- mdgc(dat$seen_obs, verbose = FALSE, maxpts = 5000L, 
                        n_threads = 4L, maxit = 25L, use_aprx = TRUE))
    
    dat_pass <- dat$seen_obs
    is_cat <- sapply(dat_pass, function(x) is.logical(x) | is.ordered(x))
    dat_pass[, is_cat] <- lapply(dat_pass[, is_cat], as.integer)
    mixedgc_time <- 
      system.time(mixedgc_res <- impute_mixedgc(dat_pass, eps = 1e-4))
    
    miss_forest_arg <- dat$seen_obs
    is_log <- sapply(miss_forest_arg, is.logical)
    miss_forest_arg[, is_log] <- lapply(
      miss_forest_arg[, is_log], as.factor)
    sink(tempfile())
    miss_time <- system.time(
      miss_res <- missForest(miss_forest_arg, verbose = FALSE))
    sink()
    
    miss_res$ximp[, is_log] <- lapply(
      miss_res$ximp[, is_log], function(x) as.integer(x) > 1L)
    
    # impute using the other estimate
    mdgc_obj <- get_mdgc(dat$seen_obs)
    impu_mixedgc_est <- mdgc_impute(mdgc_obj, mixedgc_res$R)
    impu_mixedgc_est <- threshold(dat$seen_obs, impu_mixedgc_est)
    
    # gather output for the correlation matrix estimates
    vcov_res <- list(truth = dat$Sigma, mdgc = mdgc_res$vcov, 
                     mixedgc = mixedgc_res$R)
    get_rel_err <- function(est, truth, keep = seq_len(NROW(truth)))
      norm(truth[keep, keep] - est[keep, keep], "F") / norm(truth, "F")
    
    vcov_res <- within(vcov_res, {
      mdgc_rel_err    = get_rel_err(mdgc   , truth)
      mixedgc_rel_err = get_rel_err(mixedgc, truth)
    })
    
    # gather output for the imputation 
    mixedgc_imp_res <- as.data.frame(mixedgc_res$Ximp)
    is_bin <- sapply(dat$seen_obs, is.logical)
    mixedgc_imp_res[, is_bin] <- 
      lapply(mixedgc_imp_res[, is_bin, drop = FALSE], `>`, e2 = 0)
    is_ord <- sapply(dat$seen_obs, is.ordered)
    mixedgc_imp_res[, is_ord] <- mapply(function(x, idx)
      ordered(x, labels = levels(dat$seen_obs[[idx]])), 
      x = mixedgc_imp_res[, is_ord, drop = FALSE], 
      i = which(is_ord), SIMPLIFY = FALSE)
    
    get_bin_err <- function(x){
      . <- function(z) z[, is_bin, drop = FALSE]
      get_classif_error(
        .(x), truth = .(dat$truth_obs), observed = .(dat$seen_obs))
    }
    get_ord_err <- function(x){
      . <- function(z) z[, is_ord, drop = FALSE]
      get_classif_error(
        .(x), truth = .(dat$truth_obs), observed = .(dat$seen_obs))
    }
          
    err <- list(
      mdgc_bin = get_bin_err(mdgc_res$ximp), 
      mixedgc_bin = get_bin_err(mixedgc_imp_res), 
      mixed_bin = get_bin_err(impu_mixedgc_est),
      missForest_bin = get_bin_err(miss_res$ximp),
      
      mdgc_class = get_ord_err(mdgc_res$ximp), 
      mixedgc_class = get_ord_err(mixedgc_imp_res), 
      mixed_class = get_ord_err(impu_mixedgc_est),
      missForest_class = get_ord_err(miss_res$ximp),
      
      mdgc_rmse = get_rmse(
        mdgc_res$ximp, truth = dat$truth_obs, observed = dat$seen_obs),
      mixedgc_rmse = get_rmse(
        mixedgc_imp_res, truth = dat$truth_obs, observed = dat$seen_obs),
      mixed_rmse = get_rmse(
        impu_mixedgc_est, truth = dat$truth_obs, observed = dat$seen_obs), 
      missForest_rmse = get_rmse(
        miss_res$ximp, truth = dat$truth_obs, observed = dat$seen_obs))
    
    # gather the times
    times <- list(mdgc = mdgc_time, mixedgc = mixedgc_time, 
                  missForest = miss_time)
    
    # save stats to check convergence
    conv_stats <- list(mdgc = mdgc_res$logLik, 
                       mixedgc = mixedgc_res$loglik)
    
    # save output 
    out <- list(vcov_res = vcov_res, err = err, times = times, 
                conv_stats = conv_stats)
    saveRDS(out, file_name)
  }
  
  # print summary stat to the console while knitting
  . <- function(x)
    message(paste(sprintf("%8.3f", x), collapse = " "))
  with(out, {
    message(paste(
      "mdgc    logLik", 
      paste(sprintf("%.2f", conv_stats$mdgc), collapse = " ")))
    message(paste(
      "mixedgc logLik", 
      paste(sprintf("%.2f", conv_stats$mixedgc), collapse = " ")))
    message(sprintf(
      "Relative correlation matrix estimate errors are %.4f %.4f", 
      vcov_res$mdgc_rel_err, vcov_res$mixedgc_rel_err))
    message(sprintf(
      "Times are %.2f %.2f %.2f", 
      times$mdgc["elapsed"], times$mixedgc["elapsed"], 
      times$missForest["elapsed"]))
    
    message(sprintf(
      "Binary classifcation errors are %.2f %.2f %.2f (%.2f)", 
      mean(err$mdgc_bin), mean(err$mixedgc_bin), 
      mean(err$missForest_bin), mean(err$mixed_bin)))
    # .(err$mdgc_bin)
    # .(err$mixedgc_bin)
    # .(err$mixed_bin)
    
    message(sprintf(
      "Ordinal classifcation errors are %.2f %.2f %.2f (%.2f)", 
      mean(err$mdgc_class), mean(err$mixedgc_class), 
      mean(err$missForest_class), mean(err$mixed_class)))
    # .(err$mdgc_class)
    # .(err$mixedgc_class)
    # .(err$mixed_class)
    
    message(sprintf(
      "Mean RMSEs are %.2f %.2f %.2f (%.2f)",
      mean(err$mdgc_rmse), mean(err$mixedgc_rmse), 
      mean(err$missForest_rmse), mean(err$mixed_rmse)))
    # .(err$mdgc_rmse)
    # .(err$mixedgc_rmse)
    # .(err$mixed_rmse)
    message("")
  })
  
  out  
})
```

The difference in computation time is given below:

``` r
# assign function to show the summary stats
show_sim_stats <- function(v1, v2, v3, what, sub_ele = NULL){
  vals <- sapply(res, function(x) 
    do.call(rbind, x[[what]][c(v1, v2, v3)]), 
    simplify = "array")
  if(!is.null(sub_ele))
    vals <- vals[, sub_ele, , drop = FALSE]
    
  cat("Means and standard errors:\n")
  mea_se <- function(x)
    c(mean = mean(x), SE = sd(x) / sqrt(length(x)))
  print(t(apply(vals, 1L, mea_se)))
  
  cat("\nDifference:\n")
  print(t(apply(
    c(vals[v1, , ]) - 
      aperm(vals[c(v2, v3), , , drop = FALSE], c(3L, 2L, 1L)), 
    3L, mea_se)))
}

# compare estimation time
show_sim_stats(1L, 2L, 3L, "times", "elapsed")
#> Means and standard errors:
#>             mean    SE
#> mdgc        5.39 0.110
#> mixedgc    20.87 0.150
#> missForest 44.85 0.604
#> 
#> Difference:
#>             mean    SE
#> mixedgc    -15.5 0.174
#> missForest -39.5 0.596
```

The summary stats for the relative Frobenius norm between the estimated
and true correlation matrix is given below:

``` r
# relative norms
show_sim_stats("mixedgc_rel_err", "mdgc_rel_err", NULL, "vcov_res")
#> Means and standard errors:
#>                   mean       SE
#> mixedgc_rel_err 0.1187 0.001230
#> mdgc_rel_err    0.0866 0.000969
#> 
#> Difference:
#>                mean       SE
#> mdgc_rel_err 0.0321 0.000909
```

Finally, here are the results for the classification error for the
binary and ordinal outcomes and the root mean square error:

``` r
# the binary variables
show_sim_stats("mdgc_bin", "mixedgc_bin", "missForest_bin", "err")
#> Means and standard errors:
#>                 mean      SE
#> mdgc_bin       0.245 0.00187
#> mixedgc_bin    0.252 0.00186
#> missForest_bin 0.294 0.00188
#> 
#> Difference:
#>                    mean      SE
#> mixedgc_bin    -0.00712 0.00255
#> missForest_bin -0.04974 0.00255

# the ordinal variables
show_sim_stats("mdgc_class", "mixedgc_class", "missForest_class", "err")
#> Means and standard errors:
#>                   mean      SE
#> mdgc_class       0.590 0.00212
#> mixedgc_class    0.623 0.00245
#> missForest_class 0.658 0.00173
#> 
#> Difference:
#>                     mean      SE
#> mixedgc_class    -0.0332 0.00314
#> missForest_class -0.0680 0.00269

# the continous variables
show_sim_stats("mdgc_rmse", "mixedgc_rmse", "missForest_rmse", "err")
#> Means and standard errors:
#>                  mean      SE
#> mdgc_rmse       0.760 0.00402
#> mixedgc_rmse    0.767 0.00404
#> missForest_rmse 0.850 0.00340
#> 
#> Difference:
#>                     mean      SE
#> mixedgc_rmse    -0.00741 0.00576
#> missForest_rmse -0.09019 0.00526
```

It is important to emphasize that missForest is not estimating the true
model.

### Detailed Example

In this section, we show a more detailed example where we use some of
the non-exported functions package. This section is mainly included to
give an idea of what is going on under the hood and to illustrate how
the methods can be used in a user defined optimization function.

``` r
# assign function to evalute the log marginal likelihood
log_ml <- function(...)
  mdgc_log_ml(ptr = log_ml_ptr, ...)

# print the approximate log marginal likelihood at the true parameters
set.seed(1)
print(log_ml(dat$Sigma, n_threads = 1L), digits = 7)
#> [1] -23382.54
print(log_ml(dat$Sigma, n_threads = 1L, use_aprx = TRUE), digits = 7)
#> [1] -23382.69

# check standard error
sd(replicate(20, log_ml(dat$Sigma, n_threads = 4L)))
#> [1] 0.099
sd(replicate(20, log_ml(dat$Sigma, n_threads = 4L, use_aprx = TRUE)))
#> [1] 0.0964

# without reordering
print(log_ml(dat$Sigma, n_threads = 4L, do_reorder = FALSE), digits = 7)
#> [1] -23382.73

# check standard error
sd(replicate(20, log_ml(dat$Sigma, n_threads = 4L, do_reorder = FALSE)))
#> [1] 0.613

# check computation time
mark(
  `1 thread                  ` = 
    log_ml(dat$Sigma                    , n_threads = 1L), 
  `2 threads                 ` = 
    log_ml(dat$Sigma                    , n_threads = 2L),
  `4 threads                 ` = 
    log_ml(dat$Sigma                    , n_threads = 4L), 
  `4 threads (w/ approx)     ` = 
    log_ml(dat$Sigma                    , n_threads = 4L, 
           use_aprx = TRUE),
  `4 threads (w/o reordering)` = 
    log_ml(dat$Sigma, do_reorder = FALSE, n_threads = 4L), 
  min_iterations = 5, check = FALSE)
#> # A tibble: 5 x 6
#>   expression                      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                 <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 1 thread                   841.73ms 853.31ms     1.18     33.9KB        0
#> 2 2 threads                  460.34ms 466.53ms     2.11     33.9KB        0
#> 3 4 threads                  237.98ms 245.13ms     4.10     33.9KB        0
#> 4 4 threads (w/ approx)      145.56ms  154.2ms     6.53     33.9KB        0
#> 5 4 threads (w/o reordering)    3.82s    3.86s     0.257    33.9KB        0

#####
# we can also get an approximation of the gradient
t1 <- log_ml(dat$Sigma, comp_derivs = TRUE, n_threads = 1L, rel_eps = 1e-3)
t2 <- log_ml(dat$Sigma, comp_derivs = TRUE, n_threads = 4L, rel_eps = 1e-3)
all.equal(t1, t2, tolerance = 1e-2)
#> [1] TRUE
  
mark(
  `1 thread                  ` = 
    log_ml(dat$Sigma, comp_derivs = TRUE                    , n_threads = 1L), 
  `2 threads                 ` = 
    log_ml(dat$Sigma, comp_derivs = TRUE                    , n_threads = 2L),
  `4 threads                 ` = 
    log_ml(dat$Sigma, comp_derivs = TRUE                    , n_threads = 4L), 
  `4 threads (w/ approx)     ` = 
    log_ml(dat$Sigma, comp_derivs = TRUE                    , n_threads = 4L, 
           use_aprx = TRUE), 
  `4 threads (w/o reordering)` = 
    log_ml(dat$Sigma, comp_derivs = TRUE, do_reorder = FALSE, n_threads = 4L), 
  min_iterations = 5, check = FALSE)
#> # A tibble: 5 x 6
#>   expression                      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                 <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 1 thread                      8.67s    8.73s    0.114     35.7KB        0
#> 2 2 threads                     4.58s    4.64s    0.215     35.7KB        0
#> 3 4 threads                     2.46s    2.49s    0.398     35.7KB        0
#> 4 4 threads (w/ approx)         1.75s    1.78s    0.562     35.7KB        0
#> 5 4 threads (w/o reordering)    10.5s   11.22s    0.0893    35.7KB        0

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
                 Sigma = S_ex, maxvls = 100000L, abs_eps = -1, 
                 rel_eps = 1e-5, derivs = derivs)
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
#> [1] 9.64e-09
all.equal(c(use_mvtnorm()), c(use_this_pkg()), tolerance = 1e-5)
#> [1] TRUE
mark(mvtnorm = use_mvtnorm(), mdgc = use_this_pkg(), 
     min_iterations = 25, check = FALSE)
#> # A tibble: 2 x 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 mvtnorm       1.1ms   2.94ms      246.    4.43KB        0
#> 2 mdgc         3.07ms   7.72ms      134.    2.49KB        0

sd(replicate(25, use_mvtnorm()))
#> [1] 2.99e-09
sd(replicate(25, use_this_pkg()))
#> [1] 2.86e-09

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
#   rel_eps: relative error for each term.
#   indices: integer vector with which terms to include. 
par_fn <- function(par, seed = NULL, comp_derivs = FALSE, 
                   n_threads = 1L, rel_eps = 1e-2, 
                   indices = 0:(NROW(dat$seen_obs) - 1L)){
  if(!is.null(seed))
    set.seed(seed)
  Arg <- get_lchol_inv(par)
  
  res <- log_ml(Arg, comp_derivs = comp_derivs, indices = indices,
                n_threads = n_threads, rel_eps = rel_eps)
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
r1 <- par_fn(lSig, comp_derivs = TRUE, n_threads = 4L, rel_eps = 1e-3, 
             indices = 1:100)
r2 <- numDeriv::jacobian(par_fn, lSig, seed = 1L, n_threads = 6L, 
                         rel_eps = 1e-3, indices = 1:100)
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
                                  maxit = 25L, eps = 1e-3, seed = 1L, 
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
                                         maxit = 25L, eps = 1e-2))
#>    user  system elapsed 
#>   227.9     0.0    59.2

# compare estimates with truth
norm(res$result - dat$Sigma)
#> [1] 0.572
res$result
#>           [,1]     [,2]    [,3]    [,4]    [,5]    [,6]    [,7]     [,8]
#>  [1,]  1.00000  0.25539 -0.3507  0.1999 -0.0876  0.1413  0.0226 -0.00395
#>  [2,]  0.25539  1.00000 -0.2898 -0.0153 -0.5045 -0.2883  0.1429  0.12470
#>  [3,] -0.35074 -0.28980  1.0000  0.3397  0.3827 -0.0701  0.2222 -0.07865
#>  [4,]  0.19985 -0.01526  0.3397  1.0000  0.1292  0.1015  0.0207 -0.46382
#>  [5,] -0.08759 -0.50447  0.3827  0.1292  1.0000  0.1684 -0.2665 -0.21495
#>  [6,]  0.14128 -0.28830 -0.0701  0.1015  0.1684  1.0000 -0.1013 -0.25842
#>  [7,]  0.02258  0.14293  0.2222  0.0207 -0.2665 -0.1013  1.0000  0.08110
#>  [8,] -0.00395  0.12470 -0.0786 -0.4638 -0.2149 -0.2584  0.0811  1.00000
#>  [9,]  0.24371  0.00529 -0.0973  0.1546  0.1211 -0.0858 -0.2445  0.03323
#> [10,] -0.01942  0.24051  0.0348 -0.0919 -0.0600 -0.4552  0.2419  0.53846
#> [11,] -0.45717 -0.02523  0.1543 -0.3548 -0.0449 -0.4255  0.1492 -0.23124
#> [12,] -0.00237  0.08878 -0.2268 -0.1632  0.0980  0.3280 -0.0468 -0.10100
#> [13,]  0.24702 -0.02689 -0.1881 -0.1090  0.2037 -0.0516  0.1435  0.18150
#> [14,]  0.04648 -0.10524  0.3487  0.2525  0.3927  0.2806 -0.1925 -0.31162
#> [15,]  0.49867 -0.00661 -0.3776  0.0200 -0.0345 -0.3535  0.1560 -0.19918
#>           [,9]   [,10]   [,11]    [,12]   [,13]    [,14]    [,15]
#>  [1,]  0.24371 -0.0194 -0.4572 -0.00237  0.2470  0.04648  0.49867
#>  [2,]  0.00529  0.2405 -0.0252  0.08878 -0.0269 -0.10524 -0.00661
#>  [3,] -0.09725  0.0348  0.1543 -0.22683 -0.1881  0.34868 -0.37762
#>  [4,]  0.15458 -0.0919 -0.3548 -0.16317 -0.1090  0.25253  0.02000
#>  [5,]  0.12108 -0.0600 -0.0449  0.09796  0.2037  0.39271 -0.03451
#>  [6,] -0.08583 -0.4552 -0.4255  0.32802 -0.0516  0.28059 -0.35355
#>  [7,] -0.24452  0.2419  0.1492 -0.04683  0.1435 -0.19252  0.15597
#>  [8,]  0.03323  0.5385 -0.2312 -0.10100  0.1815 -0.31162 -0.19918
#>  [9,]  1.00000  0.0150 -0.0403 -0.33297  0.3207  0.00359  0.19357
#> [10,]  0.01503  1.0000 -0.0533 -0.36161  0.0751 -0.40749  0.10794
#> [11,] -0.04030 -0.0533  1.0000 -0.27225 -0.0972 -0.02600  0.19153
#> [12,] -0.33297 -0.3616 -0.2722  1.00000 -0.0968 -0.09584 -0.21871
#> [13,]  0.32072  0.0751 -0.0972 -0.09683  1.0000  0.34908  0.31284
#> [14,]  0.00359 -0.4075 -0.0260 -0.09584  0.3491  1.00000 -0.30625
#> [15,]  0.19357  0.1079  0.1915 -0.21871  0.3128 -0.30625  1.00000
dat$Sigma
#>           [,1]      [,2]    [,3]    [,4]    [,5]    [,6]    [,7]     [,8]
#>  [1,]  1.00000  0.269549 -0.3697  0.1637 -0.1449  0.1322 -0.0147 -0.03518
#>  [2,]  0.26955  1.000000 -0.3143 -0.0395 -0.5366 -0.3259  0.1208  0.11829
#>  [3,] -0.36970 -0.314321  1.0000  0.3465  0.4054 -0.0265  0.2285 -0.03318
#>  [4,]  0.16370 -0.039486  0.3465  1.0000  0.1247  0.1262 -0.0507 -0.44168
#>  [5,] -0.14492 -0.536557  0.4054  0.1247  1.0000  0.1882 -0.2547 -0.17724
#>  [6,]  0.13219 -0.325944 -0.0265  0.1262  0.1882  1.0000 -0.0924 -0.23529
#>  [7,] -0.01469  0.120819  0.2285 -0.0507 -0.2547 -0.0924  1.0000  0.12393
#>  [8,] -0.03518  0.118288 -0.0332 -0.4417 -0.1772 -0.2353  0.1239  1.00000
#>  [9,]  0.21173  0.013117 -0.1315  0.1793  0.0963 -0.0950 -0.2879 -0.00969
#> [10,] -0.02762  0.201126  0.0502 -0.1142 -0.0134 -0.3572  0.3594  0.53025
#> [11,] -0.47207  0.001810  0.1409 -0.4012 -0.0329 -0.4305  0.1754 -0.20624
#> [12,]  0.00643  0.084866 -0.2180 -0.1354  0.0764  0.3226 -0.0697 -0.08190
#> [13,]  0.21763 -0.000919 -0.1747 -0.1148  0.1956 -0.0125  0.1919  0.16004
#> [14,]  0.03350 -0.067011  0.3727  0.2441  0.3658  0.2234 -0.1878 -0.26035
#> [15,]  0.50698  0.016654 -0.4246 -0.0257 -0.0796 -0.3252  0.1624 -0.25868
#>           [,9]     [,10]     [,11]    [,12]     [,13]    [,14]   [,15]
#>  [1,]  0.21173 -0.027616 -0.472070  0.00643  0.217634  0.03350  0.5070
#>  [2,]  0.01312  0.201126  0.001810  0.08487 -0.000919 -0.06701  0.0167
#>  [3,] -0.13150  0.050198  0.140925 -0.21801 -0.174669  0.37274 -0.4246
#>  [4,]  0.17928 -0.114211 -0.401232 -0.13535 -0.114824  0.24410 -0.0257
#>  [5,]  0.09628 -0.013389 -0.032926  0.07643  0.195601  0.36583 -0.0796
#>  [6,] -0.09496 -0.357185 -0.430495  0.32258 -0.012487  0.22336 -0.3252
#>  [7,] -0.28793  0.359388  0.175371 -0.06967  0.191945 -0.18783  0.1624
#>  [8,] -0.00969  0.530252 -0.206243 -0.08190  0.160038 -0.26035 -0.2587
#>  [9,]  1.00000 -0.086428 -0.071441 -0.31859  0.369579  0.03261  0.2275
#> [10,] -0.08643  1.000000 -0.000806 -0.38475  0.043988 -0.40294  0.1112
#> [11,] -0.07144 -0.000806  1.000000 -0.28584 -0.094194 -0.00869  0.1529
#> [12,] -0.31859 -0.384747 -0.285837  1.00000 -0.071469 -0.06700 -0.2364
#> [13,]  0.36958  0.043988 -0.094194 -0.07147  1.000000  0.37581  0.2932
#> [14,]  0.03261 -0.402944 -0.008686 -0.06700  0.375812  1.00000 -0.3396
#> [15,]  0.22750  0.111167  0.152916 -0.23642  0.293212 -0.33962  1.0000

# or plot both of them and compare
do_plot(res$result, dat$Sigma, "Estimates")
```

<img src="man/figures/README-detailed-1.png" width="100%" />

``` r

res$fun_vals # log marginal likelihood estimates at each iteration
#>  [1] -25892 -24069 -23567 -23479 -23429 -23402 -23390 -23382 -23380 -23380
#> [11] -23357 -23355 -23349 -23347 -23346 -23343 -23342 -23341 -23340 -23339
#> [21] -23337 -23337 -23336 -23335 -23335

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
adam <- function(val, batch_size, n_threads = 4L, maxit = 25L, 
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
  val = start_val, alpha = 1e-2, maxit = 25L, batch_size = 100L))
#>    user  system elapsed 
#> 188.571   0.045  55.628

# compare estimates with the truth
norm(res_adam$result - dat$Sigma)
#> [1] 0.543
res_adam$result
#>           [,1]     [,2]    [,3]    [,4]    [,5]    [,6]    [,7]     [,8]
#>  [1,]  1.00000  0.25373 -0.3548  0.2073 -0.0863  0.1424  0.0161  0.00126
#>  [2,]  0.25373  1.00000 -0.2875 -0.0165 -0.4993 -0.2957  0.1453  0.12089
#>  [3,] -0.35479 -0.28746  1.0000  0.3302  0.3809 -0.0766  0.2313 -0.07174
#>  [4,]  0.20732 -0.01654  0.3302  1.0000  0.1280  0.0994  0.0119 -0.45903
#>  [5,] -0.08629 -0.49933  0.3809  0.1280  1.0000  0.1649 -0.2666 -0.20810
#>  [6,]  0.14236 -0.29572 -0.0766  0.0994  0.1649  1.0000 -0.1072 -0.27645
#>  [7,]  0.01610  0.14533  0.2313  0.0119 -0.2666 -0.1072  1.0000  0.08685
#>  [8,]  0.00126  0.12089 -0.0717 -0.4590 -0.2081 -0.2764  0.0869  1.00000
#>  [9,]  0.24082  0.00680 -0.0954  0.1517  0.1278 -0.0911 -0.2477  0.03644
#> [10,] -0.02329  0.23872  0.0364 -0.0906 -0.0467 -0.4611  0.2507  0.53682
#> [11,] -0.45629 -0.03081  0.1595 -0.3553 -0.0471 -0.4236  0.1635 -0.24617
#> [12,] -0.01068  0.08541 -0.2275 -0.1625  0.0982  0.3243 -0.0493 -0.09172
#> [13,]  0.24734 -0.02741 -0.1897 -0.1036  0.1993 -0.0457  0.1520  0.18700
#> [14,]  0.04407 -0.10215  0.3512  0.2490  0.3974  0.2658 -0.1976 -0.30657
#> [15,]  0.49250 -0.00646 -0.3740  0.0270 -0.0282 -0.3608  0.1592 -0.20501
#>           [,9]   [,10]   [,11]   [,12]   [,13]    [,14]    [,15]
#>  [1,]  0.24082 -0.0233 -0.4563 -0.0107  0.2473  0.04407  0.49250
#>  [2,]  0.00680  0.2387 -0.0308  0.0854 -0.0274 -0.10215 -0.00646
#>  [3,] -0.09544  0.0364  0.1595 -0.2275 -0.1897  0.35121 -0.37400
#>  [4,]  0.15173 -0.0906 -0.3553 -0.1625 -0.1036  0.24897  0.02702
#>  [5,]  0.12785 -0.0467 -0.0471  0.0982  0.1993  0.39741 -0.02819
#>  [6,] -0.09106 -0.4611 -0.4236  0.3243 -0.0457  0.26577 -0.36079
#>  [7,] -0.24774  0.2507  0.1635 -0.0493  0.1520 -0.19763  0.15916
#>  [8,]  0.03644  0.5368 -0.2462 -0.0917  0.1870 -0.30657 -0.20501
#>  [9,]  1.00000  0.0118 -0.0249 -0.3304  0.3228  0.00364  0.19079
#> [10,]  0.01183  1.0000 -0.0473 -0.3658  0.0771 -0.40175  0.12400
#> [11,] -0.02488 -0.0473  1.0000 -0.2761 -0.1028 -0.01529  0.19627
#> [12,] -0.33036 -0.3658 -0.2761  1.0000 -0.0930 -0.11912 -0.23213
#> [13,]  0.32284  0.0771 -0.1028 -0.0930  1.0000  0.35091  0.31526
#> [14,]  0.00364 -0.4018 -0.0153 -0.1191  0.3509  1.00000 -0.30037
#> [15,]  0.19079  0.1240  0.1963 -0.2321  0.3153 -0.30037  1.00000
dat$Sigma
#>           [,1]      [,2]    [,3]    [,4]    [,5]    [,6]    [,7]     [,8]
#>  [1,]  1.00000  0.269549 -0.3697  0.1637 -0.1449  0.1322 -0.0147 -0.03518
#>  [2,]  0.26955  1.000000 -0.3143 -0.0395 -0.5366 -0.3259  0.1208  0.11829
#>  [3,] -0.36970 -0.314321  1.0000  0.3465  0.4054 -0.0265  0.2285 -0.03318
#>  [4,]  0.16370 -0.039486  0.3465  1.0000  0.1247  0.1262 -0.0507 -0.44168
#>  [5,] -0.14492 -0.536557  0.4054  0.1247  1.0000  0.1882 -0.2547 -0.17724
#>  [6,]  0.13219 -0.325944 -0.0265  0.1262  0.1882  1.0000 -0.0924 -0.23529
#>  [7,] -0.01469  0.120819  0.2285 -0.0507 -0.2547 -0.0924  1.0000  0.12393
#>  [8,] -0.03518  0.118288 -0.0332 -0.4417 -0.1772 -0.2353  0.1239  1.00000
#>  [9,]  0.21173  0.013117 -0.1315  0.1793  0.0963 -0.0950 -0.2879 -0.00969
#> [10,] -0.02762  0.201126  0.0502 -0.1142 -0.0134 -0.3572  0.3594  0.53025
#> [11,] -0.47207  0.001810  0.1409 -0.4012 -0.0329 -0.4305  0.1754 -0.20624
#> [12,]  0.00643  0.084866 -0.2180 -0.1354  0.0764  0.3226 -0.0697 -0.08190
#> [13,]  0.21763 -0.000919 -0.1747 -0.1148  0.1956 -0.0125  0.1919  0.16004
#> [14,]  0.03350 -0.067011  0.3727  0.2441  0.3658  0.2234 -0.1878 -0.26035
#> [15,]  0.50698  0.016654 -0.4246 -0.0257 -0.0796 -0.3252  0.1624 -0.25868
#>           [,9]     [,10]     [,11]    [,12]     [,13]    [,14]   [,15]
#>  [1,]  0.21173 -0.027616 -0.472070  0.00643  0.217634  0.03350  0.5070
#>  [2,]  0.01312  0.201126  0.001810  0.08487 -0.000919 -0.06701  0.0167
#>  [3,] -0.13150  0.050198  0.140925 -0.21801 -0.174669  0.37274 -0.4246
#>  [4,]  0.17928 -0.114211 -0.401232 -0.13535 -0.114824  0.24410 -0.0257
#>  [5,]  0.09628 -0.013389 -0.032926  0.07643  0.195601  0.36583 -0.0796
#>  [6,] -0.09496 -0.357185 -0.430495  0.32258 -0.012487  0.22336 -0.3252
#>  [7,] -0.28793  0.359388  0.175371 -0.06967  0.191945 -0.18783  0.1624
#>  [8,] -0.00969  0.530252 -0.206243 -0.08190  0.160038 -0.26035 -0.2587
#>  [9,]  1.00000 -0.086428 -0.071441 -0.31859  0.369579  0.03261  0.2275
#> [10,] -0.08643  1.000000 -0.000806 -0.38475  0.043988 -0.40294  0.1112
#> [11,] -0.07144 -0.000806  1.000000 -0.28584 -0.094194 -0.00869  0.1529
#> [12,] -0.31859 -0.384747 -0.285837  1.00000 -0.071469 -0.06700 -0.2364
#> [13,]  0.36958  0.043988 -0.094194 -0.07147  1.000000  0.37581  0.2932
#> [14,]  0.03261 -0.402944 -0.008686 -0.06700  0.375812  1.00000 -0.3396
#> [15,]  0.22750  0.111167  0.152916 -0.23642  0.293212 -0.33962  1.0000

# use plot instead
do_plot(res_adam$result, dat$Sigma, "Estimates (ADAM)")
```

<img src="man/figures/README-detailed-2.png" width="100%" />

``` r

# look at the maximum log marginal likelihood both at the end and after 
# each iteration
log_ml(res_adam$result)
#> [1] -23334
funvals_adam_org <- 
  apply(res_adam$estimates, 2L, function(x) log_ml(get_lchol_inv(x)))
funvals_adam_org
#>  [1] -24334 -23664 -23437 -23376 -23357 -23349 -23344 -23341 -23339 -23338
#> [11] -23337 -23336 -23336 -23335 -23335 -23335 -23335 -23334 -23334 -23334
#> [21] -23334 -23334 -23334 -23334 -23334
res_adam$fun_vals # likely lower bounds on the log-marginal likelihood
#>  [1] -25083 -23983 -23556 -23430 -23396 -23384 -23378 -23374 -23372 -23370
#> [11] -23369 -23369 -23368 -23368 -23367 -23367 -23367 -23367 -23366 -23366
#> [21] -23366 -23366 -23366 -23366 -23366
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
svrg <- function(val, batch_size, n_threads = 4L, maxit = 25L, 
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
  val = start_val, lr = 1e-3, maxit = 25L, batch_size = 100L))
#>    user  system elapsed 
#> 343.703   0.028 102.594

# compare estimates with the truth
norm(res_svrg$result - dat$Sigma)
#> [1] 0.554
res_svrg$result
#>           [,1]      [,2]    [,3]    [,4]    [,5]    [,6]    [,7]    [,8]
#>  [1,]  1.00000  2.56e-01 -0.3490  0.1981 -0.0891  0.1413  0.0224 -0.0035
#>  [2,]  0.25601  1.00e+00 -0.2859 -0.0193 -0.5061 -0.2899  0.1436  0.1228
#>  [3,] -0.34900 -2.86e-01  1.0000  0.3410  0.3801 -0.0712  0.2229 -0.0780
#>  [4,]  0.19807 -1.93e-02  0.3410  1.0000  0.1297  0.1005  0.0198 -0.4653
#>  [5,] -0.08912 -5.06e-01  0.3801  0.1297  1.0000  0.1688 -0.2675 -0.2139
#>  [6,]  0.14131 -2.90e-01 -0.0712  0.1005  0.1688  1.0000 -0.0950 -0.2705
#>  [7,]  0.02244  1.44e-01  0.2229  0.0198 -0.2675 -0.0950  1.0000  0.0805
#>  [8,] -0.00350  1.23e-01 -0.0780 -0.4653 -0.2139 -0.2705  0.0805  1.0000
#>  [9,]  0.24445  6.94e-03 -0.0970  0.1559  0.1218 -0.0879 -0.2470  0.0334
#> [10,] -0.01843  2.41e-01  0.0357 -0.0904 -0.0565 -0.4561  0.2458  0.5418
#> [11,] -0.45966 -3.08e-02  0.1547 -0.3564 -0.0415 -0.4249  0.1511 -0.2358
#> [12,] -0.00264  8.91e-02 -0.2264 -0.1617  0.0989  0.3293 -0.0453 -0.0961
#> [13,]  0.24455 -2.90e-02 -0.1878 -0.1108  0.2056 -0.0497  0.1478  0.1847
#> [14,]  0.04790 -9.97e-02  0.3484  0.2543  0.3894  0.2806 -0.1969 -0.3140
#> [15,]  0.50213 -7.91e-05 -0.3786  0.0216 -0.0374 -0.3607  0.1584 -0.2050
#>           [,9]    [,10]   [,11]    [,12]   [,13]    [,14]     [,15]
#>  [1,]  0.24445 -0.01843 -0.4597 -0.00264  0.2445  0.04790  5.02e-01
#>  [2,]  0.00694  0.24100 -0.0308  0.08906 -0.0290 -0.09967 -7.91e-05
#>  [3,] -0.09703  0.03571  0.1547 -0.22645 -0.1878  0.34841 -3.79e-01
#>  [4,]  0.15591 -0.09042 -0.3564 -0.16174 -0.1108  0.25426  2.16e-02
#>  [5,]  0.12181 -0.05649 -0.0415  0.09888  0.2056  0.38941 -3.74e-02
#>  [6,] -0.08785 -0.45608 -0.4249  0.32925 -0.0497  0.28065 -3.61e-01
#>  [7,] -0.24698  0.24581  0.1511 -0.04529  0.1478 -0.19687  1.58e-01
#>  [8,]  0.03341  0.54176 -0.2358 -0.09611  0.1847 -0.31402 -2.05e-01
#>  [9,]  1.00000  0.00876 -0.0370 -0.33496  0.3230 -0.00169  1.93e-01
#> [10,]  0.00876  1.00000 -0.0515 -0.36825  0.0728 -0.41018  1.10e-01
#> [11,] -0.03697 -0.05146  1.0000 -0.27173 -0.0967 -0.02500  1.91e-01
#> [12,] -0.33496 -0.36825 -0.2717  1.00000 -0.0978 -0.09861 -2.22e-01
#> [13,]  0.32297  0.07281 -0.0967 -0.09783  1.0000  0.34894  3.12e-01
#> [14,] -0.00169 -0.41018 -0.0250 -0.09861  0.3489  1.00000 -3.07e-01
#> [15,]  0.19293  0.10987  0.1908 -0.22170  0.3117 -0.30668  1.00e+00
dat$Sigma
#>           [,1]      [,2]    [,3]    [,4]    [,5]    [,6]    [,7]     [,8]
#>  [1,]  1.00000  0.269549 -0.3697  0.1637 -0.1449  0.1322 -0.0147 -0.03518
#>  [2,]  0.26955  1.000000 -0.3143 -0.0395 -0.5366 -0.3259  0.1208  0.11829
#>  [3,] -0.36970 -0.314321  1.0000  0.3465  0.4054 -0.0265  0.2285 -0.03318
#>  [4,]  0.16370 -0.039486  0.3465  1.0000  0.1247  0.1262 -0.0507 -0.44168
#>  [5,] -0.14492 -0.536557  0.4054  0.1247  1.0000  0.1882 -0.2547 -0.17724
#>  [6,]  0.13219 -0.325944 -0.0265  0.1262  0.1882  1.0000 -0.0924 -0.23529
#>  [7,] -0.01469  0.120819  0.2285 -0.0507 -0.2547 -0.0924  1.0000  0.12393
#>  [8,] -0.03518  0.118288 -0.0332 -0.4417 -0.1772 -0.2353  0.1239  1.00000
#>  [9,]  0.21173  0.013117 -0.1315  0.1793  0.0963 -0.0950 -0.2879 -0.00969
#> [10,] -0.02762  0.201126  0.0502 -0.1142 -0.0134 -0.3572  0.3594  0.53025
#> [11,] -0.47207  0.001810  0.1409 -0.4012 -0.0329 -0.4305  0.1754 -0.20624
#> [12,]  0.00643  0.084866 -0.2180 -0.1354  0.0764  0.3226 -0.0697 -0.08190
#> [13,]  0.21763 -0.000919 -0.1747 -0.1148  0.1956 -0.0125  0.1919  0.16004
#> [14,]  0.03350 -0.067011  0.3727  0.2441  0.3658  0.2234 -0.1878 -0.26035
#> [15,]  0.50698  0.016654 -0.4246 -0.0257 -0.0796 -0.3252  0.1624 -0.25868
#>           [,9]     [,10]     [,11]    [,12]     [,13]    [,14]   [,15]
#>  [1,]  0.21173 -0.027616 -0.472070  0.00643  0.217634  0.03350  0.5070
#>  [2,]  0.01312  0.201126  0.001810  0.08487 -0.000919 -0.06701  0.0167
#>  [3,] -0.13150  0.050198  0.140925 -0.21801 -0.174669  0.37274 -0.4246
#>  [4,]  0.17928 -0.114211 -0.401232 -0.13535 -0.114824  0.24410 -0.0257
#>  [5,]  0.09628 -0.013389 -0.032926  0.07643  0.195601  0.36583 -0.0796
#>  [6,] -0.09496 -0.357185 -0.430495  0.32258 -0.012487  0.22336 -0.3252
#>  [7,] -0.28793  0.359388  0.175371 -0.06967  0.191945 -0.18783  0.1624
#>  [8,] -0.00969  0.530252 -0.206243 -0.08190  0.160038 -0.26035 -0.2587
#>  [9,]  1.00000 -0.086428 -0.071441 -0.31859  0.369579  0.03261  0.2275
#> [10,] -0.08643  1.000000 -0.000806 -0.38475  0.043988 -0.40294  0.1112
#> [11,] -0.07144 -0.000806  1.000000 -0.28584 -0.094194 -0.00869  0.1529
#> [12,] -0.31859 -0.384747 -0.285837  1.00000 -0.071469 -0.06700 -0.2364
#> [13,]  0.36958  0.043988 -0.094194 -0.07147  1.000000  0.37581  0.2932
#> [14,]  0.03261 -0.402944 -0.008686 -0.06700  0.375812  1.00000 -0.3396
#> [15,]  0.22750  0.111167  0.152916 -0.23642  0.293212 -0.33962  1.0000

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
#>  [1] -24218 -23653 -23472 -23405 -23374 -23359 -23350 -23345 -23341 -23339
#> [11] -23337 -23336 -23335 -23334 -23334 -23333 -23333 -23333 -23332 -23332
#> [21] -23332 -23332 -23332 -23331 -23331

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
#> 1 R version     64.3ms   70.7ms      13.9    1020KB     23.6
#> 2 C++ verison    189µs  278.4µs    3405.      234KB     16.0

# then we can compute an approximation of the covariance matrix as follows
system.time(chat <- cov2cor(cov(t(tmp), use = "pairwise.complete.obs")))
#>    user  system elapsed 
#>   0.002   0.000   0.002

# the starting value is already quite close
norm(chat - dat$Sigma)
#> [1] 0.953
do_plot(chat, dat$Sigma, "Starting value")
```

<img src="man/figures/README-detailed-4.png" width="100%" />

``` r

# run ADAM again 
start_val <- get_lchol(chat)
set.seed(1)
system.time(res_adam  <- adam(
  val = start_val, alpha = 1e-2, maxit = 25L, batch_size = 100L))
#>    user  system elapsed 
#> 209.692   0.012  63.145

# for comparisons, we also run the code using one thread
set.seed(1)
system.time(res_adam_ser  <- adam(
  val = start_val, alpha = 1e-2, maxit = 25L, batch_size = 100L, 
  n_threads = 1L))
#>    user  system elapsed 
#>     175       0     175

# we get (roughly) the same
norm(res_adam$result - res_adam_ser$result)
#> [1] 0.00447

# plot estimate
norm(res_adam$result - dat$Sigma)
#> [1] 0.542
do_plot(res_adam$result, dat$Sigma, "Estimates (ADAM)")
```

<img src="man/figures/README-detailed-5.png" width="100%" />

``` r

# check log marginal likelihood like before
log_ml(res_adam$result)
#> [1] -23333
funvals_adam <- 
  apply(res_adam$estimates, 2L, function(x) log_ml(get_lchol_inv(x)))
funvals_adam
#>  [1] -23395 -23358 -23347 -23342 -23339 -23338 -23337 -23336 -23335 -23335
#> [11] -23334 -23334 -23334 -23334 -23334 -23334 -23334 -23334 -23334 -23333
#> [21] -23333 -23334 -23333 -23333 -23333
res_adam$fun_vals # likely lower bounds on the log-marginal likelihood
#>  [1] -23499 -23397 -23383 -23374 -23372 -23369 -23368 -23367 -23367 -23367
#> [11] -23366 -23366 -23366 -23366 -23366 -23366 -23365 -23365 -23365 -23365
#> [21] -23365 -23365 -23365 -23365 -23365

# do the same with SVRG
set.seed(1)
system.time(res_svrg  <- svrg(
  val = start_val, lr = 1e-3, maxit = 25L, batch_size = 100L))
#>    user  system elapsed 
#> 366.702   0.028 110.666

# compare estimates with the truth
norm(res_svrg$result - dat$Sigma)
#> [1] 0.552
res_svrg$result
#>           [,1]      [,2]    [,3]    [,4]    [,5]    [,6]    [,7]     [,8]
#>  [1,]  1.00000  0.255939 -0.3490  0.1981 -0.0892  0.1413  0.0225 -0.00345
#>  [2,]  0.25594  1.000000 -0.2859 -0.0192 -0.5061 -0.2903  0.1437  0.12279
#>  [3,] -0.34899 -0.285919  1.0000  0.3409  0.3801 -0.0712  0.2229 -0.07803
#>  [4,]  0.19805 -0.019239  0.3409  1.0000  0.1297  0.0999  0.0199 -0.46560
#>  [5,] -0.08920 -0.506093  0.3801  0.1297  1.0000  0.1688 -0.2675 -0.21382
#>  [6,]  0.14131 -0.290331 -0.0712  0.0999  0.1688  1.0000 -0.0943 -0.27215
#>  [7,]  0.02246  0.143665  0.2229  0.0199 -0.2675 -0.0943  1.0000  0.08067
#>  [8,] -0.00345  0.122788 -0.0780 -0.4656 -0.2138 -0.2722  0.0807  1.00000
#>  [9,]  0.24450  0.006981 -0.0970  0.1560  0.1221 -0.0881 -0.2475  0.03316
#> [10,] -0.01826  0.241087  0.0357 -0.0902 -0.0561 -0.4557  0.2457  0.54202
#> [11,] -0.45974 -0.030728  0.1547 -0.3565 -0.0415 -0.4249  0.1512 -0.23615
#> [12,] -0.00258  0.089057 -0.2265 -0.1616  0.0990  0.3291 -0.0453 -0.09555
#> [13,]  0.24435 -0.029009 -0.1879 -0.1109  0.2056 -0.0495  0.1481  0.18516
#> [14,]  0.04787 -0.099457  0.3485  0.2542  0.3894  0.2807 -0.1972 -0.31408
#> [15,]  0.50213 -0.000253 -0.3784  0.0215 -0.0373 -0.3614  0.1582 -0.20560
#>           [,9]    [,10]   [,11]    [,12]   [,13]    [,14]     [,15]
#>  [1,]  0.24450 -0.01826 -0.4597 -0.00258  0.2444  0.04787  0.502129
#>  [2,]  0.00698  0.24109 -0.0307  0.08906 -0.0290 -0.09946 -0.000253
#>  [3,] -0.09700  0.03569  0.1547 -0.22648 -0.1879  0.34854 -0.378431
#>  [4,]  0.15599 -0.09024 -0.3565 -0.16162 -0.1109  0.25423  0.021479
#>  [5,]  0.12206 -0.05614 -0.0415  0.09903  0.2056  0.38945 -0.037331
#>  [6,] -0.08807 -0.45569 -0.4249  0.32909 -0.0495  0.28071 -0.361360
#>  [7,] -0.24749  0.24572  0.1512 -0.04530  0.1481 -0.19718  0.158186
#>  [8,]  0.03316  0.54202 -0.2362 -0.09555  0.1852 -0.31408 -0.205596
#>  [9,]  1.00000  0.00744 -0.0362 -0.33527  0.3234 -0.00237  0.192831
#> [10,]  0.00744  1.00000 -0.0514 -0.36916  0.0726 -0.41036  0.109878
#> [11,] -0.03624 -0.05142  1.0000 -0.27159 -0.0967 -0.02468  0.190951
#> [12,] -0.33527 -0.36916 -0.2716  1.00000 -0.0978 -0.09884 -0.222037
#> [13,]  0.32343  0.07262 -0.0967 -0.09776  1.0000  0.34895  0.311900
#> [14,] -0.00237 -0.41036 -0.0247 -0.09884  0.3490  1.00000 -0.306726
#> [15,]  0.19283  0.10988  0.1910 -0.22204  0.3119 -0.30673  1.000000
dat$Sigma
#>           [,1]      [,2]    [,3]    [,4]    [,5]    [,6]    [,7]     [,8]
#>  [1,]  1.00000  0.269549 -0.3697  0.1637 -0.1449  0.1322 -0.0147 -0.03518
#>  [2,]  0.26955  1.000000 -0.3143 -0.0395 -0.5366 -0.3259  0.1208  0.11829
#>  [3,] -0.36970 -0.314321  1.0000  0.3465  0.4054 -0.0265  0.2285 -0.03318
#>  [4,]  0.16370 -0.039486  0.3465  1.0000  0.1247  0.1262 -0.0507 -0.44168
#>  [5,] -0.14492 -0.536557  0.4054  0.1247  1.0000  0.1882 -0.2547 -0.17724
#>  [6,]  0.13219 -0.325944 -0.0265  0.1262  0.1882  1.0000 -0.0924 -0.23529
#>  [7,] -0.01469  0.120819  0.2285 -0.0507 -0.2547 -0.0924  1.0000  0.12393
#>  [8,] -0.03518  0.118288 -0.0332 -0.4417 -0.1772 -0.2353  0.1239  1.00000
#>  [9,]  0.21173  0.013117 -0.1315  0.1793  0.0963 -0.0950 -0.2879 -0.00969
#> [10,] -0.02762  0.201126  0.0502 -0.1142 -0.0134 -0.3572  0.3594  0.53025
#> [11,] -0.47207  0.001810  0.1409 -0.4012 -0.0329 -0.4305  0.1754 -0.20624
#> [12,]  0.00643  0.084866 -0.2180 -0.1354  0.0764  0.3226 -0.0697 -0.08190
#> [13,]  0.21763 -0.000919 -0.1747 -0.1148  0.1956 -0.0125  0.1919  0.16004
#> [14,]  0.03350 -0.067011  0.3727  0.2441  0.3658  0.2234 -0.1878 -0.26035
#> [15,]  0.50698  0.016654 -0.4246 -0.0257 -0.0796 -0.3252  0.1624 -0.25868
#>           [,9]     [,10]     [,11]    [,12]     [,13]    [,14]   [,15]
#>  [1,]  0.21173 -0.027616 -0.472070  0.00643  0.217634  0.03350  0.5070
#>  [2,]  0.01312  0.201126  0.001810  0.08487 -0.000919 -0.06701  0.0167
#>  [3,] -0.13150  0.050198  0.140925 -0.21801 -0.174669  0.37274 -0.4246
#>  [4,]  0.17928 -0.114211 -0.401232 -0.13535 -0.114824  0.24410 -0.0257
#>  [5,]  0.09628 -0.013389 -0.032926  0.07643  0.195601  0.36583 -0.0796
#>  [6,] -0.09496 -0.357185 -0.430495  0.32258 -0.012487  0.22336 -0.3252
#>  [7,] -0.28793  0.359388  0.175371 -0.06967  0.191945 -0.18783  0.1624
#>  [8,] -0.00969  0.530252 -0.206243 -0.08190  0.160038 -0.26035 -0.2587
#>  [9,]  1.00000 -0.086428 -0.071441 -0.31859  0.369579  0.03261  0.2275
#> [10,] -0.08643  1.000000 -0.000806 -0.38475  0.043988 -0.40294  0.1112
#> [11,] -0.07144 -0.000806  1.000000 -0.28584 -0.094194 -0.00869  0.1529
#> [12,] -0.31859 -0.384747 -0.285837  1.00000 -0.071469 -0.06700 -0.2364
#> [13,]  0.36958  0.043988 -0.094194 -0.07147  1.000000  0.37581  0.2932
#> [14,]  0.03261 -0.402944 -0.008686 -0.06700  0.375812  1.00000 -0.3396
#> [15,]  0.22750  0.111167  0.152916 -0.23642  0.293212 -0.33962  1.0000

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
#>  [1] -23406 -23375 -23359 -23350 -23345 -23341 -23339 -23337 -23336 -23335
#> [11] -23334 -23334 -23333 -23333 -23333 -23332 -23332 -23332 -23332 -23332
#> [21] -23331 -23331 -23331 -23331 -23331

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

<div id="ref-Genz02">

Genz, Alan, and Frank Bretz. 2002. “Comparison of Methods for the
Computation of Multivariate T Probabilities.” *Journal of Computational
and Graphical Statistics* 11 (4). Taylor & Francis: 950–71.
<https://doi.org/10.1198/106186002394>.

</div>

<div id="ref-hoff07">

Hoff, Peter D. 2007. “Extending the Rank Likelihood for Semiparametric
Copula Estimation.” *Ann. Appl. Stat.* 1 (1). The Institute of
Mathematical Statistics: 265–83. <https://doi.org/10.1214/07-AOAS107>.

</div>

<div id="ref-zhao19">

Zhao, Yuxuan, and Madeleine Udell. 2019. “Missing Value Imputation for
Mixed Data via Gaussian Copula.”

</div>

</div>
