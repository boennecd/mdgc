
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
#> 1 Setup time   16.7ms   17.6ms      55.9    8.86MB     12.1

# fit the model using two different methods
set.seed(60941821)
system.time(
  fit_adam <- mdgc_fit(
    ptr = log_ml_ptr, vcov = start_val, n_threads = 4L, 
    lr = 1e-3, maxit = 25L, batch_size = 100L, method = "adam", 
     rel_eps = 1e-3, maxpts = 5000L))
#>    user  system elapsed 
#>   37.52    0.00    9.38
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
#> Log marginal likelihood approximation is    -23332.86
#> Previous approximate gradient norm was         272.13
#> 
#> End if iteration   18 with learning rate 0.00070932
#> Log marginal likelihood approximation is    -23332.65
#> Previous approximate gradient norm was         248.33
#> 
#> End if iteration   19 with learning rate 0.00069514
#> Log marginal likelihood approximation is    -23332.45
#> Previous approximate gradient norm was         243.03
#> 
#> End if iteration   20 with learning rate 0.00068123
#> Log marginal likelihood approximation is    -23332.28
#> Previous approximate gradient norm was         238.29
#> 
#> End if iteration   21 with learning rate 0.00066761
#> Log marginal likelihood approximation is    -23332.31
#> Previous approximate gradient norm was         245.98
#> 
#> End if iteration   22 with learning rate 0.00065426
#> Log marginal likelihood approximation is    -23332.12
#> Previous approximate gradient norm was         234.63
#> 
#> End if iteration   23 with learning rate 0.00064117
#> Log marginal likelihood approximation is    -23331.99
#> Previous approximate gradient norm was         229.97
#> 
#> End if iteration   24 with learning rate 0.00062835
#> Log marginal likelihood approximation is    -23331.72
#> Previous approximate gradient norm was         235.85
#> 
#> End if iteration   25 with learning rate 0.00061578
#> Log marginal likelihood approximation is    -23331.65
#> Previous approximate gradient norm was         227.54
#>    user  system elapsed 
#>  77.267   0.007  19.324

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
#> Previous approximate gradient norm was         302.68
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
#> Log marginal likelihood approximation is    -23332.86
#> Previous approximate gradient norm was         272.13
#> 
#> End if iteration   18 with learning rate 0.00070932
#> Log marginal likelihood approximation is    -23332.65
#> Previous approximate gradient norm was         248.33
#> 
#> End if iteration   19 with learning rate 0.00069514
#> Log marginal likelihood approximation is    -23332.45
#> Previous approximate gradient norm was         243.03
#> 
#> End if iteration   20 with learning rate 0.00068123
#> Log marginal likelihood approximation is    -23332.28
#> Previous approximate gradient norm was         238.29
#> 
#> End if iteration   21 with learning rate 0.00066761
#> Log marginal likelihood approximation is    -23332.31
#> Previous approximate gradient norm was         245.98
#> 
#> End if iteration   22 with learning rate 0.00065426
#> Log marginal likelihood approximation is    -23332.12
#> Previous approximate gradient norm was         234.63
#> 
#> End if iteration   23 with learning rate 0.00064117
#> Log marginal likelihood approximation is    -23331.99
#> Previous approximate gradient norm was         229.97
#> 
#> End if iteration   24 with learning rate 0.00062835
#> Log marginal likelihood approximation is    -23331.72
#> Previous approximate gradient norm was         235.85
#> 
#> End if iteration   25 with learning rate 0.00061578
#> Log marginal likelihood approximation is    -23331.65
#> Previous approximate gradient norm was         227.54
#>    user  system elapsed 
#>  56.725   0.008  14.184
norm(fit_svrg_aprx$result - fit_svrg$result, "F") # essentially the same
#> [1] 9.54e-08

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
#>   17.22    0.00    4.71

# look at the result for one of the observations
imp_res[2L]
#> [[1]]
#> [[1]]$X1
#> [1] 0.142
#> 
#> [[1]]$X2
#> [1] 2.08
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
#> 0.194 0.806 
#> 
#> [[1]]$X8
#> FALSE  TRUE 
#>     0     1 
#> 
#> [[1]]$X9
#> FALSE  TRUE 
#> 0.805 0.195 
#> 
#> [[1]]$X10
#> FALSE  TRUE 
#> 0.248 0.752 
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
#> 1 0.237 0.693 0.798 0.0666 1.210369 FALSE FALSE FALSE FALSE  TRUE   E   C   B
#> 2 0.142 2.085 0.249 0.0927 0.000152 FALSE  TRUE  TRUE FALSE  TRUE   E   B   A
#> 3 1.301 0.748 0.629 0.4280 0.571830 FALSE  TRUE  TRUE  TRUE  TRUE   E   A   E
#> 4 2.702 0.779 1.137 2.1776 1.700870 FALSE  TRUE  TRUE  TRUE  TRUE   A   B   D
#> 5 0.925 0.913 0.205 0.6046 0.171311  TRUE  TRUE FALSE FALSE FALSE   E   B   D
#> 6 0.115 1.069 1.341 0.2374 0.233156 FALSE  TRUE  TRUE FALSE  TRUE   E   A   A
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
#> 0.276 0.292 0.225 0.318 0.288 0.563 0.640 0.623 0.599 0.548

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
#>  47.032   0.024  47.056

# turn binary variables back to logicals
miss_res$ximp[, is_log] <- lapply(
  miss_res$ximp[, is_log], function(x) as.integer(x) > 1L)

rbind(mdgc       = get_classif_error(thresh_dat),
      missForest = get_classif_error(miss_res$ximp))
#>               X6    X7    X8    X9   X10   X11   X12   X13   X14   X15
#> mdgc       0.276 0.292 0.225 0.318 0.288 0.563 0.640 0.623 0.599 0.548
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
#>  15.360   0.001   4.302

# compare the estimated correlation matrix with the truth
norm(dat$Sigma - res$vcov, "F") / norm(dat$Sigma, "F")
#> [1] 0.0956

# compute the classifcation error and RMSE
get_classif_error(res$ximp)
#>    X6    X7    X8    X9   X10   X11   X12   X13   X14   X15 
#> 0.274 0.285 0.226 0.314 0.288 0.574 0.638 0.626 0.605 0.543
get_rmse(res$ximp)
#>    X1    X2    X3    X4    X5 
#> 0.644 0.783 0.651 0.796 0.747
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
#>      20       0      20

# compare the estimated correlation matrix with the truth
get_rel_err <- function(est, keep = seq_len(NROW(truth)), truth = dat$Sigma)
  norm(truth[keep, keep] - est[keep, keep], "F") / 
  norm(truth, "F")

c(mdgc                     = get_rel_err(res$vcov), 
  mixedgcImp               = get_rel_err(imp_apr_em$R), 
  `mdgc bin/ordered`       = get_rel_err(res$vcov    , is_cat),
  `mixedgcImp bin/ordered` = get_rel_err(imp_apr_em$R, is_cat),
  `mdgc continuous`        = get_rel_err(res$vcov    , !is_cat),
  `mixedgcImp continuous`  = get_rel_err(imp_apr_em$R, !is_cat))
#>                   mdgc             mixedgcImp       mdgc bin/ordered 
#>                 0.0956                 0.1243                 0.0742 
#> mixedgcImp bin/ordered        mdgc continuous  mixedgcImp continuous 
#>                 0.1083                 0.0241                 0.0257

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
#> mdgc       0.274 0.285 0.226 0.314 0.288 0.574 0.638 0.626 0.605 0.543
#> mixedgcImp 0.281 0.328 0.232 0.320 0.288 0.626 0.694 0.688 0.609 0.556
rbind(mdgc       = get_rmse(res$ximp),
      mixedgcImp = get_rmse(imp_apr_res))
#>               X1    X2    X3    X4    X5
#> mdgc       0.644 0.783 0.651 0.796 0.747
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
  out <- readRDS(file_name)
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
    
    message(sprintf(
      "Ordinal classifcation errors are %.2f %.2f %.2f (%.2f)", 
      mean(err$mdgc_class), mean(err$mixedgc_class), 
      mean(err$missForest_class), mean(err$mixed_class)))
    
    message(sprintf(
      "Mean RMSEs are %.2f %.2f %.2f (%.2f)",
      mean(err$mdgc_rmse), mean(err$mixedgc_rmse), 
      mean(err$missForest_rmse), mean(err$mixed_rmse)))
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
#> mdgc        5.69 0.117
#> mixedgc    21.45 0.164
#> missForest 44.66 0.594
#> 
#> Difference:
#>             mean    SE
#> mixedgc    -15.8 0.195
#> missForest -39.0 0.590
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
#> mdgc_bin       0.244 0.00187
#> mixedgc_bin    0.252 0.00186
#> missForest_bin 0.294 0.00188
#> 
#> Difference:
#>                   mean      SE
#> mixedgc_bin    -0.0072 0.00255
#> missForest_bin -0.0498 0.00255

# the ordinal variables
show_sim_stats("mdgc_class", "mixedgc_class", "missForest_class", "err")
#> Means and standard errors:
#>                   mean      SE
#> mdgc_class       0.590 0.00215
#> mixedgc_class    0.623 0.00245
#> missForest_class 0.658 0.00173
#> 
#> Difference:
#>                     mean      SE
#> mixedgc_class    -0.0331 0.00317
#> missForest_class -0.0678 0.00273

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
#> 1 1 thread                   847.18ms 849.77ms     1.17     33.9KB        0
#> 2 2 threads                  463.09ms  464.8ms     2.13     33.9KB        0
#> 3 4 threads                  249.03ms 261.54ms     3.77     33.9KB        0
#> 4 4 threads (w/ approx)      156.69ms  166.8ms     6.05     33.9KB        0
#> 5 4 threads (w/o reordering)    3.99s    4.05s     0.247    33.9KB        0

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
#> 1 1 thread                      1.74m    1.76m   0.00950    35.7KB        0
#> 2 2 threads                     54.1s   54.64s   0.0183     35.7KB        0
#> 3 4 threads                    29.87s   30.52s   0.0328     35.7KB        0
#> 4 4 threads (w/ approx)         21.3s   21.42s   0.0467     35.7KB        0
#> 5 4 threads (w/o reordering)    31.5s   31.78s   0.0314     35.7KB        0

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
#> 1 mvtnorm      1.08ms   3.53ms      259.    4.43KB        0
#> 2 mdgc         2.91ms   7.44ms      145.    2.49KB        0

sd(replicate(25, use_mvtnorm()))
#> [1] 3.19e-09
sd(replicate(25, use_this_pkg()))
#> [1] 4.61e-09

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
#>     user   system  elapsed 
#> 3804.843    0.008  961.145

# compare estimates with truth
norm(res$result - dat$Sigma)
#> [1] 0.572
res$result
#>           [,1]     [,2]    [,3]    [,4]    [,5]    [,6]    [,7]     [,8]
#>  [1,]  1.00000  0.25581 -0.3502  0.1990 -0.0883  0.1414  0.0226 -0.00368
#>  [2,]  0.25581  1.00000 -0.2885 -0.0168 -0.5052 -0.2872  0.1425  0.12445
#>  [3,] -0.35019 -0.28845  1.0000  0.3408  0.3823 -0.0704  0.2223 -0.07851
#>  [4,]  0.19901 -0.01677  0.3408  1.0000  0.1295  0.1021  0.0204 -0.46344
#>  [5,] -0.08827 -0.50523  0.3823  0.1295  1.0000  0.1682 -0.2663 -0.21460
#>  [6,]  0.14137 -0.28722 -0.0704  0.1021  0.1682  1.0000 -0.1011 -0.25685
#>  [7,]  0.02260  0.14251  0.2223  0.0204 -0.2663 -0.1011  1.0000  0.08094
#>  [8,] -0.00368  0.12445 -0.0785 -0.4634 -0.2146 -0.2568  0.0809  1.00000
#>  [9,]  0.24362  0.00572 -0.0974  0.1543  0.1203 -0.0860 -0.2440  0.03350
#> [10,] -0.01922  0.24025  0.0351 -0.0918 -0.0604 -0.4546  0.2413  0.53719
#> [11,] -0.45837 -0.02766  0.1549 -0.3546 -0.0432 -0.4247  0.1484 -0.23045
#> [12,] -0.00205  0.08932 -0.2271 -0.1634  0.0975  0.3274 -0.0468 -0.10129
#> [13,]  0.24612 -0.02815 -0.1873 -0.1092  0.2044 -0.0511  0.1431  0.18102
#> [14,]  0.04805 -0.10320  0.3478  0.2530  0.3912  0.2799 -0.1921 -0.31115
#> [15,]  0.50045 -0.00315 -0.3788  0.0204 -0.0367 -0.3533  0.1555 -0.19858
#>           [,9]   [,10]   [,11]    [,12]   [,13]    [,14]    [,15]
#>  [1,]  0.24362 -0.0192 -0.4584 -0.00205  0.2461  0.04805  0.50045
#>  [2,]  0.00572  0.2403 -0.0277  0.08932 -0.0282 -0.10320 -0.00315
#>  [3,] -0.09741  0.0351  0.1549 -0.22709 -0.1873  0.34784 -0.37884
#>  [4,]  0.15430 -0.0918 -0.3546 -0.16339 -0.1092  0.25304  0.02041
#>  [5,]  0.12035 -0.0604 -0.0432  0.09751  0.2044  0.39117 -0.03667
#>  [6,] -0.08601 -0.4546 -0.4247  0.32739 -0.0511  0.27993 -0.35331
#>  [7,] -0.24402  0.2413  0.1484 -0.04678  0.1431 -0.19207  0.15549
#>  [8,]  0.03350  0.5372 -0.2304 -0.10129  0.1810 -0.31115 -0.19858
#>  [9,]  1.00000  0.0157 -0.0405 -0.33263  0.3201  0.00369  0.19380
#> [10,]  0.01569  1.0000 -0.0540 -0.36023  0.0749 -0.40676  0.10827
#> [11,] -0.04047 -0.0540  1.0000 -0.27250 -0.0965 -0.02666  0.19002
#> [12,] -0.33263 -0.3602 -0.2725  1.00000 -0.0966 -0.09595 -0.21818
#> [13,]  0.32011  0.0749 -0.0965 -0.09655  1.0000  0.34902  0.31161
#> [14,]  0.00369 -0.4068 -0.0267 -0.09595  0.3490  1.00000 -0.30590
#> [15,]  0.19380  0.1083  0.1900 -0.21818  0.3116 -0.30590  1.00000
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
#>  [1] -25892 -24070 -23568 -23479 -23429 -23402 -23390 -23382 -23380 -23380
#> [11] -23357 -23355 -23349 -23347 -23346 -23343 -23342 -23341 -23340 -23339
#> [21] -23337 -23337 -23336 -23336 -23336

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
#> 2860.06    0.02  764.43

# compare estimates with the truth
norm(res_adam$result - dat$Sigma)
#> [1] 0.541
res_adam$result
#>           [,1]     [,2]    [,3]    [,4]    [,5]    [,6]    [,7]     [,8]
#>  [1,]  1.00000  0.25381 -0.3548  0.2074 -0.0864  0.1423  0.0159  0.00142
#>  [2,]  0.25381  1.00000 -0.2877 -0.0163 -0.4993 -0.2946  0.1448  0.12035
#>  [3,] -0.35482 -0.28768  1.0000  0.3302  0.3811 -0.0766  0.2310 -0.07155
#>  [4,]  0.20737 -0.01632  0.3302  1.0000  0.1279  0.0989  0.0117 -0.45815
#>  [5,] -0.08641 -0.49927  0.3811  0.1279  1.0000  0.1643 -0.2656 -0.20734
#>  [6,]  0.14234 -0.29460 -0.0766  0.0989  0.1643  1.0000 -0.1072 -0.27474
#>  [7,]  0.01594  0.14478  0.2310  0.0117 -0.2656 -0.1072  1.0000  0.08640
#>  [8,]  0.00142  0.12035 -0.0716 -0.4582 -0.2073 -0.2747  0.0864  1.00000
#>  [9,]  0.23995  0.00670 -0.0949  0.1511  0.1275 -0.0906 -0.2477  0.03580
#> [10,] -0.02307  0.23814  0.0364 -0.0900 -0.0466 -0.4598  0.2502  0.53589
#> [11,] -0.45638 -0.03093  0.1595 -0.3557 -0.0470 -0.4232  0.1633 -0.24644
#> [12,] -0.01047  0.08582 -0.2276 -0.1624  0.0981  0.3239 -0.0491 -0.09178
#> [13,]  0.24712 -0.02763 -0.1897 -0.1036  0.1993 -0.0455  0.1520  0.18653
#> [14,]  0.04416 -0.10189  0.3513  0.2488  0.3973  0.2650 -0.1974 -0.30587
#> [15,]  0.49241 -0.00656 -0.3742  0.0268 -0.0281 -0.3609  0.1590 -0.20515
#>           [,9]   [,10]   [,11]   [,12]   [,13]    [,14]    [,15]
#>  [1,]  0.23995 -0.0231 -0.4564 -0.0105  0.2471  0.04416  0.49241
#>  [2,]  0.00670  0.2381 -0.0309  0.0858 -0.0276 -0.10189 -0.00656
#>  [3,] -0.09489  0.0364  0.1595 -0.2276 -0.1897  0.35128 -0.37415
#>  [4,]  0.15114 -0.0900 -0.3557 -0.1624 -0.1036  0.24876  0.02684
#>  [5,]  0.12749 -0.0466 -0.0470  0.0981  0.1993  0.39728 -0.02815
#>  [6,] -0.09063 -0.4598 -0.4232  0.3239 -0.0455  0.26500 -0.36092
#>  [7,] -0.24773  0.2502  0.1633 -0.0491  0.1520 -0.19740  0.15900
#>  [8,]  0.03580  0.5359 -0.2464 -0.0918  0.1865 -0.30587 -0.20515
#>  [9,]  1.00000  0.0118 -0.0246 -0.3303  0.3218  0.00384  0.18988
#> [10,]  0.01178  1.0000 -0.0477 -0.3653  0.0764 -0.40111  0.12347
#> [11,] -0.02460 -0.0477  1.0000 -0.2759 -0.1028 -0.01560  0.19609
#> [12,] -0.33034 -0.3653 -0.2759  1.0000 -0.0929 -0.11951 -0.23200
#> [13,]  0.32179  0.0764 -0.1028 -0.0929  1.0000  0.35145  0.31538
#> [14,]  0.00384 -0.4011 -0.0156 -0.1195  0.3515  1.00000 -0.30034
#> [15,]  0.18988  0.1235  0.1961 -0.2320  0.3154 -0.30034  1.00000
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
#>  [1] -24335 -23665 -23437 -23376 -23357 -23349 -23344 -23341 -23339 -23338
#> [11] -23337 -23336 -23336 -23336 -23335 -23335 -23335 -23335 -23334 -23334
#> [21] -23334 -23334 -23334 -23334 -23334
res_adam$fun_vals # likely lower bounds on the log-marginal likelihood
#>  [1] -25084 -23983 -23556 -23430 -23396 -23384 -23378 -23374 -23372 -23371
#> [11] -23369 -23369 -23368 -23368 -23367 -23367 -23367 -23366 -23366 -23366
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
#>     user   system  elapsed 
#> 5588.044    0.043 1498.003

# compare estimates with the truth
norm(res_svrg$result - dat$Sigma)
#> [1] 0.553
res_svrg$result
#>           [,1]      [,2]    [,3]    [,4]    [,5]    [,6]    [,7]     [,8]
#>  [1,]  1.00000  0.256030 -0.3490  0.1981 -0.0892  0.1413  0.0224 -0.00338
#>  [2,]  0.25603  1.000000 -0.2861 -0.0192 -0.5061 -0.2895  0.1432  0.12258
#>  [3,] -0.34899 -0.286111  1.0000  0.3410  0.3803 -0.0715  0.2228 -0.07789
#>  [4,]  0.19810 -0.019245  0.3410  1.0000  0.1296  0.1000  0.0197 -0.46495
#>  [5,] -0.08923 -0.506140  0.3803  0.1296  1.0000  0.1685 -0.2670 -0.21372
#>  [6,]  0.14135 -0.289516 -0.0715  0.1000  0.1685  1.0000 -0.0951 -0.26977
#>  [7,]  0.02241  0.143211  0.2228  0.0197 -0.2670 -0.0951  1.0000  0.08004
#>  [8,] -0.00338  0.122584 -0.0779 -0.4649 -0.2137 -0.2698  0.0800  1.00000
#>  [9,]  0.24429  0.006695 -0.0968  0.1557  0.1218 -0.0880 -0.2468  0.03319
#> [10,] -0.01825  0.240650  0.0357 -0.0900 -0.0563 -0.4555  0.2456  0.54095
#> [11,] -0.45969 -0.030714  0.1546 -0.3568 -0.0416 -0.4244  0.1508 -0.23549
#> [12,] -0.00251  0.089377 -0.2265 -0.1619  0.0987  0.3287 -0.0452 -0.09580
#> [13,]  0.24437 -0.029059 -0.1878 -0.1106  0.2056 -0.0496  0.1471  0.18469
#> [14,]  0.04807 -0.099502  0.3485  0.2542  0.3893  0.2804 -0.1966 -0.31357
#> [15,]  0.50226 -0.000171 -0.3786  0.0217 -0.0374 -0.3605  0.1581 -0.20498
#>           [,9]   [,10]   [,11]    [,12]   [,13]    [,14]     [,15]
#>  [1,]  0.24429 -0.0183 -0.4597 -0.00251  0.2444  0.04807  0.502258
#>  [2,]  0.00670  0.2407 -0.0307  0.08938 -0.0291 -0.09950 -0.000171
#>  [3,] -0.09678  0.0357  0.1546 -0.22655 -0.1878  0.34852 -0.378644
#>  [4,]  0.15574 -0.0900 -0.3568 -0.16190 -0.1106  0.25418  0.021682
#>  [5,]  0.12175 -0.0563 -0.0416  0.09868  0.2056  0.38928 -0.037428
#>  [6,] -0.08800 -0.4555 -0.4244  0.32869 -0.0496  0.28036 -0.360468
#>  [7,] -0.24682  0.2456  0.1508 -0.04522  0.1471 -0.19657  0.158112
#>  [8,]  0.03319  0.5410 -0.2355 -0.09580  0.1847 -0.31357 -0.204979
#>  [9,]  1.00000  0.0087 -0.0368 -0.33456  0.3227 -0.00182  0.192754
#> [10,]  0.00870  1.0000 -0.0520 -0.36786  0.0725 -0.40921  0.109304
#> [11,] -0.03684 -0.0520  1.0000 -0.27169 -0.0968 -0.02514  0.190500
#> [12,] -0.33456 -0.3679 -0.2717  1.00000 -0.0975 -0.09901 -0.221769
#> [13,]  0.32269  0.0725 -0.0968 -0.09752  1.0000  0.34909  0.311547
#> [14,] -0.00182 -0.4092 -0.0251 -0.09901  0.3491  1.00000 -0.306617
#> [15,]  0.19275  0.1093  0.1905 -0.22177  0.3115 -0.30662  1.000000
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
#>  [1] -24218 -23654 -23472 -23405 -23375 -23359 -23350 -23345 -23341 -23339
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
#> 1 R version     70.8ms   77.6ms      12.7    1020KB     21.6
#> 2 C++ verison  194.1µs  282.9µs    3354.      234KB     16.0

# then we can compute an approximation of the covariance matrix as follows
system.time(chat <- cov2cor(cov(t(tmp), use = "pairwise.complete.obs")))
#>    user  system elapsed 
#>   0.003   0.000   0.003

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
#>     user   system  elapsed 
#> 3019.570    0.044  805.350

# for comparisons, we also run the code using one thread
set.seed(1)
system.time(res_adam_ser  <- adam(
  val = start_val, alpha = 1e-2, maxit = 25L, batch_size = 100L, 
  n_threads = 1L))
#>     user   system  elapsed 
#> 2609.191    0.008 2609.324

# we get (roughly) the same
norm(res_adam$result - res_adam_ser$result)
#> [1] 0.000531

# plot estimate
norm(res_adam$result - dat$Sigma)
#> [1] 0.539
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
#>  [1] -23395 -23358 -23347 -23342 -23339 -23338 -23336 -23336 -23335 -23335
#> [11] -23335 -23334 -23334 -23334 -23334 -23334 -23334 -23333 -23334 -23334
#> [21] -23334 -23333 -23333 -23333 -23333
res_adam$fun_vals # likely lower bounds on the log-marginal likelihood
#>  [1] -23499 -23397 -23383 -23374 -23372 -23370 -23368 -23368 -23367 -23367
#> [11] -23366 -23366 -23366 -23366 -23366 -23365 -23365 -23365 -23365 -23365
#> [21] -23365 -23365 -23365 -23365 -23365

# do the same with SVRG
set.seed(1)
system.time(res_svrg  <- svrg(
  val = start_val, lr = 1e-3, maxit = 25L, batch_size = 100L))
#>     user   system  elapsed 
#> 5905.350    0.052 1578.056

# compare estimates with the truth
norm(res_svrg$result - dat$Sigma)
#> [1] 0.551
res_svrg$result
#>           [,1]     [,2]    [,3]    [,4]    [,5]    [,6]    [,7]    [,8]
#>  [1,]  1.00000  0.25601 -0.3489  0.1981 -0.0893  0.1413  0.0224 -0.0034
#>  [2,]  0.25601  1.00000 -0.2861 -0.0193 -0.5061 -0.2898  0.1433  0.1224
#>  [3,] -0.34894 -0.28606  1.0000  0.3410  0.3802 -0.0715  0.2228 -0.0779
#>  [4,]  0.19807 -0.01927  0.3410  1.0000  0.1296  0.0996  0.0198 -0.4651
#>  [5,] -0.08929 -0.50611  0.3802  0.1296  1.0000  0.1685 -0.2671 -0.2136
#>  [6,]  0.14131 -0.28984 -0.0715  0.0996  0.1685  1.0000 -0.0945 -0.2712
#>  [7,]  0.02242  0.14335  0.2228  0.0198 -0.2671 -0.0945  1.0000  0.0802
#>  [8,] -0.00340  0.12235 -0.0779 -0.4651 -0.2136 -0.2712  0.0802  1.0000
#>  [9,]  0.24433  0.00674 -0.0967  0.1559  0.1220 -0.0881 -0.2471  0.0330
#> [10,] -0.01812  0.24073  0.0357 -0.0898 -0.0559 -0.4554  0.2457  0.5410
#> [11,] -0.45972 -0.03070  0.1546 -0.3569 -0.0416 -0.4244  0.1511 -0.2360
#> [12,] -0.00257  0.08931 -0.2265 -0.1618  0.0988  0.3286 -0.0452 -0.0953
#> [13,]  0.24425 -0.02905 -0.1879 -0.1107  0.2056 -0.0493  0.1476  0.1851
#> [14,]  0.04798 -0.09932  0.3486  0.2542  0.3893  0.2804 -0.1971 -0.3136
#> [15,]  0.50224 -0.00029 -0.3786  0.0216 -0.0374 -0.3611  0.1582 -0.2058
#>           [,9]    [,10]   [,11]    [,12]   [,13]    [,14]    [,15]
#>  [1,]  0.24433 -0.01812 -0.4597 -0.00257  0.2442  0.04798  0.50224
#>  [2,]  0.00674  0.24073 -0.0307  0.08931 -0.0291 -0.09932 -0.00029
#>  [3,] -0.09669  0.03572  0.1546 -0.22650 -0.1879  0.34858 -0.37860
#>  [4,]  0.15589 -0.08976 -0.3569 -0.16176 -0.1107  0.25419  0.02160
#>  [5,]  0.12201 -0.05591 -0.0416  0.09882  0.2056  0.38927 -0.03737
#>  [6,] -0.08813 -0.45536 -0.4244  0.32864 -0.0493  0.28038 -0.36109
#>  [7,] -0.24713  0.24566  0.1511 -0.04518  0.1476 -0.19708  0.15824
#>  [8,]  0.03301  0.54101 -0.2360 -0.09529  0.1851 -0.31362 -0.20579
#>  [9,]  1.00000  0.00779 -0.0363 -0.33486  0.3231 -0.00254  0.19247
#> [10,]  0.00779  1.00000 -0.0516 -0.36864  0.0724 -0.40959  0.10925
#> [11,] -0.03635 -0.05164  1.0000 -0.27158 -0.0969 -0.02492  0.19063
#> [12,] -0.33486 -0.36864 -0.2716  1.00000 -0.0976 -0.09931 -0.22207
#> [13,]  0.32311  0.07237 -0.0969 -0.09755  1.0000  0.34919  0.31171
#> [14,] -0.00254 -0.40959 -0.0249 -0.09931  0.3492  1.00000 -0.30672
#> [15,]  0.19247  0.10925  0.1906 -0.22207  0.3117 -0.30672  1.00000
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
