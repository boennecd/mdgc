
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
  n_rep <- floor((p + 3 - 1) / 3)
  type <- rep(1:3, each = n_rep)[1:p]
  is_con <- type == 1L
  is_bin <- type == 2L
  is_ord <- type == 3L
  col_nam <- c(outer(1:n_rep, c("C", "B", "O"), 
                     function(x, y) paste0(y, x)))[1:p]
  
  # sample which are masked data 
  is_mask <- matrix(runif(n * p) < .3, n)
  
  # make sure we have no rows with all missing data
  while(any(all_nans <- rowSums(is_mask) == NCOL(is_mask)))
    is_mask[all_nans, ] <- runif(sum(all_nans) * p) < .3
  
  # create observed data
  truth_obs <- data.frame(truth)
  colnames(truth_obs) <- col_nam
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
#>      C1    C2    C3     C4       C5    B1    B2    B3    B4    B5   O1   O2
#> 1 0.237 0.693 0.798 0.0666       NA FALSE FALSE FALSE FALSE  TRUE    E    C
#> 2 0.142    NA    NA 0.0927 0.000152 FALSE    NA  TRUE    NA    NA    E    B
#> 3    NA 0.748 0.629 0.4280       NA    NA  TRUE    NA    NA  TRUE <NA>    A
#> 4 2.702    NA    NA 2.1776 1.700870 FALSE  TRUE  TRUE    NA  TRUE    A <NA>
#> 5 0.925    NA 0.205 0.6046 0.171311  TRUE  TRUE FALSE FALSE FALSE    E    B
#> 6 0.115    NA 1.341     NA       NA FALSE  TRUE  TRUE    NA    NA    E <NA>
#>     O3 O4   O5
#> 1    B  B <NA>
#> 2    A  A    C
#> 3 <NA>  C    E
#> 4    D  B <NA>
#> 5 <NA>  D <NA>
#> 6    A  B <NA>

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
#> 1 Setup time   21.7ms   23.2ms      42.9    9.19MB     13.4

# fit the model using three different methods
set.seed(60941821)
system.time(
  fit_Lagran_start <- mdgc_fit(
    ptr = log_ml_ptr, vcov = start_val, n_threads = 4L, 
    maxit = 25L, method = "aug_Lagran", rel_eps = 1e-3, maxpts = 200L, 
    verbose = FALSE))
#>    user  system elapsed 
#>   141.4     0.0    35.4
system.time(
  fit_Lagran <- mdgc_fit(
    ptr = log_ml_ptr, vcov = fit_Lagran_start$result, n_threads = 4L, 
    maxit = 25L, method = "aug_Lagran", rel_eps = 1e-3, maxpts = 5000L, 
    verbose = FALSE, mu = fit_Lagran_start$mu, 
    lambda = fit_Lagran_start$lambda))
#>    user  system elapsed 
#>   31.99    0.00    8.11

system.time(
  fit_adam <- mdgc_fit(
    ptr = log_ml_ptr, vcov = start_val, n_threads = 4L, 
    lr = 1e-3, maxit = 25L, batch_size = 100L, method = "adam", 
     rel_eps = 1e-3, maxpts = 5000L))
#>    user  system elapsed 
#>   38.99    0.00    9.75

set.seed(fit_seed <- 19570958L)
system.time(
  fit_svrg <- mdgc_fit(
    ptr = log_ml_ptr, vcov = start_val, n_threads = 4L, 
    lr = 1e-3, maxit = 25L, batch_size = 100L, method = "svrg", 
    verbose = TRUE, rel_eps = 1e-3, maxpts = 5000L))
#> End of iteration    1 with learning rate 0.00100000
#> Log marginal likelihood approximation is    -23440.91
#> Previous approximate gradient norm was        3399.64
#> 
#> End of iteration    2 with learning rate 0.00098000
#> Log marginal likelihood approximation is    -23389.86
#> Previous approximate gradient norm was        1701.56
#> 
#> End of iteration    3 with learning rate 0.00096040
#> Log marginal likelihood approximation is    -23367.41
#> Previous approximate gradient norm was        1146.69
#> 
#> End of iteration    4 with learning rate 0.00094119
#> Log marginal likelihood approximation is    -23355.64
#> Previous approximate gradient norm was         844.93
#> 
#> End of iteration    5 with learning rate 0.00092237
#> Log marginal likelihood approximation is    -23348.67
#> Previous approximate gradient norm was         667.57
#> 
#> End of iteration    6 with learning rate 0.00090392
#> Log marginal likelihood approximation is    -23344.30
#> Previous approximate gradient norm was         554.26
#> 
#> End of iteration    7 with learning rate 0.00088584
#> Log marginal likelihood approximation is    -23341.28
#> Previous approximate gradient norm was         474.18
#> 
#> End of iteration    8 with learning rate 0.00086813
#> Log marginal likelihood approximation is    -23339.24
#> Previous approximate gradient norm was         419.90
#> 
#> End of iteration    9 with learning rate 0.00085076
#> Log marginal likelihood approximation is    -23337.69
#> Previous approximate gradient norm was         377.09
#> 
#> End of iteration   10 with learning rate 0.00083375
#> Log marginal likelihood approximation is    -23336.49
#> Previous approximate gradient norm was         343.38
#> 
#> End of iteration   11 with learning rate 0.00081707
#> Log marginal likelihood approximation is    -23335.62
#> Previous approximate gradient norm was         320.76
#> 
#> End of iteration   12 with learning rate 0.00080073
#> Log marginal likelihood approximation is    -23334.89
#> Previous approximate gradient norm was         302.69
#> 
#> End of iteration   13 with learning rate 0.00078472
#> Log marginal likelihood approximation is    -23334.32
#> Previous approximate gradient norm was         288.15
#> 
#> End of iteration   14 with learning rate 0.00076902
#> Log marginal likelihood approximation is    -23333.86
#> Previous approximate gradient norm was         276.06
#> 
#> End of iteration   15 with learning rate 0.00075364
#> Log marginal likelihood approximation is    -23333.45
#> Previous approximate gradient norm was         263.95
#> 
#> End of iteration   16 with learning rate 0.00073857
#> Log marginal likelihood approximation is    -23333.13
#> Previous approximate gradient norm was         253.99
#> 
#> End of iteration   17 with learning rate 0.00072380
#> Log marginal likelihood approximation is    -23332.86
#> Previous approximate gradient norm was         272.13
#> 
#> End of iteration   18 with learning rate 0.00070932
#> Log marginal likelihood approximation is    -23332.65
#> Previous approximate gradient norm was         248.33
#> 
#> End of iteration   19 with learning rate 0.00069514
#> Log marginal likelihood approximation is    -23332.45
#> Previous approximate gradient norm was         243.03
#> 
#> End of iteration   20 with learning rate 0.00068123
#> Log marginal likelihood approximation is    -23332.28
#> Previous approximate gradient norm was         238.29
#> 
#> End of iteration   21 with learning rate 0.00066761
#> Log marginal likelihood approximation is    -23332.31
#> Previous approximate gradient norm was         245.98
#> 
#> End of iteration   22 with learning rate 0.00065426
#> Log marginal likelihood approximation is    -23332.12
#> Previous approximate gradient norm was         234.63
#> 
#> End of iteration   23 with learning rate 0.00064117
#> Log marginal likelihood approximation is    -23331.99
#> Previous approximate gradient norm was         229.97
#> 
#> End of iteration   24 with learning rate 0.00062835
#> Log marginal likelihood approximation is    -23331.72
#> Previous approximate gradient norm was         235.85
#> 
#> End of iteration   25 with learning rate 0.00061578
#> Log marginal likelihood approximation is    -23331.65
#> Previous approximate gradient norm was         227.54
#>    user  system elapsed 
#>    80.1     0.0    20.0

# compare the log marginal likelihood 
print(rbind(
  `Augmented Lagrangian` = 
    mdgc_log_ml(vcov = fit_Lagran$result, ptr = log_ml_ptr, rel_eps = 1e-3),
  ADAM = 
    mdgc_log_ml(vcov = fit_adam$result  , ptr = log_ml_ptr, rel_eps = 1e-3),
  SVRG =
    mdgc_log_ml(vcov = fit_svrg$result  , ptr = log_ml_ptr, rel_eps = 1e-3)), digits = 10)
#>                              [,1]
#> Augmented Lagrangian -23331.10416
#> ADAM                 -23350.31261
#> SVRG                 -23331.61267

# we can use an approximation in the method
set.seed(fit_seed)
system.time(
  fit_svrg_aprx <- mdgc_fit(
    ptr = log_ml_ptr, vcov = start_val, n_threads = 4L, 
    lr = 1e-3, maxit = 25L, batch_size = 100L, method = "svrg", 
    verbose = TRUE, rel_eps = 1e-3, maxpts = 5000L, 
    use_aprx = TRUE))
#> End of iteration    1 with learning rate 0.00100000
#> Log marginal likelihood approximation is    -23440.91
#> Previous approximate gradient norm was        3399.64
#> 
#> End of iteration    2 with learning rate 0.00098000
#> Log marginal likelihood approximation is    -23389.86
#> Previous approximate gradient norm was        1701.56
#> 
#> End of iteration    3 with learning rate 0.00096040
#> Log marginal likelihood approximation is    -23367.41
#> Previous approximate gradient norm was        1146.69
#> 
#> End of iteration    4 with learning rate 0.00094119
#> Log marginal likelihood approximation is    -23355.64
#> Previous approximate gradient norm was         844.93
#> 
#> End of iteration    5 with learning rate 0.00092237
#> Log marginal likelihood approximation is    -23348.67
#> Previous approximate gradient norm was         667.57
#> 
#> End of iteration    6 with learning rate 0.00090392
#> Log marginal likelihood approximation is    -23344.30
#> Previous approximate gradient norm was         554.26
#> 
#> End of iteration    7 with learning rate 0.00088584
#> Log marginal likelihood approximation is    -23341.28
#> Previous approximate gradient norm was         474.18
#> 
#> End of iteration    8 with learning rate 0.00086813
#> Log marginal likelihood approximation is    -23339.24
#> Previous approximate gradient norm was         419.90
#> 
#> End of iteration    9 with learning rate 0.00085076
#> Log marginal likelihood approximation is    -23337.69
#> Previous approximate gradient norm was         377.09
#> 
#> End of iteration   10 with learning rate 0.00083375
#> Log marginal likelihood approximation is    -23336.49
#> Previous approximate gradient norm was         343.38
#> 
#> End of iteration   11 with learning rate 0.00081707
#> Log marginal likelihood approximation is    -23335.62
#> Previous approximate gradient norm was         320.76
#> 
#> End of iteration   12 with learning rate 0.00080073
#> Log marginal likelihood approximation is    -23334.89
#> Previous approximate gradient norm was         302.68
#> 
#> End of iteration   13 with learning rate 0.00078472
#> Log marginal likelihood approximation is    -23334.32
#> Previous approximate gradient norm was         288.15
#> 
#> End of iteration   14 with learning rate 0.00076902
#> Log marginal likelihood approximation is    -23333.86
#> Previous approximate gradient norm was         276.06
#> 
#> End of iteration   15 with learning rate 0.00075364
#> Log marginal likelihood approximation is    -23333.45
#> Previous approximate gradient norm was         263.95
#> 
#> End of iteration   16 with learning rate 0.00073857
#> Log marginal likelihood approximation is    -23333.13
#> Previous approximate gradient norm was         253.99
#> 
#> End of iteration   17 with learning rate 0.00072380
#> Log marginal likelihood approximation is    -23332.86
#> Previous approximate gradient norm was         272.13
#> 
#> End of iteration   18 with learning rate 0.00070932
#> Log marginal likelihood approximation is    -23332.65
#> Previous approximate gradient norm was         248.33
#> 
#> End of iteration   19 with learning rate 0.00069514
#> Log marginal likelihood approximation is    -23332.45
#> Previous approximate gradient norm was         243.03
#> 
#> End of iteration   20 with learning rate 0.00068123
#> Log marginal likelihood approximation is    -23332.28
#> Previous approximate gradient norm was         238.29
#> 
#> End of iteration   21 with learning rate 0.00066761
#> Log marginal likelihood approximation is    -23332.31
#> Previous approximate gradient norm was         245.98
#> 
#> End of iteration   22 with learning rate 0.00065426
#> Log marginal likelihood approximation is    -23332.12
#> Previous approximate gradient norm was         234.63
#> 
#> End of iteration   23 with learning rate 0.00064117
#> Log marginal likelihood approximation is    -23331.99
#> Previous approximate gradient norm was         229.97
#> 
#> End of iteration   24 with learning rate 0.00062835
#> Log marginal likelihood approximation is    -23331.72
#> Previous approximate gradient norm was         235.85
#> 
#> End of iteration   25 with learning rate 0.00061578
#> Log marginal likelihood approximation is    -23331.65
#> Previous approximate gradient norm was         227.54
#>    user  system elapsed 
#>  54.716   0.004  13.682
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

do_plot(fit_Lagran$result, dat$Sigma, "Estimates (Aug. Lagrangian)")
```

<img src="man/figures/README-sim_dat-1.png" width="100%" />

``` r
do_plot(fit_adam  $result, dat$Sigma, "Estimates (ADAM)")
```

<img src="man/figures/README-sim_dat-2.png" width="100%" />

``` r
do_plot(fit_svrg  $result, dat$Sigma, "Estimates (SVRG)")
```

<img src="man/figures/README-sim_dat-3.png" width="100%" />

``` r

norm(fit_Lagran$result - dat$Sigma, "F")
#> [1] 0.486
norm(fit_adam  $result - dat$Sigma, "F")
#> [1] 0.501
norm(fit_svrg  $result - dat$Sigma, "F")
#> [1] 0.488

# perform the imputation
system.time(
  imp_res <- mdgc_impute(mdgc_obj, fit_svrg$result, rel_eps = 1e-3,
                         maxit = 10000L, n_threads = 4L))
#>    user  system elapsed 
#>   16.78    0.00    4.61

# look at the result for one of the observations
imp_res[2L]
#> [[1]]
#> [[1]]$C1
#> [1] 0.142
#> 
#> [[1]]$C2
#> [1] 2.08
#> 
#> [[1]]$C3
#> [1] 0.249
#> 
#> [[1]]$C4
#> [1] 0.0927
#> 
#> [[1]]$C5
#> [1] 0.000152
#> 
#> [[1]]$B1
#> FALSE  TRUE 
#>     1     0 
#> 
#> [[1]]$B2
#> FALSE  TRUE 
#> 0.194 0.806 
#> 
#> [[1]]$B3
#> FALSE  TRUE 
#>     0     1 
#> 
#> [[1]]$B4
#> FALSE  TRUE 
#> 0.805 0.195 
#> 
#> [[1]]$B5
#> FALSE  TRUE 
#> 0.248 0.752 
#> 
#> [[1]]$O1
#> A B C D E 
#> 0 0 0 0 1 
#> 
#> [[1]]$O2
#> A B C D E 
#> 0 1 0 0 0 
#> 
#> [[1]]$O3
#> A B C D E 
#> 1 0 0 0 0 
#> 
#> [[1]]$O4
#> A B C D E 
#> 1 0 0 0 0 
#> 
#> [[1]]$O5
#> A B C D E 
#> 0 0 1 0 0

# compare with the observed and true data
rbind(truth = dat$truth_obs[2L, ], observed = dat$seen_obs[2L, ])
#>             C1   C2    C3     C4       C5    B1   B2   B3    B4   B5 O1 O2 O3
#> truth    0.142 2.63 0.338 0.0927 0.000152 FALSE TRUE TRUE FALSE TRUE  E  B  A
#> observed 0.142   NA    NA 0.0927 0.000152 FALSE   NA TRUE    NA   NA  E  B  A
#>          O4 O5
#> truth     A  C
#> observed  A  C

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
#>      C1    C2    C3     C4       C5    B1    B2    B3    B4    B5 O1 O2 O3 O4
#> 1 0.237 0.693 0.798 0.0666 1.210369 FALSE FALSE FALSE FALSE  TRUE  E  C  B  B
#> 2 0.142 2.085 0.249 0.0927 0.000152 FALSE  TRUE  TRUE FALSE  TRUE  E  B  A  A
#> 3 1.301 0.748 0.629 0.4280 0.571830 FALSE  TRUE  TRUE  TRUE  TRUE  E  A  E  C
#> 4 2.702 0.779 1.137 2.1776 1.700870 FALSE  TRUE  TRUE  TRUE  TRUE  A  B  D  B
#> 5 0.925 0.913 0.205 0.6046 0.171311  TRUE  TRUE FALSE FALSE FALSE  E  B  D  D
#> 6 0.115 1.069 1.341 0.2374 0.233156 FALSE  TRUE  TRUE FALSE  TRUE  E  A  A  B
#>   O5
#> 1  D
#> 2  C
#> 3  E
#> 4  E
#> 5  E
#> 6  B
head(dat$seen_obs)  # observed data
#>      C1    C2    C3     C4       C5    B1    B2    B3    B4    B5   O1   O2
#> 1 0.237 0.693 0.798 0.0666       NA FALSE FALSE FALSE FALSE  TRUE    E    C
#> 2 0.142    NA    NA 0.0927 0.000152 FALSE    NA  TRUE    NA    NA    E    B
#> 3    NA 0.748 0.629 0.4280       NA    NA  TRUE    NA    NA  TRUE <NA>    A
#> 4 2.702    NA    NA 2.1776 1.700870 FALSE  TRUE  TRUE    NA  TRUE    A <NA>
#> 5 0.925    NA 0.205 0.6046 0.171311  TRUE  TRUE FALSE FALSE FALSE    E    B
#> 6 0.115    NA 1.341     NA       NA FALSE  TRUE  TRUE    NA    NA    E <NA>
#>     O3 O4   O5
#> 1    B  B <NA>
#> 2    A  A    C
#> 3 <NA>  C    E
#> 4    D  B <NA>
#> 5 <NA>  D <NA>
#> 6    A  B <NA>
head(dat$truth_obs) # true data
#>      C1    C2    C3     C4       C5    B1    B2    B3    B4    B5 O1 O2 O3 O4
#> 1 0.237 0.693 0.798 0.0666 0.950476 FALSE FALSE FALSE FALSE  TRUE  E  C  B  B
#> 2 0.142 2.630 0.338 0.0927 0.000152 FALSE  TRUE  TRUE FALSE  TRUE  E  B  A  A
#> 3 2.864 0.748 0.629 0.4280 1.341650 FALSE  TRUE  TRUE FALSE  TRUE  C  A  D  C
#> 4 2.702 1.153 0.459 2.1776 1.700870 FALSE  TRUE  TRUE  TRUE  TRUE  A  C  D  B
#> 5 0.925 0.365 0.205 0.6046 0.171311  TRUE  TRUE FALSE FALSE FALSE  E  B  B  D
#> 6 0.115 0.563 1.341 0.7184 0.306274 FALSE  TRUE  TRUE FALSE  TRUE  E  A  A  B
#>   O5
#> 1  E
#> 2  C
#> 3  E
#> 4  E
#> 5  E
#> 6  B

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
#>    B1    B2    B3    B4    B5    O1    O2    O3    O4    O5 
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
#>    C1    C2    C3    C4    C5 
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
#>  46.941   0.024  46.966

# turn binary variables back to logicals
miss_res$ximp[, is_log] <- lapply(
  miss_res$ximp[, is_log], function(x) as.integer(x) > 1L)

rbind(mdgc       = get_classif_error(thresh_dat),
      missForest = get_classif_error(miss_res$ximp))
#>               B1    B2    B3    B4    B5    O1    O2    O3    O4    O5
#> mdgc       0.276 0.292 0.225 0.318 0.288 0.563 0.640 0.623 0.599 0.548
#> missForest 0.315 0.340 0.304 0.371 0.319 0.651 0.726 0.680 0.673 0.612
rbind(mdgc       = get_rmse(thresh_dat),
      missForest = get_rmse(miss_res$ximp))
#>               C1    C2    C3    C4    C5
#> mdgc       0.644 0.783 0.651 0.796 0.746
#> missForest 0.806 0.848 0.755 0.845 0.842
```

#### An Even Shorter Example

Here is an example where we use the `mdgc` function to do the model
estimation and the imputation:

``` r
# have a data set with missing continous, binary, and ordinal variables
head(dat$seen_obs)
#>      C1    C2    C3     C4       C5    B1    B2    B3    B4    B5   O1   O2
#> 1 0.237 0.693 0.798 0.0666       NA FALSE FALSE FALSE FALSE  TRUE    E    C
#> 2 0.142    NA    NA 0.0927 0.000152 FALSE    NA  TRUE    NA    NA    E    B
#> 3    NA 0.748 0.629 0.4280       NA    NA  TRUE    NA    NA  TRUE <NA>    A
#> 4 2.702    NA    NA 2.1776 1.700870 FALSE  TRUE  TRUE    NA  TRUE    A <NA>
#> 5 0.925    NA 0.205 0.6046 0.171311  TRUE  TRUE FALSE FALSE FALSE    E    B
#> 6 0.115    NA 1.341     NA       NA FALSE  TRUE  TRUE    NA    NA    E <NA>
#>     O3 O4   O5
#> 1    B  B <NA>
#> 2    A  A    C
#> 3 <NA>  C    E
#> 4    D  B <NA>
#> 5 <NA>  D <NA>
#> 6    A  B <NA>
  
# perform the estimation and imputation
set.seed(1)
system.time(res <- mdgc(dat$seen_obs, verbose = TRUE, maxpts = 5000L, 
                        n_threads = 4L, maxit = 25L, use_aprx = TRUE))
#> Estimating the model...
#> End of iteration    1 with learning rate 0.00100000
#> Log marginal likelihood approximation is    -23442.63
#> Previous approximate gradient norm was        3393.96
#> 
#> End of iteration    2 with learning rate 0.00098000
#> Log marginal likelihood approximation is    -23390.91
#> Previous approximate gradient norm was        1694.62
#> 
#> End of iteration    3 with learning rate 0.00096040
#> Log marginal likelihood approximation is    -23368.15
#> Previous approximate gradient norm was        1138.76
#> 
#> End of iteration    4 with learning rate 0.00094119
#> Log marginal likelihood approximation is    -23356.07
#> Previous approximate gradient norm was         841.78
#> 
#> End of iteration    5 with learning rate 0.00092237
#> Log marginal likelihood approximation is    -23348.94
#> Previous approximate gradient norm was         658.80
#> 
#> End of iteration    6 with learning rate 0.00090392
#> Log marginal likelihood approximation is    -23344.40
#> Previous approximate gradient norm was         543.07
#> 
#> End of iteration    7 with learning rate 0.00088584
#> Log marginal likelihood approximation is    -23341.44
#> Previous approximate gradient norm was         466.07
#> 
#> End of iteration    8 with learning rate 0.00086813
#> Log marginal likelihood approximation is    -23339.27
#> Previous approximate gradient norm was         411.44
#> 
#> End of iteration    9 with learning rate 0.00085076
#> Log marginal likelihood approximation is    -23337.75
#> Previous approximate gradient norm was         370.94
#> 
#> End of iteration   10 with learning rate 0.00083375
#> Log marginal likelihood approximation is    -23336.53
#> Previous approximate gradient norm was         338.69
#> 
#> End of iteration   11 with learning rate 0.00081707
#> Log marginal likelihood approximation is    -23335.66
#> Previous approximate gradient norm was         316.54
#> 
#> Performing imputation...
#>    user  system elapsed 
#>   15.34    0.00    4.31

# compare the estimated correlation matrix with the truth
norm(dat$Sigma - res$vcov, "F") / norm(dat$Sigma, "F")
#> [1] 0.0956

# compute the classifcation error and RMSE
get_classif_error(res$ximp)
#>    B1    B2    B3    B4    B5    O1    O2    O3    O4    O5 
#> 0.274 0.285 0.226 0.314 0.288 0.574 0.638 0.626 0.605 0.543
get_rmse(res$ximp)
#>    C1    C2    C3    C4    C5 
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
#>               B1    B2    B3    B4    B5    O1    O2    O3    O4    O5
#> mdgc       0.274 0.285 0.226 0.314 0.288 0.574 0.638 0.626 0.605 0.543
#> mixedgcImp 0.281 0.328 0.232 0.320 0.288 0.626 0.694 0.688 0.609 0.556
rbind(mdgc       = get_rmse(res$ximp),
      mixedgcImp = get_rmse(imp_apr_res))
#>               C1    C2    C3    C4    C5
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
log_ml <- function(..., maxpts = 10000L)
  mdgc_log_ml(ptr = log_ml_ptr, ..., maxpts = maxpts)

# print the approximate log marginal likelihood at the true parameters
set.seed(1)
print(log_ml(dat$Sigma, n_threads = 1L), digits = 7)
#> [1] -23382.4
print(log_ml(dat$Sigma, n_threads = 1L, use_aprx = TRUE), digits = 7)
#> [1] -23382.66

# check standard error
sd(replicate(20, log_ml(dat$Sigma, n_threads = 4L)))
#> [1] 0.105
sd(replicate(20, log_ml(dat$Sigma, n_threads = 4L, use_aprx = TRUE)))
#> [1] 0.0988

# without reordering
print(log_ml(dat$Sigma, n_threads = 4L, do_reorder = FALSE), digits = 7)
#> [1] -23383.55

# check standard error
sd(replicate(20, log_ml(dat$Sigma, n_threads = 4L, do_reorder = FALSE)))
#> [1] 0.822

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
#> 1 1 thread                    837.2ms  841.9ms     1.18     33.9KB        0
#> 2 2 threads                  441.62ms  454.6ms     2.21     33.9KB        0
#> 3 4 threads                  226.03ms  238.6ms     4.24     33.9KB        0
#> 4 4 threads (w/ approx)      145.43ms  149.8ms     6.63     33.9KB        0
#> 5 4 threads (w/o reordering)    1.09s     1.1s     0.908    33.9KB        0

#####
# we can also get an approximation of the gradient
t1 <- log_ml(dat$Sigma, comp_derivs = TRUE, n_threads = 1L, rel_eps = 1e-3)
t2 <- log_ml(dat$Sigma, comp_derivs = TRUE, n_threads = 4L, rel_eps = 1e-3)
all.equal(t1, t2, tolerance = 1e-2)
#> [1] "Attributes: < Component \"grad\": Mean relative difference: 0.0127 >"
  
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
#> 1 1 thread                      8.22s    8.23s     0.121    35.7KB        0
#> 2 2 threads                     4.17s     4.2s     0.238    35.7KB        0
#> 3 4 threads                     2.15s    2.17s     0.460    35.7KB        0
#> 4 4 threads (w/ approx)         1.54s    1.56s     0.644    35.7KB        0
#> 5 4 threads (w/o reordering)    2.25s    2.26s     0.436    35.7KB        0

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
#> [1] 0.00136
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
#> 1 mvtnorm      1.08ms   4.45ms      264.    4.43KB     2.01
#> 2 mdgc         2.93ms   7.38ms      150.    2.49KB     0

sd(replicate(25, use_mvtnorm()))
#> [1] 4.63e-09
sd(replicate(25, use_this_pkg()))
#> [1] 4.19e-09

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
#> 331.162   0.052  83.437

# compare estimates with truth
norm(res$result - dat$Sigma)
#> [1] 0.572
res$result
#>           [,1]     [,2]    [,3]    [,4]    [,5]    [,6]    [,7]    [,8]
#>  [1,]  1.00000  0.25541 -0.3508  0.2001 -0.0876  0.1414  0.0226 -0.0037
#>  [2,]  0.25541  1.00000 -0.2899 -0.0152 -0.5044 -0.2880  0.1428  0.1247
#>  [3,] -0.35078 -0.28990  1.0000  0.3396  0.3826 -0.0700  0.2221 -0.0785
#>  [4,]  0.20011 -0.01521  0.3396  1.0000  0.1289  0.1017  0.0204 -0.4636
#>  [5,] -0.08758 -0.50444  0.3826  0.1289  1.0000  0.1684 -0.2664 -0.2146
#>  [6,]  0.14135 -0.28804 -0.0700  0.1017  0.1684  1.0000 -0.1009 -0.2586
#>  [7,]  0.02259  0.14277  0.2221  0.0204 -0.2664 -0.1009  1.0000  0.0813
#>  [8,] -0.00370  0.12474 -0.0785 -0.4636 -0.2146 -0.2586  0.0813  1.0000
#>  [9,]  0.24372  0.00535 -0.0974  0.1548  0.1211 -0.0860 -0.2448  0.0334
#> [10,] -0.01938  0.24037  0.0348 -0.0916 -0.0597 -0.4552  0.2417  0.5380
#> [11,] -0.45695 -0.02519  0.1540 -0.3552 -0.0449 -0.4258  0.1491 -0.2311
#> [12,] -0.00236  0.08877 -0.2270 -0.1631  0.0980  0.3279 -0.0468 -0.1011
#> [13,]  0.24721 -0.02686 -0.1881 -0.1092  0.2038 -0.0516  0.1434  0.1818
#> [14,]  0.04642 -0.10523  0.3488  0.2530  0.3925  0.2805 -0.1925 -0.3119
#> [15,]  0.49841 -0.00680 -0.3775  0.0208 -0.0344 -0.3538  0.1558 -0.1993
#>           [,9]   [,10]   [,11]    [,12]   [,13]    [,14]   [,15]
#>  [1,]  0.24372 -0.0194 -0.4570 -0.00236  0.2472  0.04642  0.4984
#>  [2,]  0.00535  0.2404 -0.0252  0.08877 -0.0269 -0.10523 -0.0068
#>  [3,] -0.09737  0.0348  0.1540 -0.22696 -0.1881  0.34878 -0.3775
#>  [4,]  0.15478 -0.0916 -0.3552 -0.16312 -0.1092  0.25302  0.0208
#>  [5,]  0.12109 -0.0597 -0.0449  0.09798  0.2038  0.39253 -0.0344
#>  [6,] -0.08599 -0.4552 -0.4258  0.32789 -0.0516  0.28054 -0.3538
#>  [7,] -0.24481  0.2417  0.1491 -0.04684  0.1434 -0.19255  0.1558
#>  [8,]  0.03337  0.5380 -0.2311 -0.10108  0.1818 -0.31190 -0.1993
#>  [9,]  1.00000  0.0156 -0.0401 -0.33292  0.3206  0.00349  0.1936
#> [10,]  0.01561  1.0000 -0.0534 -0.36122  0.0749 -0.40744  0.1079
#> [11,] -0.04011 -0.0534  1.0000 -0.27224 -0.0968 -0.02640  0.1917
#> [12,] -0.33292 -0.3612 -0.2722  1.00000 -0.0968 -0.09600 -0.2187
#> [13,]  0.32056  0.0749 -0.0968 -0.09679  1.0000  0.34899  0.3125
#> [14,]  0.00349 -0.4074 -0.0264 -0.09600  0.3490  1.00000 -0.3064
#> [15,]  0.19358  0.1079  0.1917 -0.21868  0.3125 -0.30640  1.0000
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
#> 238.750   0.012  59.938

# compare estimates with the truth
norm(res_adam$result - dat$Sigma)
#> [1] 0.541
res_adam$result
#>           [,1]     [,2]    [,3]    [,4]    [,5]    [,6]    [,7]     [,8]
#>  [1,]  1.00000  0.25374 -0.3548  0.2074 -0.0862  0.1424  0.0160  0.00134
#>  [2,]  0.25374  1.00000 -0.2876 -0.0164 -0.4993 -0.2954  0.1452  0.12079
#>  [3,] -0.35484 -0.28757  1.0000  0.3301  0.3810 -0.0768  0.2314 -0.07170
#>  [4,]  0.20738 -0.01638  0.3301  1.0000  0.1279  0.0991  0.0116 -0.45883
#>  [5,] -0.08623 -0.49926  0.3810  0.1279  1.0000  0.1646 -0.2663 -0.20794
#>  [6,]  0.14243 -0.29540 -0.0768  0.0991  0.1646  1.0000 -0.1071 -0.27543
#>  [7,]  0.01597  0.14520  0.2314  0.0116 -0.2663 -0.1071  1.0000  0.08723
#>  [8,]  0.00134  0.12079 -0.0717 -0.4588 -0.2079 -0.2754  0.0872  1.00000
#>  [9,]  0.24085  0.00705 -0.0954  0.1520  0.1279 -0.0909 -0.2492  0.03657
#> [10,] -0.02309  0.23863  0.0364 -0.0905 -0.0469 -0.4601  0.2506  0.53612
#> [11,] -0.45648 -0.03077  0.1596 -0.3554 -0.0472 -0.4238  0.1635 -0.24678
#> [12,] -0.01055  0.08551 -0.2275 -0.1623  0.0982  0.3241 -0.0492 -0.09207
#> [13,]  0.24712 -0.02764 -0.1897 -0.1038  0.1992 -0.0456  0.1518  0.18637
#> [14,]  0.04414 -0.10199  0.3512  0.2490  0.3973  0.2653 -0.1978 -0.30616
#> [15,]  0.49232 -0.00632 -0.3740  0.0270 -0.0284 -0.3611  0.1592 -0.20542
#>           [,9]   [,10]   [,11]   [,12]   [,13]   [,14]    [,15]
#>  [1,]  0.24085 -0.0231 -0.4565 -0.0105  0.2471  0.0441  0.49232
#>  [2,]  0.00705  0.2386 -0.0308  0.0855 -0.0276 -0.1020 -0.00632
#>  [3,] -0.09536  0.0364  0.1596 -0.2275 -0.1897  0.3512 -0.37401
#>  [4,]  0.15196 -0.0905 -0.3554 -0.1623 -0.1038  0.2490  0.02700
#>  [5,]  0.12788 -0.0469 -0.0472  0.0982  0.1992  0.3973 -0.02837
#>  [6,] -0.09094 -0.4601 -0.4238  0.3241 -0.0456  0.2653 -0.36112
#>  [7,] -0.24924  0.2506  0.1635 -0.0492  0.1518 -0.1978  0.15923
#>  [8,]  0.03657  0.5361 -0.2468 -0.0921  0.1864 -0.3062 -0.20542
#>  [9,]  1.00000  0.0115 -0.0252 -0.3303  0.3228  0.0037  0.19012
#> [10,]  0.01153  1.0000 -0.0477 -0.3654  0.0767 -0.4015  0.12337
#> [11,] -0.02523 -0.0477  1.0000 -0.2759 -0.1029 -0.0156  0.19627
#> [12,] -0.33030 -0.3654 -0.2759  1.0000 -0.0932 -0.1194 -0.23193
#> [13,]  0.32284  0.0767 -0.1029 -0.0932  1.0000  0.3515  0.31539
#> [14,]  0.00370 -0.4015 -0.0156 -0.1194  0.3515  1.0000 -0.30013
#> [15,]  0.19012  0.1234  0.1963 -0.2319  0.3154 -0.3001  1.00000
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
#> [11] -23337 -23336 -23336 -23336 -23335 -23335 -23335 -23334 -23334 -23334
#> [21] -23334 -23334 -23334 -23334 -23334
res_adam$fun_vals # likely lower bounds on the log-marginal likelihood
#>  [1] -25084 -23983 -23556 -23430 -23396 -23384 -23377 -23374 -23372 -23371
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
#>    user  system elapsed 
#> 473.659   0.008 118.880

# compare estimates with the truth
norm(res_svrg$result - dat$Sigma)
#> [1] 0.553
res_svrg$result
#>           [,1]      [,2]    [,3]    [,4]    [,5]    [,6]    [,7]     [,8]
#>  [1,]  1.00000  2.56e-01 -0.3490  0.1981 -0.0891  0.1414  0.0225 -0.00346
#>  [2,]  0.25596  1.00e+00 -0.2860 -0.0193 -0.5061 -0.2898  0.1435  0.12284
#>  [3,] -0.34901 -2.86e-01  1.0000  0.3409  0.3802 -0.0715  0.2229 -0.07799
#>  [4,]  0.19810 -1.93e-02  0.3409  1.0000  0.1296  0.1002  0.0197 -0.46524
#>  [5,] -0.08912 -5.06e-01  0.3802  0.1296  1.0000  0.1686 -0.2674 -0.21404
#>  [6,]  0.14138 -2.90e-01 -0.0715  0.1002  0.1686  1.0000 -0.0951 -0.26978
#>  [7,]  0.02249  1.43e-01  0.2229  0.0197 -0.2674 -0.0951  1.0000  0.08049
#>  [8,] -0.00346  1.23e-01 -0.0780 -0.4652 -0.2140 -0.2698  0.0805  1.00000
#>  [9,]  0.24469  7.01e-03 -0.0970  0.1560  0.1218 -0.0882 -0.2473  0.03350
#> [10,] -0.01820  2.41e-01  0.0356 -0.0904 -0.0565 -0.4552  0.2455  0.54079
#> [11,] -0.45975 -3.07e-02  0.1546 -0.3566 -0.0417 -0.4250  0.1510 -0.23601
#> [12,] -0.00247  8.92e-02 -0.2265 -0.1618  0.0988  0.3289 -0.0453 -0.09591
#> [13,]  0.24428 -2.91e-02 -0.1879 -0.1108  0.2056 -0.0495  0.1472  0.18459
#> [14,]  0.04810 -9.95e-02  0.3484  0.2543  0.3892  0.2803 -0.1969 -0.31373
#> [15,]  0.50215 -2.92e-05 -0.3786  0.0218 -0.0376 -0.3609  0.1586 -0.20526
#>           [,9]    [,10]   [,11]    [,12]   [,13]   [,14]     [,15]
#>  [1,]  0.24469 -0.01820 -0.4598 -0.00247  0.2443  0.0481  5.02e-01
#>  [2,]  0.00701  0.24100 -0.0307  0.08919 -0.0291 -0.0995 -2.92e-05
#>  [3,] -0.09704  0.03562  0.1546 -0.22651 -0.1879  0.3484 -3.79e-01
#>  [4,]  0.15604 -0.09038 -0.3566 -0.16181 -0.1108  0.2543  2.18e-02
#>  [5,]  0.12181 -0.05646 -0.0417  0.09876  0.2056  0.3892 -3.76e-02
#>  [6,] -0.08815 -0.45520 -0.4250  0.32891 -0.0495  0.2803 -3.61e-01
#>  [7,] -0.24731  0.24551  0.1510 -0.04528  0.1472 -0.1969  1.59e-01
#>  [8,]  0.03350  0.54079 -0.2360 -0.09591  0.1846 -0.3137 -2.05e-01
#>  [9,]  1.00000  0.00851 -0.0370 -0.33441  0.3233 -0.0017  1.93e-01
#> [10,]  0.00851  1.00000 -0.0519 -0.36808  0.0725 -0.4095  1.09e-01
#> [11,] -0.03701 -0.05191  1.0000 -0.27169 -0.0967 -0.0251  1.91e-01
#> [12,] -0.33441 -0.36808 -0.2717  1.00000 -0.0976 -0.0990 -2.22e-01
#> [13,]  0.32327  0.07251 -0.0967 -0.09757  1.0000  0.3491  3.12e-01
#> [14,] -0.00170 -0.40953 -0.0251 -0.09901  0.3491  1.0000 -3.06e-01
#> [15,]  0.19277  0.10926  0.1907 -0.22167  0.3116 -0.3064  1.00e+00
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
categorical_arg <- replicate(NCOL(mdgc_obj$lower), matrix(0L, 0, 0), 
                             simplify = FALSE) 
all.equal(tmp, mdgc:::get_z_hat(
  mdgc_obj$lower, mdgc_obj$upper, mdgc_obj$code, 
  categorical = categorical_arg, n_threads = 4L))
#> [1] TRUE

# the latter is faster but both are relatively fast
mark(
  `R version  ` = get_z_hat(mdgc_obj$lower, mdgc_obj$upper, mdgc_obj$code), 
  `C++ verison` = mdgc:::get_z_hat(
  mdgc_obj$lower, mdgc_obj$upper, mdgc_obj$code, 
  categorical = categorical_arg, n_threads = 4L), 
  min_iterations = 10)
#> Warning: Some expressions had a GC in every iteration; so filtering is disabled.
#> # A tibble: 2 x 6
#>   expression       min   median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>  <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 R version     66.9ms   73.1ms      13.4    1020KB     21.4
#> 2 C++ verison  279.5µs  357.4µs    2609.      234KB     12.0

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
#> 245.015   0.024  61.550

# for comparisons, we also run the code using one thread
set.seed(1)
system.time(res_adam_ser  <- adam(
  val = start_val, alpha = 1e-2, maxit = 25L, batch_size = 100L, 
  n_threads = 1L))
#>    user  system elapsed 
#> 204.228   0.004 204.233

# we get (roughly) the same
norm(res_adam$result - res_adam_ser$result)
#> [1] 0.00605

# plot estimate
norm(res_adam$result - dat$Sigma)
#> [1] 0.538
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
#>  [1] -23395 -23358 -23347 -23342 -23339 -23337 -23337 -23336 -23335 -23335
#> [11] -23335 -23334 -23334 -23334 -23334 -23334 -23334 -23334 -23333 -23333
#> [21] -23333 -23333 -23333 -23333 -23333
res_adam$fun_vals # likely lower bounds on the log-marginal likelihood
#>  [1] -23499 -23397 -23383 -23374 -23372 -23370 -23368 -23368 -23367 -23367
#> [11] -23366 -23366 -23366 -23366 -23366 -23365 -23365 -23365 -23365 -23365
#> [21] -23365 -23365 -23365 -23365 -23365

# do the same with SVRG
set.seed(1)
system.time(res_svrg  <- svrg(
  val = start_val, lr = 1e-3, maxit = 25L, batch_size = 100L))
#>    user  system elapsed 
#> 480.017   0.012 120.478

# compare estimates with the truth
norm(res_svrg$result - dat$Sigma)
#> [1] 0.551
res_svrg$result
#>           [,1]      [,2]    [,3]    [,4]    [,5]    [,6]    [,7]     [,8]
#>  [1,]  1.00000  0.255942 -0.3490  0.1981 -0.0892  0.1413  0.0225 -0.00349
#>  [2,]  0.25594  1.000000 -0.2860 -0.0193 -0.5061 -0.2901  0.1436  0.12261
#>  [3,] -0.34896 -0.285982  1.0000  0.3409  0.3802 -0.0716  0.2229 -0.07799
#>  [4,]  0.19807 -0.019318  0.3409  1.0000  0.1296  0.0998  0.0197 -0.46538
#>  [5,] -0.08919 -0.506125  0.3802  0.1296  1.0000  0.1686 -0.2675 -0.21396
#>  [6,]  0.14134 -0.290140 -0.0716  0.0998  0.1686  1.0000 -0.0944 -0.27126
#>  [7,]  0.02249  0.143620  0.2229  0.0197 -0.2675 -0.0944  1.0000  0.08066
#>  [8,] -0.00349  0.122610 -0.0780 -0.4654 -0.2140 -0.2713  0.0807  1.00000
#>  [9,]  0.24474  0.007078 -0.0969  0.1562  0.1221 -0.0883 -0.2477  0.03334
#> [10,] -0.01806  0.241096  0.0356 -0.0901 -0.0561 -0.4550  0.2456  0.54084
#> [11,] -0.45978 -0.030665  0.1546 -0.3567 -0.0416 -0.4251  0.1513 -0.23653
#> [12,] -0.00252  0.089133 -0.2265 -0.1616  0.0989  0.3289 -0.0452 -0.09537
#> [13,]  0.24414 -0.029101 -0.1879 -0.1109  0.2056 -0.0491  0.1477  0.18505
#> [14,]  0.04800 -0.099333  0.3484  0.2543  0.3892  0.2804 -0.1974 -0.31380
#> [15,]  0.50217 -0.000161 -0.3785  0.0218 -0.0376 -0.3615  0.1587 -0.20603
#>           [,9]    [,10]   [,11]    [,12]   [,13]    [,14]     [,15]
#>  [1,]  0.24474 -0.01806 -0.4598 -0.00252  0.2441  0.04800  0.502175
#>  [2,]  0.00708  0.24110 -0.0307  0.08913 -0.0291 -0.09933 -0.000161
#>  [3,] -0.09694  0.03565  0.1546 -0.22647 -0.1879  0.34840 -0.378530
#>  [4,]  0.15620 -0.09013 -0.3567 -0.16163 -0.1109  0.25429  0.021759
#>  [5,]  0.12207 -0.05611 -0.0416  0.09889  0.2056  0.38919 -0.037601
#>  [6,] -0.08829 -0.45504 -0.4251  0.32885 -0.0491  0.28037 -0.361486
#>  [7,] -0.24766  0.24562  0.1513 -0.04519  0.1477 -0.19743  0.158680
#>  [8,]  0.03334  0.54084 -0.2365 -0.09537  0.1850 -0.31380 -0.206029
#>  [9,]  1.00000  0.00751 -0.0365 -0.33475  0.3238 -0.00247  0.192494
#> [10,]  0.00751  1.00000 -0.0516 -0.36883  0.0724 -0.40994  0.109233
#> [11,] -0.03652 -0.05159  1.0000 -0.27159 -0.0968 -0.02493  0.190749
#> [12,] -0.33475 -0.36883 -0.2716  1.00000 -0.0976 -0.09926 -0.221956
#> [13,]  0.32376  0.07240 -0.0968 -0.09761  1.0000  0.34916  0.311786
#> [14,] -0.00247 -0.40994 -0.0249 -0.09926  0.3492  1.00000 -0.306548
#> [15,]  0.19249  0.10923  0.1907 -0.22196  0.3118 -0.30655  1.000000
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

## Adding Multinomial Variables

We extend the model suggested by Zhao and Udell (2019) in this section.
The example is very similar to the previous one but with multinomial
variables.

``` r
# simulates a data set and mask some of the data.
#
# Args: 
#   n: number of observations. 
#   p: number of variables. 
#   n_lvls: number of levels for the ordinal and multinomial variables. 
#   verbose: print status during the simulation.
# 
# Returns: 
#   Simluated masked data, the true data, and true covariance matrix. 
sim_dat <- function(n, p = 4L, n_lvls = 5L, verbose = FALSE){
  # determine the type
  n_rep <- floor((p + 4 - 1) / 4)
  type <- rep(1:4, n_rep)[1:p]
  is_con  <- type == 1L
  is_bin  <- type == 2L
  is_ord  <- type == 3L
  is_mult <- type == 4L
  
  col_nam <- c(outer(c("C", "B", "O", "M"), 1:n_rep, paste0))[1:p]
  idx <- head(cumsum(c(1L, ifelse(type == 4, n_lvls, 1L))), -1L)
  
  # get the covariance matrix
  n_latent <- p + (n_lvls - 1L) * (p %/% 4)

  Sig <- cov2cor(drop(rWishart(1L, 2 * n_latent, diag(n_latent))))
  
  # Find the closest matrix which satisfy the zero constraint
  type_augmented <- unlist(sapply(seq_along(type), function(x)
    if(x %% 4 == 0) rep(x, n_lvls) else x))
  zero_indices <- do.call(rbind, tapply(
    seq_along(type_augmented)[type_augmented %% 4L == 0L], 
    type_augmented[type_augmented %% 4L == 0L], function(x)
      subset(expand.grid(row = x, col = x), row > col)))
  zero_indices <- as.matrix(zero_indices) - 1L
  
  # setup objects and functions to run augmented Lagrangian method
  par <- t(chol(Sig))
  par <- par[lower.tri(par, TRUE)]
  
  # you can skip these function definitions
  com_vec <- mdgc:::get_commutation_vec(n_latent, n_latent, FALSE)
  get_mat_from_cholesky <- function(x){
    out <- matrix(0., n_latent, n_latent)
    out[lower.tri(out, TRUE)] <- x
    tcrossprod(out)
  }
  fn <- function(x, lambda, mu){
    sig_hat <- get_mat_from_cholesky(x)
    consts <- mdgc:::lower_tri_inner(
      x = x, idx = zero_indices, jacob = FALSE, rhs = numeric())
    
    norm(sig_hat - Sig, "F")^2 / 2 - sum(lambda * consts) +
      mu / 2 * sum(consts^2)
  }
  gr <- function(x, lambda, mu){
    sig_hat <- get_mat_from_cholesky(x)
    
    L <- matrix(0., n_latent, n_latent)
    L[lower.tri(L, TRUE)] <- x
    
    d_Sig <- sig_hat - Sig
    d_Sig <- ((d_Sig + d_Sig[com_vec]) %*% L)[lower.tri(L, TRUE)]
    
    consts <- mdgc:::lower_tri_inner(
      x = x, idx = zero_indices, jacob = FALSE, rhs = numeric())
    rhs <- -lambda + mu * consts
    d_aug <- mdgc:::lower_tri_inner(
      x = x, idx = zero_indices, jacob = TRUE, rhs = rhs)
    
    d_Sig + d_aug
  }
  
  # augmented Lagrangian method
  lambda <- numeric(NROW(zero_indices))
  mu <- 1
  
  # check gradient
  if(FALSE){
    x_test      <- runif(length(par))
    lambda_test <- runif(length(lambda))
    all.equal(numDeriv::grad(fn, x_test, lambda = lambda_test, mu = mu),
              gr(x_test, lambda = lambda_test, mu = mu))
  }
  
  for(i in 1:25){
    opt <- optim(par, fn, gr, mu = mu, lambda = lambda, method = "BFGS")
    par <- opt$par
    
    consts <- mdgc:::lower_tri_inner(
      x = par, idx = zero_indices, jacob = FALSE, rhs = numeric())
      
    all_ok <- all(abs(consts) < 1e-8)
    if(all_ok)
      break

    # update lambda, mu and par. Then repaet
    lambda <- lambda - mu * consts
    mu <- mu * 5
  }
  
  # we found a matrix. We set it
  Sig_new <- cov2cor(get_mat_from_cholesky(par))
  if(verbose)
    cat(sprintf(
      "Relative norm difference between the constrained and simulated covariance matrix is %.2f\n", 
      norm(Sig - Sig_new, "F") / norm(Sig, "F")))
  
  Sig <- Sig_new
      
  # draw the observations
  truth <- matrix(rnorm(n * n_latent), n) %*% chol(Sig)
  
  # sample which are masked data 
  is_mask <- matrix(runif(n * p) < .3, n)
  
  # make sure we have no rows with all missing data
  while(any(all_nans <- rowSums(is_mask) == NCOL(is_mask)))
    is_mask[all_nans, ] <- runif(sum(all_nans) * p) < .3
  
  # create the observed data
  truth_obs <- lapply(type, function(i) if(i == 1L) numeric(n) else integer(n))
  truth_obs <- data.frame(truth_obs)
  colnames(truth_obs) <- col_nam
  
  bs_ord <- qnorm(seq(0, 1, length.out = n_lvls + 1L))
  for(i in 1:p){
    idx_i <- idx[i]
    switch(
      type[i],
      # continous
      truth_obs[, i] <- qexp(pnorm(truth[, idx_i])),
      # binary
      truth_obs[, i] <- truth[, idx_i] > 0,
      # ordinal
      {
        truth_obs[, i] <- 
          ordered(as.integer(cut(truth[, idx_i], breaks = bs_ord)))
        levels(truth_obs[, i]) <- 
          LETTERS[seq_len(length(unique(truth_obs[, i])))]
      },
      # multinomial
      {
        truth_obs[, i] <- apply(
          truth[, idx_i + 1:n_lvls - 1L], 1L, which.max)
        truth_obs[, i] <- factor(truth_obs[, i], 
                                 labels = paste0("T", 1:n_lvls))
      }, 
      stop("Type is not implemented"))
  }

  # mask the data
  seen_obs <- truth_obs
  seen_obs[is_mask] <- NA
  
  list(truth = truth, truth_obs = truth_obs, seen_obs = seen_obs, 
       Sigma = Sig)
}

# simulate and show the data
set.seed(1)
p <- 8L
dat <- sim_dat(2000L, p = p, verbose = TRUE, n_lvls = 3)
#> Relative norm difference between the constrained and simulated covariance matrix is 0.15

# show the first rows of the observed data
head(dat$seen_obs)
#>      C1    B1   O1 M1    C2    B2   O2   M2
#> 1 2.255  TRUE <NA> T1    NA  TRUE    C <NA>
#> 2 1.514    NA    B T2 0.722  TRUE    C   T2
#> 3 0.533 FALSE    C T2 0.612  TRUE    B   T1
#> 4 0.409  TRUE    A T1    NA FALSE <NA>   T1
#> 5    NA  TRUE    A T3    NA FALSE    C   T2
#> 6 1.836 FALSE    C T1 1.973 FALSE <NA>   T3

# assign object to perform the estimation and the imputation
obj <- get_mdgc(dat$seen_obs)
ptr <- get_mdgc_log_ml(obj)

# get starting values
start_vals <- mdgc_start_value(obj)

# plot the starting values and the true values (should not match because of
# overparameterization)
par(mar = c(3, 3, 2, 1), mfcol = c(1, 2))
sc <- colorRampPalette(c("Red", "White", "Blue"))(201)
image(start_vals [, NCOL(dat$Sigma):1], zlim = c(-1, 1), col = sc, 
      main = "Starting values")
image(dat$Sigma  [, NCOL(dat$Sigma):1], zlim = c(-1, 1), col = sc,
      main = "Truth")
```

<img src="man/figures/README-mult_sim-1.png" width="100%" />

``` r
# check the log marginal likelihood at the starting values and compare with
# the true values 
mdgc_log_ml(ptr, start_vals, n_threads = 1L) # at the starting values
#> [1] -11842
mdgc_log_ml(ptr, dat$Sigma , n_threads = 1L) # at the true values
#> [1] -11856

# much better than using a diagonal matrix!
mdgc_log_ml(ptr, diag(NROW(dat$Sigma)), n_threads = 1L)
#> [1] -12085

# estimate the model
system.time(
  ests <- mdgc_fit(ptr, vcov = start_vals, method = "aug_Lagran",
                   n_threads = 4L, rel_eps = 1e-2, maxpts = 1000L, 
                   minvls = 200L, use_aprx = TRUE))
#>    user  system elapsed 
#>    44.5     0.0    11.1

# refine the estimates
system.time(
  ests <- mdgc_fit(ptr, vcov = ests$result, method = "aug_Lagran",
                   n_threads = 4L, rel_eps = 1e-3, maxpts = 10000L, 
                   minvls = 1000L, mu = ests$mu, lambda = ests$lambda, 
                   use_aprx = TRUE))
#>    user  system elapsed 
#>   19.25    0.00    4.89

# compare the estimated and the true values (should not match because of
# overparameterization)
image(ests$result[, NCOL(dat$Sigma):1], zlim = c(-1, 1), col = sc, 
      main = "Estimates")
image(dat$Sigma  [, NCOL(dat$Sigma):1], zlim = c(-1, 1), col = sc, 
      main = "Truth")
```

<img src="man/figures/README-mult_sim-2.png" width="100%" />

``` r
# perform the imputation
system.time(
  imp_res <- mdgc_impute(obj, ests$result, rel_eps = 1e-3,
                         maxit = 10000L, n_threads = 4L))
#>    user  system elapsed 
#>  11.237   0.008   2.984

# look at the result for one of the observations
imp_res[1L]
#> [[1]]
#> [[1]]$C1
#> [1] 2.25
#> 
#> [[1]]$B1
#> FALSE  TRUE 
#>     0     1 
#> 
#> [[1]]$O1
#>     A     B     C 
#> 0.456 0.331 0.213 
#> 
#> [[1]]$M1
#> T1 T2 T3 
#>  1  0  0 
#> 
#> [[1]]$C2
#> [1] 0.749
#> 
#> [[1]]$B2
#> FALSE  TRUE 
#>     0     1 
#> 
#> [[1]]$O2
#> A B C 
#> 0 0 1 
#> 
#> [[1]]$M2
#>    T1    T2    T3 
#> 0.571 0.127 0.301

# compare with the observed and true data
rbind(truth = dat$truth_obs[1L, ], observed = dat$seen_obs[1L, ])
#>            C1   B1   O1 M1    C2   B2 O2   M2
#> truth    2.25 TRUE    C T1 0.505 TRUE  C   T3
#> observed 2.25 TRUE <NA> T1    NA TRUE  C <NA>

# we can threshold the data like this
threshold <- function(org_data, imputed){
  # checks
  stopifnot(NROW(org_data) == length(imputed), 
            is.list(imputed), is.data.frame(org_data))
  
  # threshold
  is_cont <- which(sapply(org_data, is.numeric))
  is_bin  <- which(sapply(org_data, is.logical)) 
  is_ord  <- which(sapply(org_data, is.ordered))
  is_mult <- which(sapply(org_data, is.factor))
  is_mult <- setdiff(is_mult, is_ord)
  stopifnot(
    length(is_cont) + length(is_bin) + length(is_ord) + length(is_mult) == 
      NCOL(org_data))
  is_cat <- c(is_bin, is_ord, is_mult)
  
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
  out <- out[, order(c(is_cont, is_bin, is_ord, is_mult))]
  if(length(is_bin) > 0)
    out[, is_bin] <- out[, is_bin] > 1L
  if(length(is_ord) > 0)
    for(i in is_ord)
      out[[i]] <- ordered(
        unlist(out[[i]]), labels = levels(org_data[, i]))
  if(length(is_mult) > 0)
    for(i in is_mult)
      out[[i]] <- factor(
        unlist(out[[i]]), labels = levels(org_data[, i]))
  
  colnames(out) <- colnames(org_data)
  out
}
thresh_dat <- threshold(dat$seen_obs, imp_res)

# compare thresholded data with observed and true data
head(thresh_dat)
#>      C1    B1 O1 M1    C2    B2 O2 M2
#> 1 2.255  TRUE  A T1 0.749  TRUE  C T1
#> 2 1.514 FALSE  B T2 0.722  TRUE  C T2
#> 3 0.533 FALSE  C T2 0.612  TRUE  B T1
#> 4 0.409  TRUE  A T1 0.737 FALSE  A T1
#> 5 1.045  TRUE  A T3 0.511 FALSE  C T2
#> 6 1.836 FALSE  C T1 1.973 FALSE  B T3
head(dat$seen_obs)  # observed data
#>      C1    B1   O1 M1    C2    B2   O2   M2
#> 1 2.255  TRUE <NA> T1    NA  TRUE    C <NA>
#> 2 1.514    NA    B T2 0.722  TRUE    C   T2
#> 3 0.533 FALSE    C T2 0.612  TRUE    B   T1
#> 4 0.409  TRUE    A T1    NA FALSE <NA>   T1
#> 5    NA  TRUE    A T3    NA FALSE    C   T2
#> 6 1.836 FALSE    C T1 1.973 FALSE <NA>   T3
head(dat$truth_obs) # true data
#>      C1    B1 O1 M1    C2    B2 O2 M2
#> 1 2.255  TRUE  C T1 0.505  TRUE  C T3
#> 2 1.514  TRUE  B T2 0.722  TRUE  C T2
#> 3 0.533 FALSE  C T2 0.612  TRUE  B T1
#> 4 0.409  TRUE  A T1 0.298 FALSE  B T1
#> 5 0.412  TRUE  A T3 0.648 FALSE  C T2
#> 6 1.836 FALSE  C T1 1.973 FALSE  B T3

# compare correct categories
get_classif_error <- function(impu_dat, truth = dat$truth_obs, 
                              observed = dat$seen_obs){
  is_cat <- sapply(truth, function(x)
    is.logical(x) || is.factor(x))
  is_match <- impu_dat[, is_cat] == truth[, is_cat]
  is_match <- matrix(is_match, ncol = sum(is_cat))
  is_match[!is.na(observed[, is_cat])] <- NA_integer_
  setNames(1 - colMeans(is_match, na.rm = TRUE), 
           colnames(truth)[is_cat])
}
get_classif_error(thresh_dat)
#>    B1    O1    M1    B2    O2    M2 
#> 0.410 0.544 0.591 0.385 0.604 0.540

# compute RMSE
get_rmse <- function(impu_dat, truth = dat$truth_obs,
                     observed = dat$seen_obs){
  is_con <- sapply(truth, is.numeric)
  err <- as.matrix(impu_dat[, is_con] - truth[, is_con])
  err[!is.na(observed[, is_con])] <- NA_real_
  sqrt(colMeans(err^2, na.rm = TRUE))
}
get_rmse(thresh_dat)
#>    C1    C2 
#> 1.021 0.965

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
#>   missForest iteration 10 in progress...done!
#>    user  system elapsed 
#>   9.733   0.036   9.769

# turn binary variables back to logicals
miss_res$ximp[, is_log] <- lapply(
  miss_res$ximp[, is_log], function(x) as.integer(x) > 1L)

# compare errors
rbind(mdgc       = get_classif_error(thresh_dat),
      missForest = get_classif_error(miss_res$ximp))
#>              B1    O1    M1    B2    O2    M2
#> mdgc       0.41 0.544 0.591 0.385 0.604 0.540
#> missForest 0.47 0.607 0.602 0.417 0.639 0.558
rbind(mdgc       = get_rmse(thresh_dat),
      missForest = get_rmse(miss_res$ximp))
#>              C1    C2
#> mdgc       1.02 0.965
#> missForest 1.06 0.966
```

### Edgar Anderson’s Iris Data

We make a small example below were we take the iris data set and
randomly mask it. Then we compare the imputation method in this package
with missForest.

``` r
# load the iris data set
data(iris)

# assign function to produce iris data set with NAs.
# 
# Args:
#   p_na: chance of missing a value.
get_iris <- function(p_na = .2){
  is_miss <- matrix(p_na > runif(NROW(iris) * NCOL(iris)), 
                    NROW(iris), NCOL(iris))
  while(any(all_missing <- apply(is_miss, 1, all)))
    # avoid rows with all missing variables
    is_miss[all_missing, ] <- p_na > runif(sum(all_missing) * NCOL(iris))
  
  # create data set with missing values
  out <- iris
  out[is_miss] <- NA 
  out
}

# get a data set with all missing values
set.seed(68129371)
dat <- get_iris()

# use the mdgc method
system.time(
  mdgc_res <- mdgc(dat, maxpts = 5000L, n_threads = 4L, maxit = 25L, 
                   use_aprx = TRUE, method = "aug_Lagran"))
#>    user  system elapsed 
#>    5.41    0.00    1.37

# some of the impuated values
head(mdgc_res$ximp)
#>   Sepal.Length Sepal.Width Petal.Length Petal.Width Species
#> 1          5.1         3.5          1.4         0.2  setosa
#> 2          4.9         3.0          1.4         0.2  setosa
#> 3          4.7         3.2          1.3         0.2  setosa
#> 4          4.6         3.1          1.5         0.2  setosa
#> 5          5.1         3.6          1.4         0.2  setosa
#> 6          5.8         3.9          1.7         1.1  setosa

# the errors
get_classif_error(impu_dat = mdgc_res$ximp, truth = iris, 
                  observed = dat)
#> Species 
#>   0.244
get_rmse(impu_dat = mdgc_res$ximp, truth = iris, observed = dat)
#> Sepal.Length  Sepal.Width Petal.Length  Petal.Width 
#>        0.419        0.340        0.291        0.305

# compare with missForest
system.time(miss_res <- missForest(dat))
#>   missForest iteration 1 in progress...done!
#>   missForest iteration 2 in progress...done!
#>   missForest iteration 3 in progress...done!
#>   missForest iteration 4 in progress...done!
#>   missForest iteration 5 in progress...done!
#>   missForest iteration 6 in progress...done!
#>    user  system elapsed 
#>   0.365   0.004   0.369

# the errors
get_classif_error(impu_dat = miss_res$ximp, truth = iris, 
                  observed = dat)
#> Species 
#>   0.122
get_rmse(impu_dat = miss_res$ximp, truth = iris, observed = dat)
#> Sepal.Length  Sepal.Width Petal.Length  Petal.Width 
#>        0.321        0.270        0.362        0.249
```

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
