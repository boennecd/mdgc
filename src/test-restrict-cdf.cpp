#include "restrict-cdf.h"
#include <testthat.h>
#include <limits>
#include "threat-safe-random.h"

context("restrictcdf unit tests") {
  test_that("cdf<likelihood> gives similar output to R") {
/*
 set.seed(1)
 n <- 4
 mean <- rnorm(n)
 sigma <- drop(rWishart(1L, 2L * n, diag(n)))

 lower <- rep(-Inf, n)
 upper <- rep(0, n)
 mean  <- round(mean , 3)
 sigma <- round(sigma, 3)

 library(mvtnorm)
 prob <- pmvnorm(lower, upper, mean, sigma = sigma,
 algorithm = GenzBretz(abseps = 1e-9, maxpts = 1000000L))

 dput(mean)
 dput(sigma)
 dput(prob)
*/
    Rcpp::RNGScope rngScope;
    parallelrng::set_rng_seeds(1L);

    arma::vec mean;
    arma::mat sigma;

    mean << -0.626 << 0.18 << -0.836 << 1.595;

    sigma << 8.287 << -0.848 << -0.879 << -1.788 << -0.848 << 3.581
          << 2.916 << -3.957 << -0.879 << 2.916 << 7.361 << -0.648
          << -1.788 << -3.957 << -0.648 << 11.735;
    sigma.reshape(4L, 4L);
    restrictcdf::cdf<restrictcdf::likelihood>::set_working_memory(4L, 1L);

    double const abseps = std::pow(std::numeric_limits<double>::epsilon(),
                                   .33);
    double constexpr E_prop(0.0181507102495727);
    {
      auto res = restrictcdf::cdf<restrictcdf::likelihood>(
        mean, sigma).approximate(1000000L, abseps, -1);

      expect_true(res.inform == 0L);
      expect_true(res.abserr                        < 100. * abseps);
      expect_true(std::abs(res.finest[0L] - E_prop) < 100. * abseps);
    }
  }

  test_that("cdf<likelihood> gives similar output to R (with bounds)") {
    /*
     set.seed(1)
     n <- 4
     mean <- rnorm(n)
     sigma <- drop(rWishart(1L, 2L * n, diag(n)))

     lower <- c(-Inf, -1 ,  1, -Inf)
     upper <- c(   1, Inf,  3,    2)
     mean  <- round(mean , 3)
     sigma <- round(sigma, 3)

     library(mvtnorm)
     prob <- pmvnorm(lower, upper, mean, sigma = sigma,
     algorithm = GenzBretz(abseps = 1e-9, maxpts = 1000000L))

     dput(mean)
     dput(sigma)
     dput(prob)
     */
    Rcpp::RNGScope rngScope;
    parallelrng::set_rng_seeds(1L);
    constexpr double const Inf = std::numeric_limits<double>::infinity();

    arma::vec lower, upper, mean;
    arma::mat sigma;

    lower << -Inf   << -1   << 1      << -Inf;
    upper << 1      << Inf  << 3      << 2;
    mean  << -0.626 << 0.18 << -0.836 << 1.595;

    sigma << 8.287 << -0.848 << -0.879 << -1.788 << -0.848 << 3.581
          << 2.916 << -3.957 << -0.879 << 2.916 << 7.361 << -0.648
          << -1.788 << -3.957 << -0.648 << 11.735;
    sigma.reshape(4L, 4L);
    restrictcdf::cdf<restrictcdf::likelihood>::set_working_memory(4L, 1L);

    double const abseps = std::pow(std::numeric_limits<double>::epsilon(),
                                   .33);
    double constexpr E_prop(0.0693596863013216);
    {
      auto res = restrictcdf::cdf<restrictcdf::likelihood>(
        lower, upper, mean, sigma, false).approximate(1000000L, abseps, -1);

      expect_true(res.inform == 0L);
      expect_true(res.abserr                        < 100. * abseps);
      expect_true(std::abs(res.finest[0L] - E_prop) < 100. * abseps);
    }
  }

  test_that("cdf<deriv> gives similar output to R") {
    /*
     set.seed(2)
     n <- 4
    mean <- rnorm(n)
    sigma <- drop(rWishart(1L, 2L * n, diag(n)))

    lower <- rep(-Inf, n)
    upper <- rep(0, n)
    mean  <- round(mean , 3)
    sigma <- round(sigma, 3)

    library(mvtnorm)
    f <- function(par = c(mean, sigma[upper.tri(sigma, TRUE)])){
    set.seed(1)
    mean <- par[1:n]
    sigma[upper.tri(sigma, TRUE)] <- par[-(1:n)]
    sigma[lower.tri(sigma)] <- t(sigma)[lower.tri(sigma)]
    pmvnorm(lower, upper, mean, sigma = sigma,
            algorithm = GenzBretz(abseps = 1e-9, maxpts = 10000L))
    }

    dput(mean)
    dput(sigma)
    library(numDeriv)
    f_val <- f()
    jac <- drop(jacobian(f, c(mean, sigma[upper.tri(sigma, TRUE)])))
    idx <- n + setdiff(1:((n * (n + 1))/ 2L), cumsum(1:n))
    jac[idx] <- jac[idx] / 2

    dput(c(f_val, jac))
     */
    Rcpp::RNGScope rngScope;
    parallelrng::set_rng_seeds(1L);

    arma::vec mean;
    arma::mat sigma;

    mean << -0.897 << 0.185 << 1.588 << -1.13;

    sigma << 6.703 << -0.621 << -0.359 << -1.017 << -0.621 << 3.85 << 0.847
          << -1.931 << -0.359 << 0.847 << 13.438 << 6.106 << -1.017
          << -1.931 << 6.106 << 11.67;
    sigma.reshape(4L, 4L);
    restrictcdf::cdf<restrictcdf::deriv>::set_working_memory(4L, 1L);

    arma::vec expect;
    expect << 0.0724695784076199 << -0.0198615541220946 << -0.0348869861554662
           << -0.0165673832643242 << -0.01129675889845 << -0.00060013043408606
           << 0.00434208509359746 << 0.00150139459997432 << 0.00229746814495135
           << 0.00453344700742309 << -0.000140031793841558 << 0.00134120630703679
           << 0.00191440588268495 << 0.00196875385020239 << -0.00114337234043834;

    double const abseps = std::pow(
      std::numeric_limits<double>::epsilon(), .4);
    auto res = restrictcdf::cdf<restrictcdf::deriv>(
      mean, sigma).approximate(1000000L, abseps, -1);

    expect_true(res.inform == 0L);
    expect_true(res.finest.n_elem == expect.n_elem);
    for(unsigned i = 0; i < expect.n_elem; ++i)
      expect_true(std::abs(res.finest[i] - expect[i]) < 1e-5);
  }

  test_that("cdf<deriv> gives similar output to R (with bounds)") {
    /*
     set.seed(2)
     n <- 4
     mean <- rnorm(n)
     sigma <- drop(rWishart(1L, 2L * n, diag(n)))

     lower <- c(-Inf, -1 ,  1, -Inf)
     upper <- c(   1, Inf,  3,    2)
     mean  <- round(mean , 3)
     sigma <- round(sigma, 3)

     library(mvtnorm)
     f <- function(par = c(mean, sigma[upper.tri(sigma, TRUE)])){
     set.seed(1)
     mean <- par[1:n]
     sigma[upper.tri(sigma, TRUE)] <- par[-(1:n)]
     sigma[lower.tri(sigma)] <- t(sigma)[lower.tri(sigma)]
     pmvnorm(lower, upper, mean, sigma = sigma,
     algorithm = GenzBretz(abseps = 1e-9, maxpts = 10000L))
     }

     dput(mean)
     dput(sigma)
     library(numDeriv)
     f_val <- f()
     jac <- drop(jacobian(f, c(mean, sigma[upper.tri(sigma, TRUE)])))
     idx <- n + setdiff(1:((n * (n + 1))/ 2L), cumsum(1:n))
     jac[idx] <- jac[idx] / 2

     dput(c(f_val, jac))
     */
    Rcpp::RNGScope rngScope;
    parallelrng::set_rng_seeds(1L);

    arma::vec lower, upper, mean;
    arma::mat sigma;
    constexpr double const Inf = std::numeric_limits<double>::infinity();

    lower << -Inf  << -1    << 1     << -Inf;
    upper << 1     << Inf   << 3     << 2;
    mean << -0.897 << 0.185 << 1.588 << -1.13;

    sigma << 6.703 << -0.621 << -0.359 << -1.017 << -0.621 << 3.85 << 0.847
          << -1.931 << -0.359 << 0.847 << 13.438 << 6.106 << -1.017
          << -1.931 << 6.106 << 11.67;
    sigma.reshape(4L, 4L);
    restrictcdf::cdf<restrictcdf::deriv>::set_working_memory(4L, 1L);

    arma::vec expect;
    expect << 0.10798294314481     << -0.016084833560177   << 0.0207517208958157
           << 0.00542773417060216  << -0.00878426310906256 << -0.00239199279704842
           << -0.00184810009109447 << -0.00451083625704058 << -0.000415749485688216
           << 0.00120569609620797  << -0.00426153807790556 << 0.000511171228139317
           << -0.00150307105427825 << 0.000764221234825768 << -0.00178202876555804;

    double const abseps = std::pow(
      std::numeric_limits<double>::epsilon(), .4);
    auto res = restrictcdf::cdf<restrictcdf::deriv>(
      lower, upper, mean, sigma, false).approximate(1000000L, abseps, -1);

    expect_true(res.inform == 0L);
    expect_true(res.finest.n_elem == expect.n_elem);
    for(unsigned i = 0; i < expect.n_elem; ++i)
      expect_true(std::abs(res.finest[i] - expect[i]) < 1e-5);
  }
}
