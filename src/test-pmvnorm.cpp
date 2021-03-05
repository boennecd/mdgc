#include "mvtnorm-wrapper.h"
#include <testthat.h>
#include <limits>
#include "threat-safe-random.h"

context("pmvnorm unit tests") {
  test_that("pmvnorm gives similar output to R") {
/*
 set.seed(1)
 n <- 4
 lower <- rnorm(n)
 upper <- lower + abs(rnorm(n))
 mean <- rnorm(n)
 sigma <- drop(rWishart(1L, 2L * n, diag(n)))

 lower[1:3] <- c(-Inf    , lower[2], -Inf)
 upper[1:3] <- c(upper[2], Inf     ,  Inf)
 lower <- round(lower, 3)
 upper <- round(upper, 3)
 mean  <- round(mean , 3)
 sigma <- round(sigma, 3)

 library(mvtnorm)
 prob <- pmvnorm(lower, upper, mean, sigma = sigma,
 algorithm = GenzBretz(abseps = 1e-9, maxpts = 1000000L))

 dput(lower)
 dput(upper)
 dput(mean)
 dput(sigma)
 dput(prob)
*/
    Rcpp::RNGScope rngScope;
    parallelrng::set_rng_seeds(1L);
    constexpr double Inf = std::numeric_limits<double>::infinity();

    arma::vec lower, upper, mean;
    arma::mat sigma;

    lower << -Inf << 0.184 << -Inf << 1.595;
    upper << 1.004 << Inf << Inf << 2.334;
    mean << 0.576 << -0.305 << 1.512 << 0.39;

    sigma << 4.869 << -0.099 << 0.961 << 1.726 << -0.099 << 5.01 << -2.789
          << 0.132 << 0.961 << -2.789 << 6.67 << -4.177 << 1.726 << 0.132
          << -4.177 << 7.966;
    sigma.reshape(4L, 4L);

    double const abs_eps = std::sqrt(std::numeric_limits<double>::epsilon());
    {
      auto res = pmvnorm::cdf(lower, upper, mean, sigma, 1000000L, abs_eps,
                              -1);
      expect_true(res.inform == 0L);
      expect_true(res.error                                 < 100. * abs_eps);
      expect_true(std::fabs(res.value - 0.0196461341023563) < 100. * abs_eps);
    }

    auto const infin = pmvnorm::get_infin(lower, upper);
    auto const cor_vec_res = pmvnorm::get_cor_vec(sigma);

    lower /= cor_vec_res.sds;
    upper /= cor_vec_res.sds;
    mean  /= cor_vec_res.sds;

    {
      auto res = pmvnorm::cdf(lower, upper, infin, mean,
                              cor_vec_res.cor_vec, 1000000L, abs_eps, -1);
      expect_true(res.inform == 0L);
      expect_true(res.error                                 < 100. * abs_eps);
      expect_true(std::fabs(res.value - 0.0196461341023563) < 100. * abs_eps);
    }
  }
}
