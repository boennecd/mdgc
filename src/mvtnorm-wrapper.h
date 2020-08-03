#ifndef MVTNORM_WRAPPER_H
#define MVTNORM_WRAPPER_H
#include "arma-wrap.h"

namespace pmvnorm {
/**
 * @param lower The lower bounds.
 * @param upper The upper bounds.
 * @return the infin argument for the mvtdst subroutine.
 */
arma::ivec get_infin(arma::vec const &lower, arma::vec const &upper);

struct cor_vec_res {
  arma::vec cor_vec, sds;
};

/**
 * @return a struct with the correlation matrix and standard deviation. The
 * correlation matrix is stored as a upper diagonal matrix.
 */
cor_vec_res get_cor_vec(const arma::mat&);

struct cdf_res {
  double error, value;
  int inform, intvls;
};

/**
 * Approximates the multivariate normal CDF over a hyperrectangle.
 *
 * @param lower Lower bounds.
 * @param upper Upper bounds.
 * @param mean Mean vector.
 * @param cov Covariance matrix.
 * @param maxpts Maximum number of integrand evaluations.
 * @param abseps Absolute convergence threshold.
 * @param releps Relative convergence threshold.
 */
cdf_res cdf(arma::vec lower, arma::vec upper, arma::vec mean,
            arma::mat const &cov, int const maxpts = -1L,
            double const abseps = -1, double const releps = 1e-5);

/**
 * Approximates the multivariate normal CDF over a hyperrectangle.
 *
 * @param lower Lower bounds.
 * @param upper Upper bounds.
 * @param infin infin argument for the mvtdst subroutine.
 * @param mean Mean vector.
 * @param cor_vec upper triangle of the correlation matrix.
 * @param maxpts Maximum number of integrand evaluations.
 * @param abseps Absolute convergence threshold.
 * @param releps Relative convergence threshold.
 */
cdf_res cdf(arma::vec const &lower, arma::vec const &upper,
            arma::ivec const &infin, arma::vec const &mean,
            arma::vec const &cor_vec, int const maxpts = -1L,
            double const abseps = -1, double const releps = 1e-5);
}

#endif
