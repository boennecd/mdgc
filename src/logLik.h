#ifndef LOGLIK_H
#define LOGLIK_H
#include "arma-wrap.h"
#include "mvtnorm-wrapper.h"

namespace mdgc {
/***
 * class to approximate the marginal likelihood term for a given
 * observation.
 */
class log_ml_term {
  /// indices for interval observation and observed values
  arma::uvec const idx_int,
                   idx_obs;

  /// number of censored
  size_t const n_int = idx_int.n_elem,
               n_obs = idx_obs.n_elem;

  /// observed values
  arma::vec const obs_val;
  /// stores lower bounds for the CDF
  arma::vec const lower;
  /// stores upper bounds for the CDF
  arma::vec const upper;
  /// infin values to pass to CDF approximation
  arma::ivec const infin = pmvnorm::get_infin(lower, upper);

public:
  log_ml_term(arma::uvec const &idx_int, arma::uvec const &idx_obs,
              arma::vec const &obs_val, arma::vec const &lower,
              arma::vec const &upper):
  idx_int(idx_int), idx_obs(idx_obs), obs_val(obs_val), lower(lower),
  upper(upper) {
#ifdef DO_CHECKS
    if(obs_val.n_elem != n_obs)
      throw std::invalid_argument("log_ml_term::log_ml_term: invalid 'obs_val'");
    if(lower.n_elem != n_int)
      throw std::invalid_argument("log_ml_term::log_ml_term: invalid 'lower'");
    if(upper.n_elem != n_int)
      throw std::invalid_argument("log_ml_term::log_ml_term: invalid 'upper'");
#endif
  }

  /**
   * Approximates the log marginal likelihood and its gradient.
   *
   * @param vcov Covariance matrix at which to evalute the log marginal
   * likelihood term.
   * @param derivs Matrix with derivatives to increament.
   * @param maxpts Maximum number of integrand evaluations.
   * @param abseps Absolute convergence threshold.
   * @param releps Relative convergence threshold.
   * @param comp_deriv true if the derivative should be computed.
   * @param do_reorder true if the order of integrations may be reordered.
   *
   * @return the log marginal likelihood approximation.
   *
   */
  double approximate(arma::mat const &vcov, arma::mat &derivs,
                     int const maxpts, double const abseps,
                     double const releps, bool const comp_deriv,
                     bool const do_reorder) const;

  /**
   * sets the working memory. Must be called prior to calling approximate.
   *
   * @param terms Vector with terms to be approximated.
   * @param n_threads Number of threads that will be used in the
   * approximation.
   */
  static void set_working_memory(
      std::vector<log_ml_term> const &terms, size_t const n_threads);
};

} // namespace mdgc

#endif
