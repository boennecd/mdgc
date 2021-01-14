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
  /** indices for interval observations and observed values. */
  arma::uvec const idx_int, idx_obs;

  /// number of censored and observed values.
  inline arma::uword n_int() const noexcept {
     return idx_int.n_elem;
  }
  inline arma::uword n_obs() const noexcept {
     return idx_obs.n_elem;
  }

  /// observed values
  arma::vec const obs_val;
  /*
   * 3xn matrix with multinomial outcomes. The first index is the outcome,
   * the second index is the number of categories, and the third index is
   * the index of the first latent variable.
   *
   * TODO: we only need one of these as the other two are likely the same
   * for all log_ml_terms. Thus, we can reduce the storage requirements.
   */
  arma::imat const multinomial;
  /// true if there are any multinomial variables
  inline bool any_mult() const noexcept {
     return multinomial.n_cols > 0;
  }
  /// number of multinomial outcomes
  inline arma::uword n_cate() const noexcept {
     return multinomial.n_cols;
  }
  /** indices for the latent variables of the observed multinomial levels.
   *  The vector is empty if any_mult() is false. */
  arma::uvec const idx_cat_obs;
  /** indices for the latent variables that needs to be integrated out which
   *  are not in idx_cat_obs. The vector is empty if any_mult() is false. */
  arma::uvec const idx_not_cat_obs;

  /** stores lower bounds for the CDF. Notice that multinomial add the
   *  number of levels less one with negative infinity. */
  arma::vec const lower;
  /** stores upper bounds for the CDF. Notice that multinomial add the
   *  number of levels less one with zeros. */
  arma::vec const upper;

public:
  log_ml_term(arma::uvec const &idx_int, arma::uvec const &idx_obs,
              arma::vec const &obs_val, arma::vec const &lower_in,
              arma::vec const &upper_in, arma::imat const &multinomial):
  idx_int(idx_int), idx_obs(idx_obs), obs_val(obs_val),
  multinomial(multinomial),

  idx_cat_obs(([&]() -> arma::uvec {
    if(!any_mult())
      return arma::uvec();

    arma::uvec out(n_cate());
    size_t k = 0;
    for(int j = 0; j < static_cast<int>(n_int()) and k < n_cate(); ++j)
      if(idx_int[j] == static_cast<size_t>(multinomial.at(2, k))){
        int const idx_obs_lvl = multinomial.at(0, k) - 1 + j,
                   n_lvls     = multinomial.at(1, k);

        for(int l = 0; l < n_lvls; ++l, ++j)
          if(j == idx_obs_lvl)
            // found the one we want
            break;

        out[k] = j;
        ++k;
      }

    return out;
  })()),

  idx_not_cat_obs(([&]() -> arma::uvec {
    if(!any_mult())
      return arma::uvec();

    arma::uvec out(n_int() - n_cate());
    size_t k = 0;
    for(int j = 0; j < static_cast<int>(n_int()); ++j){
      if(k < n_cate() and
           idx_int[j] == static_cast<size_t>(multinomial.at(2, k))){
        int const idx_obs_lvl = multinomial.at(0, k) - 1 + j,
                   n_lvls     = multinomial.at(1, k);

        for(int l = 0; l < n_lvls; ++l, ++j)
          if(j == idx_obs_lvl)
            // found the one we do not want
            break;
          else
            out[j - k] = j;

        ++k;
        continue;
      }

      out[j - k] = j;
    }

    return out;
  })()),

  lower(([&]() -> arma::vec {
    if(lower_in.n_elem != n_int())
      throw std::invalid_argument("log_ml_term::log_ml_term: invalid lower_in");

    if(!any_mult())
      return lower_in;

    arma::vec out(n_int() - n_cate());
    size_t k(0);
    for(int j = 0; j < static_cast<int>(n_int()); ++j){
      if(static_cast<size_t>(k) < multinomial.n_cols and
           idx_int[j] == static_cast<size_t>(multinomial.at(2, k))){
          int const idx_obs_lvl = multinomial.at(0, k) - 1 + j,
                     n_lvls     = multinomial.at(1, k);

         for(int l = 0; l < n_lvls; ++l, ++j)
            if(j != idx_obs_lvl)
               out[j - k - (l >= multinomial.at(0, k))] = lower_in[j];

         --j; // there is an increment in the for-loop
         ++k;
         continue;
      }

      out[j - k] = lower_in[j];
    }

    return out;
  })()),

  upper(([&]() -> arma::vec {
    if(upper_in.n_elem != n_int())
      throw std::invalid_argument("log_ml_term::log_ml_term: invalid upper_in");

    if(!any_mult())
      return upper_in;

    arma::vec out(n_int() - n_cate());
    size_t k(0);
    for(int j = 0; j < static_cast<int>(n_int()); ++j){
      if(static_cast<size_t>(k) < multinomial.n_cols and
           idx_int[j] == static_cast<size_t>(multinomial.at(2, k))){
        int const idx_obs_lvl = multinomial.at(0, k) - 1 + j,
                  n_lvls      = multinomial.at(1, k);

        for(int l = 0; l < n_lvls; ++l, ++j)
          if(j != idx_obs_lvl)
            out[j - k - (l >= multinomial.at(0, k))] =
              upper_in[j] - upper_in[idx_obs_lvl];

        --j; // there is an increment in the for-loop
        ++k;
        continue;
      }

      out[j - k] = upper_in[j];
    }

    return out;
  })()) {
    if(obs_val.n_elem != n_obs())
      throw std::invalid_argument("log_ml_term::log_ml_term: invalid 'obs_val'");
    if(lower.n_elem != n_int() - n_cate())
      throw std::invalid_argument("log_ml_term::log_ml_term: invalid 'lower'");
    if(upper.n_elem != n_int() - n_cate())
      throw std::invalid_argument("log_ml_term::log_ml_term: invalid 'upper'");
    if(any_mult() and multinomial.n_rows != 3)
       throw std::invalid_argument("log_ml_term::log_ml_term: invalid 'multinomial'");
    if(idx_cat_obs.n_elem != n_cate())
      throw std::runtime_error("log_ml_term::log_ml_term: created invalid idx_cat_obs");
  }

  /**
   * Approximates the log marginal likelihood and its gradient.
   *
   * @param vcov Covariance matrix at which to evalute the log marginal
   * likelihood term.
   * @param mu Mean vector at which to evaluate the log marginal likelihood
   * term. Does not have any effect for continuous variables.
   * @param derivs_vcov Matrix with derivatives to increment for the
   * covariance matrix.
   * @param derivs_mea Vector with derivatives to increament for the
   * mean vector.
   * @param maxpts Maximum number of integrand evaluations.
   * @param abs_eps Absolute convergence threshold.
   * @param rel_eps Relative convergence threshold.
   * @param comp_deriv true if the derivative should be computed.
   * @param do_reorder true if the order of integrations may be reordered.
   * @param minvls Minimum number of samples.
   * @param use_aprx logical for whether to use an approximation of
   * the normal CDF.
   *
   * Mean terms for observed variables are always assumed to be zero. The
   * corresponding derivatives are set to zero.
   *
   * @return the log marginal likelihood approximation.
   *
   */
  double approximate(arma::mat const &vcov, arma::vec const &mu,
                     arma::mat &derivs_vcov, arma::vec &derivs_mea,
                     int const maxpts, double const abs_eps,
                     double const rel_eps, bool const comp_deriv,
                     bool const do_reorder, size_t const minvls,
                     bool const use_aprx) const;

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
