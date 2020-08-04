#include "logLik.h"
#include "restrict-cdf.h"
#include <algorithm>
#include "fast-commutation.h"
using std::log;

namespace mdgc {
using namespace restrictcdf;

static double const log_2_pi = log(2 * M_PI);

double log_ml_term::approximate
(arma::mat const &vcov, arma::mat &derivs, int const maxpts,
 double const abseps, double const releps, bool const comp_deriv,
 bool const do_reorder) const {
#ifdef DO_CHECKS
  {
    size_t const i1 =
      idx_obs.size() > 0 ?
      *std::max_element(idx_obs.begin(), idx_obs.end()) : 0L,
                 i2 =
      idx_int.size() > 0 ?
      *std::max_element(idx_int.begin(), idx_int.end()) : 0L;
    size_t const max_idx = std::max(i1, i2);
    if(vcov.n_rows < max_idx + 1L or vcov.n_cols < max_idx + 1L)
      throw std::invalid_argument("log_ml_term::approximate: invalid vcov");
    if(comp_deriv and (vcov.n_cols != derivs.n_cols or vcov.n_rows != derivs.n_rows))
      throw std::invalid_argument("log_ml_term::approximate: invalid derivs");
  }
#endif

  // TODO: this function be written a lot smarter...
  double out(.0);

  /* add terms from the observed outcomes */
  arma::mat S_oo;
  arma::vec obs_scaled;
  if(n_obs > 0){
    out -= static_cast<double>(n_obs) / 2. * log_2_pi;
    S_oo = vcov(idx_obs, idx_obs);
    double val, sgn;
    arma::log_det(val, sgn, S_oo);
    out -= .5 * val;
    obs_scaled = arma::solve(
      S_oo, obs_val, arma::solve_opts::likely_sympd);
    out -= arma::dot(obs_val, obs_scaled) / 2.;

    if(comp_deriv){
      derivs(idx_obs, idx_obs) += obs_scaled * obs_scaled.t() / 2.;
      // TODO: the inverse term is computed many times!
      derivs(idx_obs, idx_obs) -= S_oo.i() / 2.;
    }
  }

  if(n_int > 0){
    arma::mat V = vcov(idx_int, idx_int);
    arma::vec mea = -upper;
    arma::mat S_oi, S_oo_inv_S_oi;
    if(n_obs > 0){
      S_oi = vcov(idx_obs, idx_int);
      S_oo_inv_S_oi =  arma::solve(
        S_oo, S_oi, arma::solve_opts::likely_sympd);
      V -= S_oi.t() * S_oo_inv_S_oi;
      mea += S_oi.t() * obs_scaled;
    }

    if(comp_deriv){
      auto res =
        cdf<deriv>(mea, V, do_reorder).approximate(maxpts, abseps, releps);
      double const p_hat = res.finest[0];
      arma::vec d_mu(res.finest.memptr() + 1L, n_int, false),
                 d_V(res.finest.memptr() + 1L + n_int,
                     (n_int * (n_int + 1L)) / 2L, false);
      d_mu /= p_hat;
      d_V  /= p_hat;

      arma::vec d_V_full(n_int * n_int);
      {
        double const *val = d_V.begin();
        for(size_t c = 0; c < n_int; ++c){
          size_t const inc = c * n_int;
          for(size_t r = 0; r  <= c; ++r, ++val){
            d_V_full[r         + inc] = *val;
            d_V_full[r * n_int + c  ] = *val;
          }
        }
      }

      out += log(p_hat);

      /* handle the terms from the mean */
      arma::mat const dum_diag_mat =
        arma::diagmat(arma::vec(S_oi.n_cols, arma::fill::ones));

      if(n_obs > 0){
        arma::mat tmp = arma::kron(obs_scaled, S_oo_inv_S_oi);
        arma::mat inc = tmp * d_mu;
        inc.reshape(n_obs, n_obs);
        derivs(idx_obs, idx_obs) -= (inc + inc.t()) / 2.;

        tmp = arma::kron(obs_scaled, dum_diag_mat);
        inc = tmp * d_mu;
        inc /= 2.;
        inc.reshape(n_int, n_obs);
        derivs(idx_int, idx_obs) += inc;
        derivs(idx_obs, idx_int) += inc.t();

      }

      /* handle the terms from the covariance matrix */
      {
        arma::mat dum(d_V_full.memptr(), n_int, n_int, false);
        derivs(idx_int, idx_int) += dum;
      }

      if(n_obs > 0){
        arma::mat tmp = arma::kron(S_oo_inv_S_oi, S_oo_inv_S_oi);
        arma::mat inc = tmp * d_V_full;
        inc.reshape(n_obs, n_obs);
        derivs(idx_obs, idx_obs) += inc;

        arma::mat const K = get_commutation(
          S_oo_inv_S_oi.n_cols, S_oo_inv_S_oi.n_rows);
        tmp = arma::kron(S_oo_inv_S_oi, dum_diag_mat) +
          K.t() * arma::kron(dum_diag_mat, S_oo_inv_S_oi);
        inc = tmp * d_V_full;

        inc.reshape(n_int, n_obs);

        inc /= 2.;
        derivs(idx_int, idx_obs) -= inc;
        derivs(idx_obs, idx_int) -= inc.t();

      }

    } else {
      auto res =
        cdf<likelihood>(mea, V, do_reorder).approximate(
            maxpts, abseps, releps);
      out += log(res.finest[0]);
    }
  }

  return out;
}

void log_ml_term::set_working_memory(
    std::vector<log_ml_term> const &terms, size_t const n_threads){
  size_t max_n_int = 0L;
  for(auto const &t : terms)
    if(t.n_int > max_n_int)
      max_n_int = t.n_int;

  cdf<deriv>::set_working_memory(max_n_int, n_threads);
}

} // namespace mdgc
