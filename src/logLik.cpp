#include "logLik.h"
#include "restrict-cdf.h"
using std::log;

namespace mdgc {
using namespace restrictcdf;

static double const log_2_pi = log(2 * M_PI);

double log_ml_term::approximate
(arma::mat const &vcov, arma::mat &derivs, int const maxpts,
 double const abseps, double const releps) const {
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

  }

  if(n_int > 0){
    arma::mat V = vcov(idx_int, idx_int);
    arma::vec mea = -upper;
    if(n_obs > 0){
      arma::mat S_oi = vcov(idx_obs, idx_int);
      V -= S_oi.t() * arma::solve(
        S_oo, S_oi, arma::solve_opts::likely_sympd);
      mea += S_oi.t() * obs_scaled;
    }

    auto res = cdf<likelihood>(mea, V).approximate(maxpts, abseps, releps);
    out += log(res.finest[0]);
  }

  return out;
}

} // namespace mdgc
