#include "logLik.h"
#include "restrict-cdf.h"
#include <algorithm>
#include "fast-commutation.h"
#include "lp_utils.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using std::log;

namespace {
static cache_mem<double> log_ml_mem;

inline double * get_working_memory_log_ml() MDGC_NOEXCEPT {
  return log_ml_mem.get_mem();
}

void set_working_memory_log_ml(size_t const n_int, size_t const n_obs,
                               size_t const n_threads){
  size_t const n_int_sq = n_int * n_int,
               n_obs_sq = n_obs * n_obs,
        size_shared_mem = n_int_sq + n_obs_sq + n_obs + n_obs * n_int,
          size_temp_mem = n_int_sq + n_obs_sq + n_int * n_obs + 2 * n_int,
                max_dim = size_shared_mem + size_temp_mem;

  log_ml_mem.set_n_mem(max_dim, n_threads);
}
} // namespace

namespace mdgc {
using namespace restrictcdf;

static double const log_2_pi = log(2 * M_PI);

double log_ml_term::approximate
(arma::mat const &vcov, arma::vec const &mu, arma::mat &derivs_vcov,
 arma::vec &derivs_mea,
 int const maxpts, double const abs_eps, double const rel_eps,
 bool const comp_deriv,bool const do_reorder, size_t const minvls,
 bool const use_aprx) const {
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
    if(comp_deriv and (vcov.n_cols != derivs_vcov.n_cols or vcov.n_rows != derivs_vcov.n_rows))
      throw std::invalid_argument("log_ml_term::approximate: invalid derivs_vcov");
    if(comp_deriv and derivs_mea.n_elem != vcov.n_cols)
      throw std::invalid_argument("log_ml_term::approximate: invalid derivs_mea");
  }
#endif
  double out(.0);

  // handle memory
  size_t const n_int_sq = n_int() * n_int(),
               n_obs_sq = n_obs() * n_obs(),
  /* memory that is needed for common quantities */
        size_shared_mem = n_int_sq + n_obs_sq + n_obs() + n_obs() * n_int();
  double * const wk_mem = get_working_memory_log_ml();

  size_t shared_cur = 0L;
  auto get_shared = [&](size_t const siz){
    double *out = wk_mem + shared_cur;
    shared_cur += siz;
    return out;
  };

  // lambda function to return working memory
  size_t tmp_cur = 0L;
  double * const wk_use = wk_mem + size_shared_mem;
  auto get_temp_mem = [&](size_t const siz, bool const reset){
    if(reset)
      tmp_cur = 0L;
    double *out = wk_use + tmp_cur;
    tmp_cur += siz;
    return out;
  };

  /* add terms from the observed outcomes */
  arma::mat S_oo(get_shared(n_obs_sq), n_obs(), n_obs(), false, true);
  arma::vec obs_scaled(get_shared(n_obs()), n_obs(), false, true);
  if(n_obs() > 0){
    out -= static_cast<double>(n_obs()) / 2. * log_2_pi;
    S_oo = vcov(idx_obs, idx_obs);
    double val, sgn;
    arma::log_det(val, sgn, S_oo);
    out -= .5 * val;
    if(!arma::solve(
      obs_scaled, S_oo, obs_val, arma::solve_opts::likely_sympd))
      throw std::runtime_error("log_ml_term::approximate: solve() failed");
    out -= arma::dot(obs_val, obs_scaled) / 2.;

    if(comp_deriv){
      // TODO: the inverse term is computed many times!
      arma::mat S_00_inv(get_temp_mem(n_obs_sq, true), n_obs(), n_obs(),
                         false, true);
      if(!arma::inv_sympd(S_00_inv, S_oo))
        throw std::runtime_error("log_ml_term::approximate: inv() failed");
      for(size_t j = 0; j < n_obs(); ++j){
        size_t const jj = idx_obs[j];
        for(size_t i = 0; i < n_obs(); ++i){
          size_t const ii = idx_obs[i];
          // TODO: use that it is symmetric
          derivs_vcov.at(ii, jj) +=
            (obs_scaled[i] * obs_scaled[j] - S_00_inv.at(i , j)) / 2.;
        }
      }
    }
  }

  if(n_int() > 0){
    arma::mat V(get_temp_mem(n_int_sq, true), n_int(), n_int(), false);
    V = vcov(idx_int, idx_int);

    arma::vec mea(get_temp_mem(n_int(), false), n_int(), false);
    mea = mu(idx_int);
    arma::mat S_oo_inv_S_oi(get_shared(n_obs() * n_int()), n_obs(), n_int(),
                            false, true);

    if(n_obs() > 0){
      arma::mat S_oi(get_temp_mem(n_obs() * n_int(), false), n_obs(),
                     n_int(), false, true);
      S_oi = vcov(idx_obs, idx_int);
      if(!arma::solve(
        S_oo_inv_S_oi, S_oo, S_oi, arma::solve_opts::likely_sympd))
        throw std::runtime_error("log_ml_term::approximate: solve() failed");
      V -= S_oi.t() * S_oo_inv_S_oi;
      mea += S_oi.t() * obs_scaled;
    }

    // handle that we only need a subset of the columns
    // TODO: memory allocation and can be avoided with simple loops
    arma::mat D(n_int() - n_cate(), n_cate(), arma::fill::zeros);
    if(any_mult()){
      size_t j(0);
      for(size_t k = 0; k < D.n_cols; ++k)
        for(; j < n_int(); ++j)
          if(idx_int[j] == static_cast<size_t>(multinomial.at(2, k))){
            int const n_lvls = multinomial.at(1, k);
            for(int l = 0; l < n_lvls - 1; ++l, ++j)
              D.at(j - k, k) = 1;
            ++j;
            break;
          }

      // TODO: memory allocation
      mea = mea(idx_not_cat_obs) - D * mea(idx_cat_obs);
      // TODO: memory allocation
      arma::mat const D_sig = D * V(idx_cat_obs, idx_not_cat_obs);
      // TODO: memory allocation
      V = V(idx_not_cat_obs, idx_not_cat_obs) +
        D * V(idx_cat_obs, idx_cat_obs) * D.t() - D_sig - D_sig.t();
    }

    if(comp_deriv){
      deriv functor(mea, V);
      auto res =
        cdf<deriv>(functor, lower, upper, mea, V, do_reorder,
                   use_aprx).approximate(
                       maxpts, abs_eps, rel_eps, minvls);
      double const p_hat = res.likelihood;
      out += log(p_hat);

      arma::uword const dim_int = n_int() - n_cate();
      arma::vec d_mu(res.derivs.begin(), dim_int, false, true),
                 d_V(d_mu      .end()  ,
                     (dim_int * (dim_int + 1L)) / 2L, false, true);
      d_mu /= p_hat;
      d_V  /= p_hat;

      arma::vec d_mu_full(get_shared(n_int()), n_int(), false, true);
      arma::mat d_V_full (get_shared(n_int_sq), n_int(), n_int(), false,
                          true);
      if(!any_mult()){
        double const *val = d_V.begin();
        for(size_t c = 0; c < n_int(); ++c){
          for(size_t r = 0; r  <= c; ++r, ++val){
            d_V_full.at(r, c) = *val;
            d_V_full.at(c, r) = *val;
          }
        }

        d_mu_full = d_mu;

      } else {
        // TODO: memory allocation
        arma::mat dum(dim_int, dim_int);
        {
          double const *val = d_V.begin();
          for(size_t c = 0; c < dim_int; ++c){
            for(size_t r = 0; r  <= c; ++r, ++val){
              dum.at(r, c) = *val;
              dum.at(c, r) = *val;
            }
          }
        }

        // TODO: memory allocation
        arma::mat const dum_D = dum * D;
        d_V_full(idx_not_cat_obs, idx_not_cat_obs) = dum;
        // TODO: memory allocation
        d_V_full(idx_cat_obs    , idx_cat_obs    ) = D.t() * dum_D;
        d_V_full(idx_cat_obs    , idx_not_cat_obs) = -dum_D.t();
        d_V_full(idx_not_cat_obs, idx_cat_obs    ) = -dum_D;

        d_mu_full(idx_not_cat_obs) = d_mu;
        // TODO: memory allocation
        d_mu_full(idx_cat_obs    ) = -D.t() * d_mu;

      }

      /* handle the terms from the mean */
      derivs_mea(idx_int) += d_mu_full;

      if(n_obs() > 0){
        {
          // TODO: memory allocation
          arma::mat inc = S_oo_inv_S_oi * d_mu_full * obs_scaled.t();
          inc /= 2.;
          derivs_vcov(idx_obs, idx_obs) -= inc;
          derivs_vcov(idx_obs, idx_obs) -= inc.t();
        }

        {
          // TODO: memory allocation
          arma::mat inc = d_mu_full * obs_scaled.t();

          inc /= 2.;
          inc.reshape(n_int(), n_obs());
          derivs_vcov(idx_int, idx_obs) += inc;
          derivs_vcov(idx_obs, idx_int) += inc.t();
        }
      }

      /* handle the terms from the covariance matrix */
      derivs_vcov(idx_int, idx_int) += d_V_full;

      if(n_obs() > 0){
        {
          // TODO: memory allocation
          arma::mat const inc =
            S_oo_inv_S_oi * d_V_full * S_oo_inv_S_oi.t();
          derivs_vcov(idx_obs, idx_obs) += inc;
        }

        {
          // TODO: memory allocation
          arma::mat inc = d_V_full * S_oo_inv_S_oi.t();

          derivs_vcov(idx_int, idx_obs) -= inc;
          derivs_vcov(idx_obs, idx_int) -= inc.t();
        }
      }

    } else {
      likelihood functor;
      auto res =
        cdf<likelihood>(functor, lower, upper, mea, V, do_reorder,
                        use_aprx).approximate(maxpts, abs_eps, rel_eps);
      out += log(res.likelihood);
    }
  }

  return out;
}

void log_ml_term::set_working_memory(
    std::vector<log_ml_term> const &terms, size_t const n_threads){
  size_t max_n_int = 0L, max_n_obs = 0L;
  for(auto const &t : terms){
    if(t.n_int() > max_n_int)
      max_n_int = t.n_int();
    if(t.n_obs() > max_n_obs)
      max_n_obs = t.n_obs();
  }

  likelihood::alloc_mem(max_n_int, n_threads);
  deriv     ::alloc_mem(max_n_int, n_threads);
  set_working_memory_log_ml(max_n_int, max_n_obs, n_threads);
}

} // namespace mdgc
