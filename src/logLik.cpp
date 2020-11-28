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
static std::unique_ptr<double[]> log_ml_wk;
static size_t wk_mem_per_thread = 0L,
                current_wk_size = 0L;

double * get_working_memory_log_ml(){
#ifdef _OPENMP
  size_t const my_num = omp_get_thread_num();
#else
  size_t const my_num(0L);
#endif

  return log_ml_wk.get() + my_num * wk_mem_per_thread;
}

void set_working_memory_log_ml(size_t const n_int, size_t const n_obs,
                               size_t const n_threads){
  size_t const n_int_sq = n_int * n_int,
               n_obs_sq = n_obs * n_obs,
        size_shared_mem = n_int_sq + n_obs_sq + n_obs + n_obs * n_int,
          size_temp_mem = n_int_sq + n_obs_sq + n_int * n_obs + n_int,
                max_dim = size_shared_mem + size_temp_mem;

  constexpr size_t const cachline_size = 128L,
                                  mult = cachline_size / sizeof(double),
                              min_size = 2L * mult;

  size_t m_dim = max_dim;
  m_dim = std::max(m_dim, min_size);
  m_dim = (m_dim + mult - 1L) / mult;
  m_dim *= mult;
  wk_mem_per_thread = m_dim;

  size_t const new_size =
    std::max(n_threads, static_cast<size_t>(1L)) * m_dim;
  if(new_size > current_wk_size){
    log_ml_wk.reset(new double[new_size]);
    current_wk_size = new_size;

  }
}
} // namespace

namespace mdgc {
using namespace restrictcdf;

static double const log_2_pi = log(2 * M_PI);

double log_ml_term::approximate
(arma::mat const &vcov, arma::mat &derivs, int const maxpts,
 double const abs_eps, double const rel_eps, bool const comp_deriv,
 bool const do_reorder, size_t const minvls,
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
    if(comp_deriv and (vcov.n_cols != derivs.n_cols or vcov.n_rows != derivs.n_rows))
      throw std::invalid_argument("log_ml_term::approximate: invalid derivs");
  }
#endif
  double out(.0);

  // handle memory
  size_t const n_int_sq = n_int * n_int,
               n_obs_sq = n_obs * n_obs,
  /* memory that is needed for common quantities */
        size_shared_mem = n_int_sq + n_obs_sq + n_obs + n_obs * n_int;
  double * const wk_mem = get_working_memory_log_ml();

  size_t shared_cur = 0L;
  auto get_shared = [&](size_t const siz){
    double *out = wk_mem + shared_cur;
    shared_cur += siz;
    return out;
  };

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
  arma::mat S_oo(get_shared(n_obs_sq), n_obs, n_obs, false, true);
  arma::vec obs_scaled(get_shared(n_obs), n_obs, false, true);
  if(n_obs > 0){
    out -= static_cast<double>(n_obs) / 2. * log_2_pi;
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
      arma::mat S_00_inv(get_temp_mem(n_obs_sq, true), n_obs, n_obs,
                         false, true);
      if(!arma::inv_sympd(S_00_inv, S_oo))
        throw std::runtime_error("log_ml_term::approximate: inv() failed");
      for(size_t j = 0; j < n_obs; ++j){
        size_t const jj = idx_obs[j];
        for(size_t i = 0; i < n_obs; ++i){
          size_t const ii = idx_obs[i];
          // TODO: use that it is symmetric
          derivs.at(ii, jj) +=
            (obs_scaled[i] * obs_scaled[j] - S_00_inv.at(i , j)) / 2.;
        }
      }
    }
  }

  if(n_int > 0){
    arma::mat V(get_temp_mem(n_int_sq, true), n_int, n_int, false, true);
    V = vcov(idx_int, idx_int);
    arma::vec mea(get_temp_mem(n_int, false), n_int, false, true);
    mea.zeros();
    arma::mat S_oo_inv_S_oi(get_shared(n_obs * n_int), n_obs, n_int,
                            false, true);
    if(n_obs > 0){
      arma::mat S_oi(get_temp_mem(n_obs * n_int, false), n_obs, n_int,
                     false, true);
      S_oi = vcov(idx_obs, idx_int);
      if(!arma::solve(
        S_oo_inv_S_oi, S_oo, S_oi, arma::solve_opts::likely_sympd))
        throw std::runtime_error("log_ml_term::approximate: solve() failed");
      V -= S_oi.t() * S_oo_inv_S_oi;
      mea += S_oi.t() * obs_scaled;
    }

    if(comp_deriv){
      auto res =
        cdf<deriv>(lower, upper, mea, V, do_reorder,
                   use_aprx).approximate(
                       maxpts, abs_eps, rel_eps, minvls);
      double const p_hat = res.finest[0];
      arma::vec d_mu(res.finest.memptr() + 1L, n_int, false, true),
                 d_V(res.finest.memptr() + 1L + n_int,
                     (n_int * (n_int + 1L)) / 2L, false, true);
      d_mu /= p_hat;
      d_V  /= p_hat;

      arma::vec d_V_full(get_shared(n_int_sq), n_int_sq, false, true);
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
      if(n_obs > 0){
        {
          // TODO: memory allocation
          arma::mat inc = S_oo_inv_S_oi * d_mu * obs_scaled.t();
          inc /= 2.;
          derivs(idx_obs, idx_obs) -= inc;
          derivs(idx_obs, idx_obs) -= inc.t();
        }

        {
          // TODO: memory allocation
          arma::mat inc = d_mu * obs_scaled.t();

          inc /= 2.; // TODO: why?
          inc.reshape(n_int, n_obs);
          derivs(idx_int, idx_obs) += inc;
          derivs(idx_obs, idx_int) += inc.t();
        }
      }

      /* handle the terms from the covariance matrix */
      arma::mat dmat_V_full(d_V_full.memptr(), n_int, n_int, false, true);
      derivs(idx_int, idx_int) += dmat_V_full;

      if(n_obs > 0){
        {
          // TODO: memory allocation
          arma::mat const inc =
            S_oo_inv_S_oi * dmat_V_full * S_oo_inv_S_oi.t();
          derivs(idx_obs, idx_obs) += inc;
        }

        {
          // TODO: memory allocation
          arma::mat inc = 2 * dmat_V_full * S_oo_inv_S_oi.t();

          inc /= 2.; // TODO: why?
          derivs(idx_int, idx_obs) -= inc;
          derivs(idx_obs, idx_int) -= inc.t();
        }
      }

    } else {
      auto res =
        cdf<likelihood>(lower, upper, mea, V, do_reorder,
                        use_aprx).approximate(
                            maxpts, abs_eps, rel_eps);
      out += log(res.finest[0]);
    }
  }

  return out;
}

void log_ml_term::set_working_memory(
    std::vector<log_ml_term> const &terms, size_t const n_threads){
  size_t max_n_int = 0L, max_n_obs = 0L;
  for(auto const &t : terms){
    if(t.n_int > max_n_int)
      max_n_int = t.n_int;
    if(t.n_obs > max_n_obs)
      max_n_obs = t.n_obs;
  }

  cdf<deriv>::set_working_memory(max_n_int, n_threads);
  set_working_memory_log_ml(max_n_int, max_n_obs, n_threads);
}

} // namespace mdgc
