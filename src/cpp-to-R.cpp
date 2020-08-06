#include "logLik.h"
#include "threat-safe-random.h"
#include "pnorm.h"
#include "qnorm.h"
#include "fast-commutation.h"
#include "restrict-cdf.h"
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "lp_utils.h"

using namespace mdgc;
using namespace arma;

struct ml_terms {
  size_t n_variables;
  std::vector<log_ml_term> terms;
};


/**
 * get objects to evaluate the log marginal likelihood terms.
 *
 * @param lower Matrix with lower bounds.
 * @param upper Matrix with upper bounds. If code is zero for the entry
 * then it is the observed value.
 * @param code Matrix with type of outcomes. 0 is an observed value,
 * 1 is a missing value, and 2 is a value in an interval.
 */
// [[Rcpp::export(rng = false)]]
SEXP get_log_lm_terms(arma::mat const &lower, arma::mat const &upper,
                      arma::imat const &code){
  auto out = Rcpp::XPtr<ml_terms>(new ml_terms());

  size_t const n = lower.n_cols,
               p = lower.n_rows;
  if(upper.n_rows != p or upper.n_cols != n)
    throw std::invalid_argument("get_log_lm_terms: invalid 'upper'");
  if(code.n_rows != p or code.n_cols != n)
    throw std::invalid_argument("get_log_lm_terms: invalid 'code'");

  /* fill in the log ml terms objects */
  std::vector<log_ml_term> &terms = out->terms;
  out->n_variables = p;
  terms.reserve(n);
  uvec w_idx_int(p),
       w_idx_obs(p);
  vec w_obs_val(p), w_upper(p), w_lower(p);

  for(size_t i = 0; i < n; ++i){
    size_t n_o(0L), n_i(0L);

    for(size_t j = 0; j < p; ++j){
      if       (code.at(j, i) == 0L){
        /* observed value */
        w_idx_obs.at(n_o  ) = j;
        w_obs_val.at(n_o++) = upper.at(j, i);
      } else if(code.at(j, i) == 1L) {
        /* do nothing. The value is missing. */
      } else if(code.at(j, i) == 2L){
        /* Z is in an interval */
        w_idx_int.at(n_i  ) = j;
        w_lower  .at(n_i  ) = lower.at(j, i);
        w_upper  .at(n_i++) = upper.at(j, i);

      } else
        throw std::invalid_argument("get_log_lm_terms: invalid code");
    }

    if(n_o == 0L and n_i == 0L)
      throw std::invalid_argument("get_log_lm_terms: no observed value");

    uvec a_idx_int, a_idx_obs;
    arma::vec a_obs_val, a_upper, a_lower;
    if(n_i > 0){
      a_idx_int = w_idx_int.subvec(0L, n_i - 1L);
      a_upper   = w_upper  .subvec(0L, n_i - 1L);
      a_lower   = w_lower  .subvec(0L, n_i - 1L);

    } else {
      a_idx_int.set_size(0L);
      a_upper  .set_size(0L);
      a_lower  .set_size(0L);

    }

    if(n_o > 0){
      a_idx_obs = w_idx_obs.subvec(0L, n_o - 1L);
      a_obs_val = w_obs_val.subvec(0L, n_o - 1L);

    } else {
      a_idx_obs.set_size(0L);
      a_obs_val.set_size(0L);

    }

    terms.emplace_back(a_idx_int, a_idx_obs, a_obs_val, a_lower, a_upper);
  }

  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector eval_log_lm_terms(
    SEXP ptr, arma::ivec const &indices, arma::mat const &vcov,
    int const maxpts, double const abseps, double const releps,
    size_t const n_threads, bool const comp_derivs,
    bool const do_reorder = true){
  Rcpp::XPtr<ml_terms> obj(ptr);
  std::vector<log_ml_term> const &terms = obj->terms;

  size_t const p = obj->n_variables;
#ifdef DO_CHECKS
  if(vcov.n_cols != p or vcov.n_rows != p)
    throw std::invalid_argument("eval_log_lm_terms: invalid vcov");
  if(n_threads < 1L)
    throw std::invalid_argument("eval_log_lm_terms: invalid n_threads");
  if(indices.n_elem > terms.size())
    throw std::invalid_argument("eval_log_lm_terms: invalid indices");
#endif
  log_ml_term::set_working_memory(terms, n_threads);

  arma::mat derivs = comp_derivs ? mat(p, p, arma::fill::zeros) : mat();
#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif
  parallelrng::set_rng_seeds(n_threads);

  double out(0.);
  size_t const n_indices = indices.n_elem;
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    arma::mat my_derivs(comp_derivs ? p : 0L, comp_derivs ? p : 0L,
                        arma::fill::zeros);

#ifdef _OPENMP
#pragma omp for schedule(static) reduction(+:out)
#endif
    for(size_t i = 0; i < n_indices; ++i)
      out += terms[indices[i]].approximate(
        vcov, my_derivs, maxpts, abseps, releps, comp_derivs,
        do_reorder);

    if(comp_derivs)
#ifdef _OPENMP
#pragma omp critical(add_derivs)
#endif
      derivs += my_derivs;
  }

  Rcpp::NumericVector res(1);
  res[0] = out;
  if(comp_derivs)
    res.attr("grad") = derivs;

  return res;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix get_z_hat
(arma::mat const &lower, arma::mat const &upper, arma::imat const &code,
 unsigned const n_threads){
  size_t const p = lower.n_rows,
               n = upper.n_cols;
  if(upper.n_rows != p or upper.n_cols != n)
    throw std::invalid_argument("get_z_hat: invalid upper");
  if(code.n_rows != p or code.n_cols != n)
    throw std::invalid_argument("get_z_hat: invalid lower");
  if(n_threads < 1)
    throw std::invalid_argument("get_z_hat: invalid n_threads");

  Rcpp::NumericMatrix out(p, n);
  double * const o = &out[0];
#ifdef _OPENMP
#pragma omp for
#endif
  for(size_t j = 0; j < n; ++j){
    double * oj = o + j * p;
    for(size_t i = 0; i < p; ++i, ++oj){
      if(code.at(i, j) <= 1L){
        *oj = upper.at(i, j);
        continue;
      }

      double const l = lower.at(i, j),
                   u = upper.at(i, j),
                   a = std::isinf(l) ? 0 : pnorm_std(l, 1L, 0L),
                   b = std::isinf(u) ? 1 : pnorm_std(u, 1L, 0L);

      *oj = qnorm_w((b + a) / 2., 0., 1., 1L, 0L);
    }
  }

  return out;
}

// [[Rcpp::export("pmvnorm")]]
Rcpp::NumericVector pmvnorm_to_R
  (arma::vec const &lower, arma::vec const &upper, arma::vec const &mu,
   arma::mat const &Sigma, int const maxvls, double const abseps,
   double const releps, bool const derivs, bool const do_reorder = true){
  parallelrng::set_rng_seeds(1L);
  restrictcdf::cdf<restrictcdf::deriv>::set_working_memory(
    lower.n_elem, 1L);

  auto res = ([&](){
    if(derivs)
      return restrictcdf::cdf<restrictcdf::deriv>
      (lower, upper, mu, Sigma, do_reorder).approximate(
          maxvls, abseps, releps);

    return restrictcdf::cdf<restrictcdf::likelihood>
      (lower, upper, mu, Sigma, do_reorder).approximate(
          maxvls, abseps, releps);
  })();

  Rcpp::NumericVector out;
  out = res.finest;
  out.attr("minvls") = res.minvls;
  out.attr("inform") = res.inform;
  out.attr("abserr") = res.abserr;

  return out;
}

// [[Rcpp::export("get_commutation", rng = false)]]
arma::mat get_commutation_to_R(unsigned const n, unsigned const m){
  return get_commutation(n, m);
}

// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector get_commutation_vec
  (unsigned const n, unsigned const m, bool const transpose){
  std::unique_ptr<size_t[]> res =
    get_commutation_unequal_vec(n, m, transpose);
  Rcpp::IntegerVector out(n * m);
  for(size_t i = 0; i < n * m; ++i)
    out[i] = *(res.get() + i) + 1L;

  return out;
}

// [[Rcpp::export("x_dot_X_kron_I", rng = false)]]
arma::mat x_dot_X_kron_I_wrap
  (arma::vec const &x, arma::mat const &X, size_t const l){
  arma::mat out(1L, X.n_cols * l);
  x_dot_X_kron_I(x, X, l, out.memptr());
  return out;
}

