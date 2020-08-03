#include "logLik.h"
#include "threat-safe-random.h"
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "prof.h"

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
// [[Rcpp::export]]
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
    size_t const n_threads, bool const comp_derivs){
  profiler pp("eval_log_lm_terms");

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
        vcov, my_derivs, maxpts, abseps, releps, comp_derivs);

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
