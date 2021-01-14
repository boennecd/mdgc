#include "logLik.h"
#include "threat-safe-random.h"
#include "pnorm.h"
#include "qnorm.h"
#include "fast-commutation.h"
#include "restrict-cdf.h"
#include <limits.h>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "lp_utils.h"
#include "multinomial-probit.h"
#include <map>

using namespace mdgc;
using namespace arma;

struct ml_terms {
  // number of latent variables
  size_t n_variables;
  // vector with zero-based indices of the non-zero mean terms.
  arma::uvec idx_non_zero_mean;
  // log marginal likelihood term objects
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
 * @param multinomial list with 3xn matrix with multinomial outcomes. The
 * first index is the outcome, the second index is the number of categories,
 * and the third index is the index of the first latent variable.
 * @param idx_non_zero_mean indices for non-zero mean variables. Indices
 * should be sorted.
 *
 * Indices are zero-based except the outcome index for multinomial
 * variables.
 */
// [[Rcpp::export(rng = false)]]
SEXP get_log_lm_terms_cpp(arma::mat const &lower, arma::mat const &upper,
                          arma::imat const &code, Rcpp::List multinomial,
                          arma::uvec const idx_non_zero_mean){
  auto out = Rcpp::XPtr<ml_terms>(new ml_terms());

  size_t const n = lower.n_cols,
               p = lower.n_rows;
  if(upper.n_rows != p or upper.n_cols != n)
    throw std::invalid_argument("get_log_lm_terms: invalid 'upper'");
  if(code.n_rows != p or code.n_cols != n)
    throw std::invalid_argument("get_log_lm_terms: invalid 'code'");
  if(static_cast<size_t>(multinomial.size()) != n)
    throw std::invalid_argument("get_log_lm_terms: invalid 'multinomial'");
  if(idx_non_zero_mean.size() > 0 and idx_non_zero_mean[0] >= p)
    throw std::invalid_argument("get_log_lm_terms: invalid 'idx_non_zero_mean'");
  if(idx_non_zero_mean.size() > 1)
    for(size_t i = 1; i < idx_non_zero_mean.size(); ++i)
      if(idx_non_zero_mean[i] <= idx_non_zero_mean[i - 1] or
           idx_non_zero_mean[i] >= p)
        throw std::invalid_argument(
            "get_log_lm_terms: invalid 'idx_non_zero_mean' (index is too large, duplicate, or not sorted)");

  /* fill in the log ml terms objects */
  std::vector<log_ml_term> &terms = out->terms;
  out->idx_non_zero_mean = idx_non_zero_mean;
  out->n_variables = p;
  terms.reserve(n);
  uvec w_idx_int(p),
       w_idx_obs(p);
  vec w_obs_val(p), w_upper(p), w_lower(p);
  arma::imat cate_arg;

  for(size_t i = 0; i < n; ++i){
    size_t n_o(0L), n_i(0L);
    cate_arg = Rcpp::as<arma::imat>(multinomial[i]);

    for(size_t j = 0; j < p; ++j){
      if       (code.at(j, i) == 0L){
        // observed value
        w_idx_obs.at(n_o  ) = j;
        w_obs_val.at(n_o++) = upper.at(j, i);
      } else if(code.at(j, i) == 1L) {
        // check if we need to remove a multinomial outcome
        for(size_t k = 0; k < cate_arg.n_cols; ++k)
          if(static_cast<size_t>(cate_arg.at(2, k)) == j){
            cate_arg.shed_col(k);
            break;
          }
        // do nothing. The value is missing.
      } else if(code.at(j, i) == 2L){
        // Z is in an interval
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

    terms.emplace_back(a_idx_int, a_idx_obs, a_obs_val, a_lower, a_upper,
                       cate_arg);
  }

  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector eval_log_lm_terms(
    SEXP ptr, arma::ivec const &indices, arma::mat const &vcov,
    arma::vec const &mu,
    int const maxpts, double const abs_eps, double const rel_eps,
    size_t const n_threads, bool const comp_derivs, unsigned const minvls,
    bool const do_reorder = true, bool const use_aprx = false){
  Rcpp::XPtr<ml_terms> obj(ptr);
  std::vector<log_ml_term> const &terms = obj->terms;

  size_t const p = obj->n_variables;
  arma::uvec const &idx_non_zero_mean = obj->idx_non_zero_mean;
  if(vcov.n_cols != p or vcov.n_rows != p)
    throw std::invalid_argument("eval_log_lm_terms: invalid vcov");
  if(n_threads < 1L)
    throw std::invalid_argument("eval_log_lm_terms: invalid n_threads");
  if(mu.size() != idx_non_zero_mean.size())
    throw std::invalid_argument("eval_log_lm_terms: invalid mu");

  log_ml_term::set_working_memory(terms, n_threads);
  arma::vec mu_arg(p, arma::fill::zeros);
  mu_arg(idx_non_zero_mean) = mu;

  mat derivs_vcov =
    comp_derivs ? mat(p, p, arma::fill::zeros) : mat();
  arma::vec derivs_mea =
    comp_derivs ? vec(p   , arma::fill::zeros) : vec();

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
    arma::mat my_derivs_vcov(comp_derivs ? p : 0L, comp_derivs ? p : 0L,
                             arma::fill::zeros);
    arma::vec my_derivs_mea(comp_derivs ? p : 0L, arma::fill::zeros);

#ifdef _OPENMP
#pragma omp for schedule(static) reduction(+:out) nowait
#endif
    for(size_t i = 0; i < n_indices; ++i){
      if(static_cast<size_t>(indices[i]) >= terms.size())
        continue;
      out += terms[indices[i]].approximate(
        vcov, mu_arg, my_derivs_vcov, my_derivs_mea, maxpts, abs_eps,
        rel_eps, comp_derivs, do_reorder, minvls, use_aprx);
    }

    if(comp_derivs){
#ifdef _OPENMP
#pragma omp critical(add_derivs)
{
#endif
      derivs_vcov += my_derivs_vcov;
      derivs_mea  += my_derivs_mea;
#ifdef _OPENMP
}
#endif
    }
  }

  Rcpp::NumericVector res(1);
  res[0] = out;
  if(comp_derivs){
    res.attr("grad_vcov") = derivs_vcov;
    res.attr("grad_mea")  = vec(derivs_mea(idx_non_zero_mean));
  }

  return res;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix get_z_hat
(arma::mat const &lower, arma::mat const &upper, arma::imat const &code,
 unsigned const n_threads, Rcpp::List multinomial){
  size_t const p = lower.n_rows,
               n = upper.n_cols;
  if(upper.n_rows != p or upper.n_cols != n)
    throw std::invalid_argument("get_z_hat: invalid upper");
  if(code.n_rows != p or code.n_cols != n)
    throw std::invalid_argument("get_z_hat: invalid lower");
  if(n_threads < 1)
    throw std::invalid_argument("get_z_hat: invalid n_threads");
  if(static_cast<size_t>(multinomial.size()) != n)
    throw std::invalid_argument("get_z_hat: invalid multinomial");

  // check if there are any multinomial variables
  bool any_cate(false);
  for(auto ca : multinomial)
    if(Rcpp::as<arma::imat>(ca).n_elem > 0){
      any_cate = true;
      break;
    }

  Rcpp::NumericMatrix out(p, n);
  double * const o = &out[0];
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_threads) schedule(static) if(!any_cate)
#endif
  for(size_t j = 0; j < n; ++j){
    arma::imat multinomial_j;
    if(any_cate)
      multinomial_j = Rcpp::as<arma::imat>(multinomial[j]);

    size_t k(0);
    double * oj = o + j * p;
    for(size_t i = 0; i < p; ++i, ++oj){
      if(any_cate and k < multinomial_j.n_cols and
           i == static_cast<size_t>(multinomial_j.at(2, k))){
        /** we will set the observed level value to the median on the latent
         *  scale conditional on it being the greatest value. We set the
         *  other values to the median of a truncated normal distribution.
         */
        int const obs_lvl = multinomial_j.at(0, k) - 1,
              idx_obs_lvl = obs_lvl + i,
                   n_lvls = multinomial_j.at(1, k);
        bool const is_first = obs_lvl == 0L;

        if(code.at(i, j) == 1L)
          // the value is missing
          for(int l = 0; l < n_lvls; ++l, ++i, ++oj)
            *oj = upper.at(i, j);
        else {
          double const val_obs = qnorm_w(
            static_cast<double>(n_lvls) / (n_lvls + 1.),
            -upper.at(idx_obs_lvl, j), is_first ? 1e-8 : 1., 1L, 0L);

          *oj++ = 0; // always zero for identification
          ++i;

          for(int l = 1; l < n_lvls; ++l, ++i, ++oj)
            if(i == static_cast<size_t>(idx_obs_lvl))
              // the observed level
              *oj = val_obs;
            else {
              // median of the truncated distribution
              double const mu = -upper.at(i, j);

              *oj = qnorm_w(
                pnorm_std((val_obs - mu), 1L, 0L) / 2., mu, 1., 1L, 0L);

            }
        }

        ++k;
        --oj; // increment in the for-loop
        --i ; // increment in the for-loop
        continue;
      }

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
   arma::mat const &Sigma, int const maxvls, double const abs_eps,
   double const rel_eps, bool const derivs, bool const do_reorder = true,
   bool const use_aprx = false){
  parallelrng::set_rng_seeds(1L);

  if(derivs){
    restrictcdf::deriv::alloc_mem(lower.n_elem, 1L);
    restrictcdf::deriv functor(mu, Sigma);

    auto res = restrictcdf::cdf<restrictcdf::deriv>
    (functor, lower, upper, mu, Sigma, do_reorder, use_aprx).approximate(
        maxvls, abs_eps, rel_eps);

    Rcpp::NumericVector out(res.derivs.n_elem + 1);
    out[0] = res.likelihood;
    std::copy(res.derivs.begin(), res.derivs.end(), &out[1]);
    out.attr("minvls") = res.minvls;
    out.attr("inform") = res.inform;
    out.attr("abserr") = res.abserr;
    return out;
  }

  restrictcdf::likelihood::alloc_mem(lower.n_elem, 1L);
  restrictcdf::likelihood functor;
  auto res = restrictcdf::cdf<restrictcdf::likelihood>
    (functor, lower, upper, mu, Sigma, do_reorder, use_aprx).approximate(
        maxvls, abs_eps, rel_eps);

  Rcpp::NumericVector out = Rcpp::NumericVector::create(res.likelihood);
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

using contin = restrictcdf::imputation::contin;
using ordinal = restrictcdf::imputation::ordinal;
using known = restrictcdf::imputation::known;
using binary = restrictcdf::imputation::binary;
using multinomial_impu = restrictcdf::imputation::multinomial;
using impute_base = restrictcdf::imputation::type_base;

inline int impute_get_output_dim(contin const &x){
  return 1L;
}

inline int impute_get_output_dim(ordinal const &x){
  return x.n_bs + 1L;
}

inline int impute_get_output_dim(binary const &x){
  return 2L;
}

inline int impute_get_output_dim(multinomial_impu const &x){
  return x.n_lvls;
}

inline int impute_get_output_dim(impute_base const *type_base){
  contin const *c = dynamic_cast<contin const *>(type_base);
  if(c)
    return impute_get_output_dim(*c);
  ordinal const *o = dynamic_cast<ordinal const *>(type_base);
  if(o)
    return impute_get_output_dim(*o);
  binary const *b = dynamic_cast<binary const *>(type_base);
  if(b)
    return impute_get_output_dim(*b);
  multinomial_impu const *m =
    dynamic_cast<multinomial_impu const *>(type_base);
  if(m)
    return impute_get_output_dim(*m);

  throw std::invalid_argument("impute_get_output_dim: not implemented");
}

inline void impute_set_val
  (contin const &x, double *&res, double const *&val){
  *res++ = *val++;
}

inline void impute_set_val
  (ordinal const &x, double *&res, double const *&val){
  double *r = res;
  size_t const nvals = x.n_bs + 1L;
  double su(0);

  for(size_t i = 0; i < nvals; ++i, ++res, ++val){
    su += *val;
    *res = *val;
  }

  // normalize
  for(size_t i = 0; i < nvals; ++i, ++r)
    *r /= su;
}

inline void impute_set_val
  (multinomial_impu const &x, double *&res, double const *&val){
  double *r = res;
  size_t const nvals = x.n_lvls;
  double su(0);

  for(size_t i = 0; i < nvals; ++i, ++res, ++val){
    su += *val;
    *res = *val;
  }

  // normalize
  for(size_t i = 0; i < nvals; ++i, ++r)
    *r /= su;
}

inline void impute_set_val
  (binary const &x, double *&res, double const *&val){
  double su(0);
  su          += *val;
  *res         = *val++;
  su          += *val;
  *(res + 1L)  = *val++;
  *res++      /= su;
  *res++      /= su;
}

inline void impute_set_val(impute_base const *type_base,
                           double *&res, double const *&val){
  contin const *c = dynamic_cast<contin const *>(type_base);
  if(c)
    return impute_set_val(*c, res, val);
  ordinal const *o = dynamic_cast<ordinal const *>(type_base);
  if(o)
    return impute_set_val(*o, res, val);
  binary const *b = dynamic_cast<binary const *>(type_base);
  if(b)
    return impute_set_val(*b, res, val);
  multinomial_impu const *m =
    dynamic_cast<multinomial_impu const *>(type_base);
  if(m)
    return impute_set_val(*m, res, val);

  throw std::invalid_argument("impute_set_val: not implemented");
}

inline SEXP impute_set_val_R
  (contin const &x, double const *&val, Rcpp::CharacterVector names,
   Rcpp::Function marg, int const code, double const truth){

  Rcpp::NumericVector out(1L);
  if(code != 1L)
    out[0L] = truth;
  else {
    out[0L] = pnorm_std(*val, 1L, 0L);
    out = Rcpp::NumericVector(marg(out));
    out.attr("names") = R_NilValue;
  }
  ++val;
  return out;
}

inline SEXP impute_set_val_R
  (ordinal const &x, double const *&val, Rcpp::CharacterVector names,
   Rcpp::Function marg, int const code, double const truth){
  size_t const n_ele = x.n_bs + 1L;
  Rcpp::NumericVector out(n_ele);

  if(code != 1L){
    // the value is observed
    out[std::lround(truth) - 1] = 1.;
    val += n_ele;

  } else
    // the value is not observed
    for(size_t i = 0; i < n_ele; ++i, ++val)
      out[i] = *val;

  out.attr("names") = names;
  return out;
}

inline SEXP impute_set_val_R
  (multinomial_impu const &x, double const *&val, Rcpp::CharacterVector names,
   Rcpp::Function marg, int const code, double const truth){
  size_t const n_ele = x.n_lvls;
  Rcpp::NumericVector out(n_ele);

  if(code != 1L){
    // the value is observed
    out[std::lround(truth) - 1] = 1.;
    val += n_ele;

  } else
    // the value is not observed
    for(size_t i = 0; i < n_ele; ++i, ++val)
      out[i] = *val;

  out.attr("names") = names;
  return out;
}

inline SEXP impute_set_val_R
  (binary const &x, double const *&val, Rcpp::CharacterVector names,
   Rcpp::Function marg, int const code, double const truth){
  Rcpp::NumericVector out(2L);
  if(code != 1L){
    out[std::lround(truth)] = 1;
    val += 2L;

  } else {
    out[0] = *val++;
    out[1] = *val++;
  }
  out.attr("names") = names;
  return out;
}

inline SEXP impute_set_val_R
  (impute_base const *type_base, double const *&val,
   Rcpp::CharacterVector names, Rcpp::Function marg,
   int const code, double const truth){
  contin const *c = dynamic_cast<contin const *>(type_base);
  if(c)
    return impute_set_val_R(*c, val, names, marg, code, truth);
  ordinal const *o = dynamic_cast<ordinal const *>(type_base);
  if(o)
    return impute_set_val_R(*o, val, names, marg, code, truth);
  binary const *b = dynamic_cast<binary const *>(type_base);
  if(b)
    return impute_set_val_R(*b, val, names, marg, code, truth);
  multinomial_impu const *m =
    dynamic_cast<multinomial_impu const *>(type_base);
  if(m)
    return impute_set_val_R(*m, val, names, marg, code, truth);

  throw std::invalid_argument("impute_set_val_R: not implemented");
}

// [[Rcpp::export()]]
Rcpp::List impute
  (arma::mat const &lower, arma::mat const &upper, arma::imat const &code,
   arma::mat const &Sigma, arma::vec const &mea,
   arma::mat const &truth, Rcpp::List margs, Rcpp::List multinomial,
   double const rel_eps, double const abs_eps, unsigned const maxit,
   Rcpp::List passed_names, Rcpp::CharacterVector outer_names,
   int const n_threads, bool const do_reorder, int const minvls,
   bool const use_aprx = false){
  // setup vector to pass to RQMC method
  std::vector<std::unique_ptr<impute_base> > const type_list = ([&](){
    std::vector<std::unique_ptr<impute_base> > out;
    out.reserve(margs.size());
    for(auto o : margs){
      Rcpp::RObject ro(o),
                    bo(ro.attr("borders")),
                    mu(ro.attr("means"));

      if(!mu.isNULL())
        // the variable is multinomial
        out.emplace_back(new restrictcdf::imputation::multinomial(
            Rcpp::NumericVector(mu).size() + 1));
      else if(bo.isNULL())
        // the variable is continous and we just have to increament the
        // pointer
        out.emplace_back(new restrictcdf::imputation::contin());
      else {
        // the variable is ordinal or binary
        Rcpp::NumericVector bo_vec(bo);
        if(bo_vec.size() > 3L)
          // the variable is ordinal
          out.emplace_back(new ordinal(&bo_vec[0], bo_vec.size()));
        else
          // the variable is binary
          out.emplace_back(new binary(0));
      }
    }

    return out;
  })();

  // vector with pointer to the above
  std::vector<impute_base const *> const type_list_ptr = ([&](){
    std::vector<impute_base const *> out;
    out.reserve(type_list.size());
    for(auto &x : type_list)
      out.emplace_back(x.get());
    return out;
  })();

  // find values using QMC
  size_t const out_dim =
    restrictcdf::imputation::get_n_integrands(type_list_ptr) - 1;

  size_t const n_obs = lower.n_cols;
  arma::mat out_dat(out_dim, n_obs);
  int const n_vars = lower.n_rows;

  restrictcdf::imputation::alloc_mem(type_list_ptr, n_threads);

#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif
  parallelrng::set_rng_seeds(n_threads);

  // thread private variables
  uvec w_idx_int        (n_vars),
       w_idx_obs        (n_vars),
       w_idx_cat_obs    (n_vars),
       w_idx_cat_not_obs(n_vars);
  vec w_obs_val(n_vars),
        w_upper(n_vars),
        w_lower(n_vars),
             mu(n_vars);

  uvec a_idx_int, a_idx_obs, a_idx_cat_obs, a_idx_cat_not_obs;
  vec a_obs_val, a_upper, a_lower, mu_use;
  mat D;

  std::vector<impute_base const *> type_i;

  std::map<int, known> known_objs;
  known_objs.insert(std::pair<int, known>(1L, known(1L)));

  std::vector<imat> cate_mat;
  cate_mat.reserve(multinomial.size());
  for(auto &cate : multinomial){
    cate_mat.emplace_back(Rcpp::as<imat>(cate));
    arma::imat &new_ele = cate_mat.back();
    for(arma::uword k = 0; k < new_ele.n_cols; ++k){
      int const n_latent_k = new_ele.at(1, k) - 1L;
      known_objs.insert(
        std::pair<int, known>(n_latent_k, n_latent_k));
    }
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) \
  firstprivate(w_idx_int, w_idx_obs, w_obs_val, w_upper, w_lower, \
               mu, w_idx_cat_obs, w_idx_cat_not_obs, known_objs) \
  private(a_idx_int, a_idx_obs, a_obs_val, a_upper, a_lower, \
          mu_use, type_i, a_idx_cat_obs, a_idx_cat_not_obs, D)
#endif
  for(size_t i = 0; i < n_obs; ++i){
    // create covariance matrix and bounds to pass
    bool any_missing(false);
    int n_i_all   (0L),
           n_o    (0L),
         n_cat_obs(0L),
           n_i    (0L),
           k      (0L),
           n_D_col(0L);
    arma::imat const &cate_mat_i = cate_mat[i];
    constexpr double const inf =
      std::numeric_limits<double>::infinity();
    D.zeros(n_vars, cate_mat_i.n_cols);

    for(int j = 0; j < n_vars; ++j){
      if(static_cast<size_t>(k) < cate_mat_i.n_cols and
           j == cate_mat_i.at(2, k)){
        int const n_lvls = cate_mat_i.at(1, k);

        if(code.at(j, i) == 1L){
          // no information (missing and needs imputation)
          for(int l = 0; l < n_lvls; ++l, ++j){
            w_idx_int.at(n_i_all  ) = j;
            mu       .at(n_i_all++) = mea[j];

            w_idx_cat_not_obs.at(n_i  ) = j - n_o;
            w_lower          .at(n_i  ) = -inf;
            w_upper          .at(n_i++) =  inf;
          }
          any_missing = true;

        } else {
          // observed the level
          int const idx_obs_lvl = cate_mat_i.at(0, k) - 1 + j;

          for(int l = 0; l < n_lvls; ++l, ++j){
            w_idx_int.at(n_i_all  ) = j;
            mu       .at(n_i_all++) = mea[j];

            if(j != idx_obs_lvl){
              D                .at(n_i, n_D_col) = 1;
              w_idx_cat_not_obs.at(n_i         ) = j - n_o;
              w_lower          .at(n_i         ) = lower.at(j, i);
              w_upper          .at(n_i++       ) =
                upper.at(j, i) - upper.at(idx_obs_lvl, i);
            } else
              w_idx_cat_obs[n_cat_obs++] = j - n_o;
          }

          ++n_D_col;
        }

        ++k;
        --j; // incremented in the for-loop
        continue;
      }

      if       (code.at(j, i) == 0L){
        /* observed value */
        w_idx_obs.at(n_o  ) = j;
        w_obs_val.at(n_o++) = upper.at(j, i);

      } else if(code.at(j, i) == 1L) {
        /* no information (missing and needs imputation) */
        w_idx_int.at(n_i_all  ) = j;
        mu       .at(n_i_all++) = mea[j];

        w_idx_cat_not_obs.at(n_i  ) = j - n_o;
        w_lower          .at(n_i  ) = -inf;
        w_upper          .at(n_i++) =  inf;
        any_missing = true;

      } else if(code.at(j, i) == 2L){
        /* Z is in an interval */
        w_idx_int.at(n_i_all  ) = j;
        mu       .at(n_i_all++) = mea[j];

        w_idx_cat_not_obs.at(n_i  ) = j - n_o;
        w_lower          .at(n_i  ) = lower.at(j, i);
        w_upper          .at(n_i++) = upper.at(j, i);

      } else
        throw std::invalid_argument("impute: invalid code");
    }

    if(!any_missing)
      continue;

    if(n_i > 0){
      a_idx_int = w_idx_int.subvec(0L, n_i_all - 1L);
      a_upper   = w_upper  .subvec(0L, n_i     - 1L);
      a_lower   = w_lower  .subvec(0L, n_i     - 1L);

    } else
      throw std::runtime_error("should not be reached");

    if(n_o > 0){
      a_idx_obs = w_idx_obs.subvec(0L, n_o - 1L);
      a_obs_val = w_obs_val.subvec(0L, n_o - 1L);

    } else {
      a_idx_obs.set_size(0L);
      a_obs_val.set_size(0L);

    }

    if(cate_mat_i.n_cols > 0){
      if(n_cat_obs > 0)
        a_idx_cat_obs = w_idx_cat_obs.subvec(0L, n_cat_obs - 1L);
      else
        a_idx_cat_obs.set_size(0L);
      a_idx_cat_not_obs = w_idx_cat_not_obs.subvec(0L, n_i - 1L);

      if(n_D_col > 0)
        D = D.submat(0, 0, n_i - 1, n_D_col - 1);
      else
        D.set_size(0, 0);
    }

    mu_use = mu.subvec(0L, n_i_all - 1L);

    mat Sigma_use;
    if(n_o > 0){
      // compute conditional mean and covariance matrix
      arma::mat const S_oo = Sigma(a_idx_obs, a_idx_obs),
                      S_oi = Sigma(a_idx_obs, a_idx_int),
                      Mtmp = arma::solve(S_oo, S_oi,
                                         arma::solve_opts::likely_sympd);
      Sigma_use = Sigma(a_idx_int, a_idx_int);

      Sigma_use -= S_oi.t() * Mtmp;
      mu_use += Mtmp.t() * a_obs_val;

    } else
      // no observed continous variables to condition on
      Sigma_use = Sigma;

    if(n_D_col > 0){
      arma::mat const D_Sigma =
        D * Sigma_use(a_idx_cat_obs, a_idx_cat_not_obs);
      Sigma_use = Sigma_use(a_idx_cat_not_obs, a_idx_cat_not_obs) +
        D * Sigma_use(a_idx_cat_obs, a_idx_cat_obs) * D.t() -
        D_Sigma - D_Sigma.t();

      mu_use = mu_use(a_idx_cat_not_obs) - D * mu_use(a_idx_cat_obs);
    }

    // create the type objects to pass
    type_i.clear();
    type_i.reserve(n_i);
    {
      arma::uword const *cur_idx_obs = a_idx_obs.begin();
      int j = 0;
      k = 0;
      for(auto &type_list_it : type_list){
        if(cur_idx_obs != a_idx_obs.end() and
             static_cast<size_t>(j) == *cur_idx_obs){
          // the outcome is known and continous so we can use the value
          // as is
          ++cur_idx_obs;
          ++j;
          continue;
        }

        /* we know it is a latent variable which is multinomial, binary or
         * ordinal. Thus, there will be a latent variable regardless of
         * whether we observe the value or not. */
        if(static_cast<size_t>(k) < cate_mat_i.n_cols and
             j == cate_mat_i.at(2, k)){
          // we have a multinomial outcome
          int const n_lvls = cate_mat_i.at(1, k);

          if(code.at(j, i) == 2L) {
            // the level is known. The number of variables we need to
            // integrate out is one less
            type_i.emplace_back(&known_objs.at(n_lvls - 1L));
            ++j; // because of the n_lvls - 1L
          } else
            // the level is unknown
            type_i.emplace_back(type_list_it.get());

          j += type_i.back()->n_latent();
          continue;
        }

        if(code.at(j, i) == 2L)
          // the level is observed. Ordinal/binary so only one latent
          // variable
          type_i.emplace_back(&known_objs.at(1L));
        else
          type_i.emplace_back(type_list_it.get());
        j += type_i.back()->n_latent();
      }
    }

    restrictcdf::imputation imputer(type_i, mu_use, Sigma_use);
    auto res = restrictcdf::cdf<restrictcdf::imputation>
      (imputer, a_lower, a_upper, mu_use, Sigma_use, do_reorder,
       use_aprx).approximate(maxit, abs_eps, rel_eps, minvls);

    res.imputations /= res.likelihood;

    double *o = out_dat.colptr(i);
    double const *rval = res.imputations.begin();
    {
      arma::uword const *cur_idx_obs = a_idx_obs.begin();
      int j(0);
      for(auto &type_j : type_list){
        if(cur_idx_obs != a_idx_obs.end() and
             static_cast<size_t>(j) == *cur_idx_obs){
          // the variable is continous and known. Do nothing.
          ++cur_idx_obs;
          o += impute_get_output_dim(type_j.get());
          j += type_j->n_latent();
          continue;
        }

        if(code.at(j, i) == 2L)
          o += impute_get_output_dim(type_j.get());
        else
          impute_set_val(type_j.get(), o, rval);
        j += type_j->n_latent();
      }
    }
  }

  // format and return
  std::vector<Rcpp::Function> funcs;
  funcs.reserve(margs.size());
  for(auto m : margs)
    funcs.emplace_back(m);

  Rcpp::List out(n_obs);
  for(size_t i = 0; i < n_obs; ++i){
    double const *o = out_dat.colptr(i);
    int const *code_j = code.colptr(i);
    Rcpp::List out_i(type_list.size());
    for(size_t j = 0; j < type_list.size(); ++j){
      out_i[j] = impute_set_val_R(
        type_list[j].get(), o, passed_names[j], funcs[j],
        *code_j, truth.at(j, i));
      code_j += type_list[j]->n_latent();
    }

    out_i.attr("names") = outer_names;
    out[i] = out_i;
  }

  return out;
}

/**
 * This function computes the inner product of rows and columns of a
 * lower triangular matrix X. That is
 *
 *   c(X) = (x(i_1)^T x(j_1), ..., x(i_m)^T x(j_m))
 *
 * The functions takes a vector with the lower triangular entries. The
 * gradient with respect to rhs^Tc'(X) is returned if jacob is true.
 */
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector lower_tri_inner
  (Rcpp::NumericVector x, Rcpp::IntegerMatrix idx, bool const jacob,
   Rcpp::NumericVector rhs) {
  if(idx.nrow() < 1)
    return Rcpp::NumericVector();

  double const fdim = .5 * (std::sqrt(8 * x.size() + 1) - 1);
  int const dim = std::lround(fdim);
  if(std::abs(fdim / dim - 1) >
       std::numeric_limits<double>::epsilon() * 10)
    throw std::invalid_argument("lower_tri_outer: invalid x");
  if(idx.ncol() != 2L)
    throw std::invalid_argument("lower_tri_outer: invalid idx");
  if(jacob and rhs.size() != idx.nrow())
    throw std::invalid_argument("lower_tri_outer: invalid rhs");

  auto tri_map = [&](int const r, int const c){
    return r + c * dim - (c * (c + 1)) / 2;
  };

  if(jacob){
    Rcpp::NumericVector out(x.size());
    for(int i = 0; i < idx.nrow(); ++i){
      int const row_i = idx(i, 0),
                col_i = idx(i, 1);

      int const n_terms = std::min(row_i, col_i) + 1;
      int icol = tri_map(col_i, 0),
          irow = tri_map(row_i, 0);
      for(int j = 0; j < n_terms; ++j, icol += dim - j, irow += dim - j){
        out(icol) += x[irow] * rhs[i];
        out(irow) += x[icol] * rhs[i];
      }
    }

    return out;
  }

  Rcpp::NumericVector out(idx.nrow());
  for(int i = 0; i < idx.nrow(); ++i){
    double out_i(0.);
    int const row_i = idx(i, 0),
              col_i = idx(i, 1);
    int const n_terms = std::min(row_i, col_i) + 1;
    double const *dcol = &x[tri_map(col_i, 0)],
                 *drow = &x[tri_map(row_i, 0)];
    for(int j = 0; j < n_terms; ++j, dcol += dim - j, drow += dim - j)
      out_i += *drow * *dcol;
    out[i] = out_i;
  }

  return out;
}

// [[Rcpp::export(rng = false)]]
double eval_multinomial_prob(int const icase, arma::vec const &means){
  if(static_cast<size_t>(icase) >= means.size() + 1 or icase < 0)
    throw std::invalid_argument("eval_multinomial_prob: invalid icase");
  if(means.size() < 1)
    throw std::invalid_argument("eval_multinomial_prob: invalid means");

  return multinomial::eval(means.begin(), icase, means.size() + 1);
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector eval_multinomial_prob_gr(
    int const icase, arma::vec const &means){
  if(static_cast<size_t>(icase) >= means.size() + 1 or icase < 0)
    throw std::invalid_argument("eval_multinomial_prob: invalid icase");
  if(means.size() < 1)
    throw std::invalid_argument("eval_multinomial_prob: invalid means");

  Rcpp::NumericVector out(means.size());
  std::unique_ptr<double[]> wk(new double[means.size()]);
  out.attr("prob") = multinomial::eval_gr(
    means.begin(), &out[0], icase, means.size() + 1, wk.get());
  return out;
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector multinomial_find_means
  (arma::vec const &probs, double const rel_eps = 3.000214e-13,
   int const max_it = 100, double const c1 = .0001,
   double const c2 = .9){
  if(probs.size() < 2 or std::abs(arma::sum(probs) - 1) >= 1e-10)
    throw std::invalid_argument("multinomial_find_means: invalid probs");

  Rcpp::NumericVector mu(probs.size() - 1);
  mu.attr("info-code") = multinomial::find_means
    (probs.begin(), &mu[0], probs.size(), rel_eps, max_it, c1, c2);

  return mu;
}
