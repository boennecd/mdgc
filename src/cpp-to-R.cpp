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
SEXP get_log_lm_terms_cpp(arma::mat const &lower, arma::mat const &upper,
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
    int const maxpts, double const abs_eps, double const rel_eps,
    size_t const n_threads, bool const comp_derivs, unsigned const minvls,
    bool const do_reorder = true, bool const use_aprx = false){
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
#pragma omp for schedule(static) reduction(+:out) nowait
#endif
    for(size_t i = 0; i < n_indices; ++i){
      if(static_cast<size_t>(indices[i]) >= terms.size())
        continue;
      out += terms[indices[i]].approximate(
        vcov, my_derivs, maxpts, abs_eps, rel_eps, comp_derivs,
        do_reorder, minvls, use_aprx);
    }

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
#pragma omp parallel for num_threads(n_threads) schedule(static)
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
using impute_base = restrictcdf::imputation::type_base;


inline size_t impute_get_output_dim(contin const &x){
  return 1L;
}

inline size_t impute_get_output_dim(ordinal const &x){
  return x.n_bs + 1L;
}

inline size_t impute_get_output_dim(binary const &x){
  return 2L;
}

inline size_t impute_get_output_dim(impute_base const *type_base){
  contin const *c = dynamic_cast<contin const *>(type_base);
  if(c)
    return impute_get_output_dim(*c);
  ordinal const *o = dynamic_cast<ordinal const *>(type_base);
  if(o)
    return impute_get_output_dim(*o);
  binary const *b = dynamic_cast<binary const *>(type_base);
  if(b)
    return impute_get_output_dim(*b);

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

  throw std::invalid_argument("impute_set_val: not implemented");
}

inline SEXP impute_set_val_R
  (contin const &x, double const *&val, Rcpp::CharacterVector names,
   Rcpp::Function marg, double const lower, double const upper,
   int const code, double const truth){
  double const v_use = code != 1L ? upper : *val;
  val++;

  Rcpp::NumericVector out(1L);
  if(code != 1L)
    out[0L] = truth;
  else {
    out[0L] = pnorm_std(v_use, 1L, 0L);
    out = Rcpp::NumericVector(marg(out));
    out.attr("names") = R_NilValue;
  }
  return out;
}

inline SEXP impute_set_val_R
  (ordinal const &x, double const *&val, Rcpp::CharacterVector names,
   Rcpp::Function marg, double const lower, double const upper,
   int const code, double const truth){
  size_t const n_ele = x.n_bs + 1L;
  Rcpp::NumericVector out(n_ele);

  if(code != 1L){
    size_t i = 0;
    for(; i < n_ele - 1L; ++i)
      if(lower < *(x.borders.get() + i)){
        out[i] = 1.;
        break;
      }

    if(i == n_ele - 1L)
      out[i] = 1.;

    val += n_ele;

  } else {
    for(size_t i = 0; i < n_ele; ++i, ++val)
      out[i] = *val;
  }

  out.attr("names") = names;
  return out;
}

inline SEXP impute_set_val_R
  (binary const &x, double const *&val, Rcpp::CharacterVector names,
   Rcpp::Function marg, double const lower, double const upper,
   int const code, double const truth){
  Rcpp::NumericVector out(2L);
  if(code != 1L){
    if(lower < x.border){
      out[0] = 1.;
      out[1] = 0.;
    } else {
      out[0] = 0.;
      out[1] = 1.;
    }
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
   double const lower, double const upper, int const code,
   double const truth){
  contin const *c = dynamic_cast<contin const *>(type_base);
  if(c)
    return impute_set_val_R(*c, val, names, marg, lower, upper, code,
                            truth);
  ordinal const *o = dynamic_cast<ordinal const *>(type_base);
  if(o)
    return impute_set_val_R(*o, val, names, marg, lower, upper, code,
                            truth);
  binary const *b = dynamic_cast<binary const *>(type_base);
  if(b)
    return impute_set_val_R(*b, val, names, marg, lower, upper, code,
                            truth);

  throw std::invalid_argument("impute_set_val_R: not implemented");
}

// [[Rcpp::export()]]
Rcpp::List impute
  (arma::mat const &lower, arma::mat const &upper, arma::imat const &code,
   arma::mat const &Sigma, arma::mat const &truth, Rcpp::List margs,
   double const rel_eps, double const abs_eps, unsigned const maxit,
   Rcpp::List passed_names, Rcpp::CharacterVector outer_names,
   int const n_threads, bool const do_reorder, int const minvls,
   bool const use_aprx = false){
  // setup vector to pass to QMC method
  std::vector<std::unique_ptr<impute_base> > const type_list = ([&](){
    std::vector<std::unique_ptr<impute_base> > out;
    out.reserve(margs.size());
    for(auto o : margs){
      Rcpp::RObject ro(o),
                    bo(ro.attr("borders"));
      if(bo.isNULL())
        out.emplace_back(new restrictcdf::imputation::contin());
      else {
        Rcpp::NumericVector bo_vec(bo);
        if(bo_vec.size() > 3L)
          out.emplace_back(new ordinal(&bo_vec[0], bo_vec.size()));
        else
          out.emplace_back(new binary(bo_vec[1]));
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
  size_t const n_vars = lower.n_rows;

  restrictcdf::imputation::alloc_mem(type_list_ptr, n_vars, n_threads);

#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif
  parallelrng::set_rng_seeds(n_threads);

  // thread private variables
  uvec w_idx_int(n_vars),
       w_idx_obs(n_vars);
  vec w_obs_val(n_vars),
        w_upper(n_vars),
        w_lower(n_vars);

  uvec a_idx_int, a_idx_obs;
  vec a_obs_val, a_upper, a_lower, mu_use;

  std::vector<impute_base const *> type_i;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) \
  firstprivate(w_idx_int, w_idx_obs, w_obs_val, w_upper, w_lower) \
  private(a_idx_int, a_idx_obs, a_obs_val, a_upper, a_lower, \
          mu_use, type_i)
#endif
  for(size_t i = 0; i < n_obs; ++i){
    // create covariance matrix and bounds to pass
    bool any_missing(false);
    size_t n_o(0L),
           n_i(0L);
    for(size_t j = 0; j < n_vars; ++j){
      if       (code.at(j, i) == 0L){
        /* observed value */
        w_idx_obs.at(n_o  ) = j;
        w_obs_val.at(n_o++) = upper.at(j, i);

      } else if(code.at(j, i) == 1L) {
        /* no information (missing and needs imputation) */
        w_idx_int.at(n_i  ) = j;
        double const inf = std::numeric_limits<double>::infinity();
        w_lower  .at(n_i  ) = -inf;
        w_upper  .at(n_i++) =  inf;
        any_missing = true;

      } else if(code.at(j, i) == 2L){
        /* Z is in an interval */
        w_idx_int.at(n_i  ) = j;
        w_lower  .at(n_i  ) = lower.at(j, i);
        w_upper  .at(n_i++) = upper.at(j, i);

      } else
        throw std::invalid_argument("impute: invalid code");
    }

    if(!any_missing)
      continue;

    if(n_i > 0){
      a_idx_int = w_idx_int.subvec(0L, n_i - 1L);
      a_upper   = w_upper  .subvec(0L, n_i - 1L);
      a_lower   = w_lower  .subvec(0L, n_i - 1L);

    } else
      throw std::runtime_error("should not be reached");

    if(n_o > 0){
      a_idx_obs = w_idx_obs.subvec(0L, n_o - 1L);
      a_obs_val = w_obs_val.subvec(0L, n_o - 1L);

    } else {
      a_idx_obs.set_size(0L);
      a_obs_val.set_size(0L);

    }

    mat Sigma_use;
    if(n_o > 0){
      // compute conditional mean and covariance matrix
      arma::mat const S_oo = Sigma(a_idx_obs, a_idx_obs),
                      S_oi = Sigma(a_idx_obs, a_idx_int),
                      Mtmp = arma::solve(S_oo, S_oi,
                                         arma::solve_opts::likely_sympd);
      Sigma_use = Sigma(a_idx_int, a_idx_int);

      Sigma_use -= S_oi.t() * Mtmp;
      mu_use = Mtmp.t() * a_obs_val;

    } else {
      mu_use.zeros(n_i);
      Sigma_use = arma::mat(
        const_cast<double*>(Sigma.begin()), Sigma.n_rows, Sigma.n_cols,
        false);

    }

    // create the type objects to pass
    type_i.clear();
    type_i.reserve(n_i);
    known known_obj;
    {
      arma::uword const *cur_idx_obs = a_idx_obs.begin();
      size_t j = 0;
      for(; j < n_vars; ++j){
        if(cur_idx_obs != a_idx_obs.end() and j == *cur_idx_obs){
          ++cur_idx_obs;
          continue;
        }

        if(code.at(j, i) == 2L){
          type_i.emplace_back(&known_obj);
        } else
          type_i.emplace_back(type_list[j].get());
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
      size_t j = 0;
      for(; j < n_vars; ++j){
        auto &type_j = type_list[j];

        if(cur_idx_obs != a_idx_obs.end() and j == *cur_idx_obs){
          ++cur_idx_obs;
          o += impute_get_output_dim(type_j.get());
          continue;
        }

        if(code.at(j, i) == 2L)
          o += impute_get_output_dim(type_j.get());
        else
          impute_set_val(type_j.get(), o, rval);
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
    Rcpp::List out_i(n_vars);
    for(size_t j = 0; j < n_vars; ++j)
      out_i[j] = impute_set_val_R(
        type_list[j].get(), o, passed_names[j], funcs[j],
        lower.at(j, i), upper.at(j, i), code.at(j, i), truth.at(j, i));

    out_i.attr("names") = outer_names;
    out[i] = out_i;
  }

  return out;
}

/**
 * This function computes the outer product of rows and columns of a
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
