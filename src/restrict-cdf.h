#ifndef RESTRICT_CDF_H
#define RESTRICT_CDF_H

#include "arma-wrap.h"
#include <array>
#include <limits>
#include <memory>
#include <cmath>
#include "pnorm.h"
#include "qnorm.h"
#include "mvtnorm-wrapper.h"
#include "config.h"
#include <algorithm>
#include "new-mvt.h"
#include "norm-cdf-approx.h"
#include <cmath>
#include <stdexcept>
#include "mdgc-mem.h"

namespace restrictcdf {
extern "C"
{
  /**
   * @param N Dimension of the integral.
   * @param lower N-vector with lower bounds.
   * @param upper N-vector with upper bounds.
   * @param delta N-vector with mean.
   * @param correl N(N - 1)/2-dimensional vector with  upper triangle of the
   * correlation matrix.
   * @param infin N-dimensional vector indicating whether the bounds are
   * finite.
   * @param pivot not sure. Set it to true.
   * @param y N-dimensional vector with workig memory.
   * @param ND N unless there is double infinite regions.
   * @param A potentially permutated version of lower.
   * @param B potentially permutated version of upper.
   * @param DL potentially permutated version of delta.
   * @param cov N(N + 1)/2-dimensional vector with potentially permutated
   * Cholesky decomposition of correl.
   * @param infi potentially permutated version of infin.
   * @param inform non-zero if something went wrong.
   * @param idx N-dimensional vector with indices of applied permutation.
   * @param doscale logical for whether to scale the cholesky decomposition
   * to have ones in the diagonal.
   *
   * cov is scaled such that the diagonal entries are one. This implies that
   * it is __not__ the Cholesky decomposition of the correlation matrix.
   * A, B, and DL are scaled accordingly.
   */
  void F77_NAME(mvsort)(
      int const* /* N */, double const* /* lower */,
      double const* /* upper */, double const* /* delta */,
      double const* /* correl */, int const* /* infin */,
      double const* /* y */, int const* /* pivot */,
      int* /* ND */, double* /* A */, double* /* B */, double* /* DL */,
      double* /* cov */, int* /* infi */, int* /* inform */,
      int* /* idx */, int const* /* doscale */);
}

inline double safe_qnorm_w(double const x) noexcept {
  constexpr double const eps =
    std::numeric_limits<double>::epsilon() *
    std::numeric_limits<double>::epsilon() *
    std::numeric_limits<double>::epsilon();
  if(x <= 0)
    return  qnorm_w(eps, 0, 1, 1L, 0L);
  else if(x >= 1.)
    return -qnorm_w(eps, 0, 1, 1L, 0L);

  return qnorm_w   (x  , 0, 1, 1L, 0L);
}

inline double safe_qnorm_aprx(double const x) noexcept {
  constexpr double const eps =
    std::numeric_limits<double>::epsilon() *
    std::numeric_limits<double>::epsilon() *
    std::numeric_limits<double>::epsilon();
  if(x <= 0)
    return  qnorm_aprx(eps);
  else if(x >= 1.)
    return -qnorm_aprx(eps);

  return qnorm_aprx   (x);
}

/**
 * copies the upper upper triangular matrix.
 *
 * @param X Matrix top copy.
 * @param x Pointer to copy to.
 */
inline void copy_upper_tri
  (arma::mat const &X, double * __restrict__ x) noexcept {
  int const p = X.n_cols;
  for(int c = 0; c < p; c++)
    for(int r = 0; r <= c; r++, x++)
      *x = X.at(r, c);
}

inline void copy_lower_tri
  (arma::mat const &X, double * __restrict__ x) noexcept {
  int const p = X.n_cols;
  for(int c = 0; c < p; c++)
    for(int r = c; r < p; r++, x++)
      *x = X.at(r, c);
}

/**
 * TODO: describe what this class does.
 */
template<class T_Functor, class out_type = typename T_Functor::out_type>
class cdf {
  T_Functor &functor;
  int const ndim, n_integrands;
  double * wk_mem = nullptr;
  int * iwk_mem = nullptr;
  bool is_permutated = false,
       use_aprx;
  constexpr bool needs_last_unif() const {
    return T_Functor::needs_last_unif();
  }

  // cached memory to use
  static cache_mem<int   > imem;
  static cache_mem<double> dmem;

  arma::ivec infin;
  arma::ivec indices;

  arma::vec lower = arma::vec(dmem.get_mem() , ndim, false),
            upper = arma::vec(lower.end()    , ndim, false),
       sigma_chol = arma::vec(upper.end()    , (ndim * (ndim + 1L)) / 2L,
                           false),
                           draw = arma::vec(sigma_chol.end(), ndim, false);
  // memory that can be used
  int * const itmp_mem = indices.end();
  double * const dtmp_mem = draw.end();

public:
  /**
   * must be called perior to calling the constructor or any member
   * functions.
   */
  static void alloc_mem(int const max_ndim, int const max_threads) {
    int const n_up_tri = (max_ndim * (max_ndim + 1)) / 2;
    imem.set_n_mem(3 * max_ndim                                 ,
                   max_threads);
    dmem.set_n_mem(7 * max_ndim + n_up_tri + max_ndim * max_ndim,
                   max_threads);
  }

  /**
   * @functor T_Functor class to use.
   * @param lower_in Lower bounds in the CDF.
   * @param upper_in Upper bounds in the CDF.
   * @param mu_in Mean vector.
   * @param sigma_in Covariance matrix.
   * @param do_reorder true if the order of integrations may be reordered.
   * @param use_aprx true if an approximation of the normal CDF should
   * be used.
   */
  cdf(T_Functor &functor, arma::vec const &lower_in,
      arma::vec const &upper_in, arma::vec const &mu_in,
      arma::mat const &sigma_in, bool const do_reorder,
      bool const use_aprx):
  functor(functor),
  ndim(mu_in.n_elem),
  n_integrands(functor.get_n_integrands()),
  use_aprx(use_aprx),
  infin(([&](){
    arma::ivec out(imem.get_mem(), ndim, false);
    pmvnorm::get_infin(out, lower_in, upper_in);
    return out;
  })()),
  indices(infin.end(), ndim, false) {
    /* checks */
    if(lower_in.n_elem > 1000 or lower_in.n_elem < 1)
      throw std::invalid_argument("cdf<T_Functor, out_type>: Either dimension zero or dimension greater than 1000");

#ifdef DO_CHECKS
    if(sigma_in.n_cols != static_cast<size_t>(ndim) or
         sigma_in.n_rows != static_cast<size_t>(ndim))
      throw std::invalid_argument("cdf::cdf: invalid 'sigma_in'");
    if(n_integrands <= 0L)
      throw std::invalid_argument("cdf::cdf: invalid 'n_integrands'");
    if(lower_in.n_elem !=  static_cast<size_t>(ndim))
      throw std::invalid_argument("cdf::cdf: invalid 'lower_in'");
    if(upper_in.n_elem !=  static_cast<size_t>(ndim))
      throw std::invalid_argument("cdf::cdf: invalid 'upper_in'");
#endif

    /* re-scale */
    double * cur_dtmp_mem = dtmp_mem;
    auto get_dmem = [&](int const n_ele) -> double * {
      double * out = cur_dtmp_mem;
      cur_dtmp_mem += n_ele;
      return out;
    };

    arma::vec sds(get_dmem(ndim), ndim, false),
    mu (get_dmem(ndim), ndim, false);
    for(int i = 0; i < ndim; ++i){
      sds[i] = std::sqrt(sigma_in.at(i, i));
      mu [i] = mu_in[i] / sds[i];
    }

    lower  = lower_in;
    lower /= sds;
    lower -= mu;

    upper  = upper_in;
    upper /= sds;
    upper -= mu;

    is_permutated = false;
    {
      int *idx = indices.begin();
      for(int i = 0; i < ndim; ++i, ++idx)
        *idx = i;
    }

    if(do_reorder and ndim > 1L){
      double * const y     = draw.begin(),
             * const A     = get_dmem(ndim),
             * const B     = get_dmem(ndim),
             * const DL    = sds.begin(),
             * const delta = mu.begin();
      sds.zeros();

      auto const correl = pmvnorm::get_cor_vec(sigma_in);
      int const pivot = 1L, doscale = 1L;
      int F_inform = 0L,
             nddim = ndim;
      std::fill(delta, delta + ndim, 0.);
      arma::ivec infi(itmp_mem, ndim, false);

      F77_CALL(mvsort)(
        &ndim, lower.memptr(), upper.memptr(), delta,
        correl.cor_vec.memptr(), infin.begin(), y, &pivot, &nddim, A, B,
        DL, sigma_chol.memptr(), infi.memptr(), &F_inform, indices.begin(),
        &doscale);

      if(F_inform != 0)
        throw std::runtime_error("cdf::cdf: error in mvsort");

      for(int i = 0; i < ndim; ++i){
        if(indices[i] != i){
          is_permutated = true;
          break;
        }
      }

      if(is_permutated){
        for(int i = 0; i < ndim; ++i){
          lower[i] = *(A + i);
          upper[i] = *(B + i);
          infin[i] = infi[i];
        }

        arma::mat sigma_permu(get_dmem(ndim * ndim), ndim, ndim, false);
        for(int j = 0; j < ndim; ++j)
          for(int i = 0; i < ndim; ++i)
            sigma_permu.at(i, j) = sigma_in.at(indices[i], indices[j]);
        functor.prep_permutated(sigma_permu, indices.begin());
        return;

      } else
        for(int i = 0; i < ndim; ++i){
          lower[i] = *(A + i);
          upper[i] = *(B + i);
        }

    } else if(!do_reorder and ndim > 1L) {
      arma::mat tmp(get_dmem(ndim * ndim), ndim, ndim, false);

      tmp = sigma_in;
      tmp.each_row() /= sds.t();
      tmp.each_col() /= sds;
      if(!arma::chol(tmp, tmp)) // TODO: memory allocation
        sigma_chol.fill(std::numeric_limits<double>::infinity());
      else
        copy_upper_tri(tmp, sigma_chol.memptr());

      if(ndim > 1L){
        /* rescale such that choleksy decomposition has ones in the diagonal */
        double * sc = sigma_chol.begin();
        for(int i = 0; i < ndim; ++i){
          double const scal = *(sc + i);
          lower[i] /= scal;
          upper[i] /= scal;
          double * const sc_end = sc + i + 1L;
          for(; sc != sc_end; ++sc)
            *sc /= scal;
        }
      }
    }
  }

  /**
   * Function to be called from mvkbrv.
   *
   * @param ndim_in Passed dimension of the integral.
   * @param unifs (0, 1) variables.
   * @param n_integrands_in Passed dimension of the integrand.
   * @param integrand_val Passed integrand approximation to be updated.
   */
  void operator()(
      int const *ndim_in, double const * unifs, int const *n_integrands_in,
      double * __restrict__ integrand_val) MDGC_NOEXCEPT {
#ifdef DO_CHECKS
    if(*ndim_in         != ndim)
      throw std::invalid_argument("cdf::eval_integrand: invalid 'ndim_in'");
    if(*n_integrands_in != n_integrands)
      throw std::invalid_argument("cdf::eval_integrand: invalid 'n_integrands_in'");
#endif

    double * const __restrict__ out = integrand_val,
           * const __restrict__ dr  = draw.begin();

    double w(1.);
    double const * __restrict__ sc   = sigma_chol.begin(),
                 * __restrict__ lw   = lower.begin(),
                 * __restrict__ up   = upper.begin(),
                 * __restrict__ unif = unifs;
    int const *infin_j = infin.begin();
    /* loop over variables and transform them to truncated normal
     * variables */
    for(int j = 0; j < ndim; ++j, ++sc, ++lw, ++up, ++infin_j, ++unif){
      double su(0.);
      double const *d = dr;
      for(int i = 0; i < j; ++i, sc++, d++)
        su += *sc * *d;

      auto pnorm_use = [&](double const x){
        return use_aprx ? pnorm_approx(x) : pnorm_std(x, 1L, 0L);
      };
      double lim_l(0.),
             lim_u(1.);
      if(*infin_j == 0L)
        lim_u = pnorm_use(*up - su);
      else if(*infin_j == 1L)
        lim_l = pnorm_use(*lw - su);
      else if(*infin_j == 2L){
        lim_l = pnorm_use(*lw - su);
        lim_u = pnorm_use(*up - su);

      }

      if(lim_l < lim_u){
        w *= lim_u - lim_l;
        if(needs_last_unif() or j + 1 < ndim){
          double const quant_val = lim_l + *unif * (lim_u - lim_l);
          *(dr + j) =
            use_aprx ?
            safe_qnorm_aprx(quant_val) :
            safe_qnorm_w   (quant_val);
        }

      } else {
        w = 0.;
        std::fill(dr + j, dr + ndim, 0.);
        break;

      }
    }

    /* evaluate the integrand and weigth the result. */
    functor(dr, out, indices.begin(), is_permutated);

    double * o = out;
    for(int i = 0; i < n_integrands; ++i, ++o)
      *o *= w;
  }

  /**
   * Integrates over (-Inf, 0)^[# of random effects]
   *
   * @functor T_Functor class to use.
   * @param mu_in Mean vector.
   * @param sigma_in Covariance matrix.
   * @param do_reorder true if the order of integrations may be reordered.
   */
  cdf(T_Functor &functor, arma::vec const &mu_in, arma::mat const &sigma_in,
      bool const do_reorder):
  cdf(
    functor,
    ([&](){
      arma::vec out(mu_in.n_elem);
      out.fill(-std::numeric_limits<double>::infinity());
      return out;
    })(), arma::vec(mu_in.n_elem, arma::fill::zeros), mu_in, sigma_in,
    do_reorder, false) { }

  /**
   * Performs the approximation.
   *
   * @param maxvls Maximum number of function evaluations allowed.
   * @param abs_eps Required absolute accuracy.
   * @param rel_eps Required relative accuracy.
   * @param minvls Minimum number of samples.
   */
  out_type approximate
  (int const maxvls, double const abs_eps, double const rel_eps,
   int const minvls = 0L){
#ifdef DO_CHECKS
    if(abs_eps <= 0 and rel_eps <= 0)
      throw std::invalid_argument("cdf::approximate: no valid convergence threshold");
    if(maxvls <= 0L)
      throw std::invalid_argument("cdf::approximate: invalid 'maxvls'");
#endif

    // setup
    // needs to have at least n_integrands memory to use. The functor class
    // needs to provide this
    double * const int_apprx = functor.get_wk_mem();

    auto sampler = parallelrng::get_unif_drawer();

    if(ndim == 1L){
      /* handle the one-dimensional case as a special case */
      functor.univariate(int_apprx, lower[0], upper[0]);
      indices[0] = 0;

      return functor.get_output(int_apprx, 0, 0, 0,
                                indices.begin());

    } else if(std::isinf(*sigma_chol.begin()))
      throw std::runtime_error("std::isinf(*sigma_chol.begin())");

    /* perform the approximation */
    auto res = rand_Korobov<cdf<T_Functor> >::comp(
      *this, ndim, minvls, maxvls, n_integrands, abs_eps, rel_eps,
      int_apprx, sampler);
    return functor.get_output(int_apprx, res.minvls, res.inform,
                              res.abserr, indices.begin());
  }
};

template<class T_Functor, class out_type>
cache_mem<int   > cdf<T_Functor, out_type>::imem;
template<class T_Functor, class out_type>
cache_mem<double> cdf<T_Functor, out_type>::dmem;

/**
 * func classes used as template argument for cdf used to approximate the
 * likelihood. */
class likelihood {
  static cache_mem<double> dmen;

public:
  static void alloc_mem
  (unsigned const max_dim, unsigned const max_threads){
    rand_Korobov<cdf<likelihood> >::alloc_mem(
        max_dim, get_n_integrands(), max_threads);
    dmen.set_n_mem(1, max_threads);
    cdf<likelihood>::alloc_mem(max_dim, max_threads);
  }
  double * get_wk_mem(){
    return dmen.get_mem();
  }
  constexpr static int get_n_integrands() {
    return 1L;
  }

  inline void operator()
  (double const *, double * out, int const *, bool const)
  MDGC_NOEXCEPT {
#ifdef DO_CHECKS
    if(!out)
      throw std::invalid_argument("likelihood::operator(): invalid out");
#endif
    *out = 1;
  }

  constexpr static bool needs_last_unif() {
    return false;
  }

  inline static void univariate(double * out,
                                double const lw, double const ub) {
    double const p_ub = std::isinf(ub) ? 1 : pnorm_std(ub, 1L, 0L),
      p_lb = std::isinf(lw) ? 0 : pnorm_std(lw, 1L, 0L);
    *out = p_ub - p_lb;
  }

  struct out_type {
    /**
     * minvls Actual number of function evaluations used.
     * inform INFORM = 0 for normal exit, when
     *             ABSERR <= MAX(ABSEPS, RELEPS*||finest||)
     *          and
     *             INTVLS <= MAXCLS.
     *        INFORM = 1 If MAXVLS was too small to obtain the required
     *        accuracy. In this case a value finest is returned with
     *        estimated absolute accuracy ABSERR. */
    int minvls, inform;
    /// maximum norm of estimated absolute accuracy of finest
    double abserr;
    /// likelihood approximation
    double likelihood;
  };
  out_type get_output(double const * res, int const minvls,
                      int const inform, double const abserr,
                      int const *){
    return out_type { minvls, inform, abserr, *res };
  }
  inline void prep_permutated(arma::mat const&, int const *) { }
};

/**
 * func classes used as template argument for cdf used to approximates both
 * the probability and the derivatives w.r.t. the mean and the covariance
 * matrix. */
class deriv {
  /// working memory
  static cache_mem<double> dmem;

  int const ndim;
  /**
   * points to the upper triangular part of the inverse of the Cholesky
   * decomposition.
   */
  double * sigma_chol_inv;
  /// points to the upper triangular part of the inverse.
  double * sig_inv;
  /// working memory to be used by cdf
  double * cdf_mem;

  inline void deriv_integrand_inner_loop
    (double * __restrict__ o, double const * __restrict__ lhs,
     unsigned const c) noexcept {
    double const * const end = lhs + c + 1;
    double const mult = *(lhs + c);
    for(; lhs != end; ++o, ++lhs)
      *o = mult * * lhs;
  }

public:
  /**
   * must be be called before calling other member functions to allocate
   * working memory.
   */
  static void alloc_mem
  (unsigned const max_dim, unsigned const max_threads){
    rand_Korobov<cdf<deriv> >::alloc_mem(
        max_dim, get_n_integrands(max_dim), max_threads);
    dmem.set_n_mem(
      2 * max_dim * max_dim + max_dim * (max_dim + 1) +
        get_n_integrands(max_dim),
      max_threads);
    cdf<deriv>::alloc_mem(max_dim, max_threads);
  }

  /**
   Note that alloc_mem must have been called before calling this method.
   */
  deriv(arma::vec const &mu, arma::mat const &sig):
  ndim(mu.n_elem) {
    // create the objects we need
    arma::mat t1(dmem.get_mem(), ndim, ndim, false),
              t2(t1.end()      , ndim, ndim, false);
    if(!arma::chol(t1, sig))
      throw std::runtime_error("deriv::deriv: chol failed");
    if(!arma::inv(t2, t1))
      throw std::runtime_error("deriv::deriv: inv failed");
    sigma_chol_inv = t2.end();
    copy_upper_tri(t2, sigma_chol_inv);

    if(!arma::inv_sympd(t1, sig))
      throw std::runtime_error("deriv::deriv: inv_sympd failed");
    sig_inv = sigma_chol_inv + (ndim * (ndim + 1)) / 2;
    copy_upper_tri(t1, sig_inv);

    cdf_mem = sig_inv + (ndim * (ndim + 1)) / 2;
  }

  void prep_permutated(arma::mat const &sigma_permu, int const *) {
    // need to re-compute the inverse of the Cholesky decomposition
    arma::mat t1(dmem.get_mem(), ndim, ndim, false),
              t2(t1.end()      , ndim, ndim, false);
    if(!arma::chol(t1, sigma_permu))
      throw std::runtime_error("deriv::prep_permutated: chol failed");
    if(!arma::inv(t2, t1))
      throw std::runtime_error("deriv::prep_permutated: inv failed");
    copy_upper_tri(t2, sigma_chol_inv);

    if(!arma::inv_sympd(t1, sigma_permu))
      throw std::runtime_error("deriv::prep_permutated: inv_sympd failed");
    copy_upper_tri(t1, sig_inv);
  }

  inline static int get_n_integrands(int const x) noexcept {
    return 1 + x + (x * (x + 1)) / 2L;
  }

  inline int get_n_integrands() noexcept{
    return get_n_integrands(ndim);
  }

  double * get_wk_mem(){
    return cdf_mem;
  }

  constexpr static bool needs_last_unif() {
    return true;
  }

  inline void operator()
  (double const * __restrict__ draw, double * __restrict__ out,
   int const *indices, bool const is_permutated) {
    arma::uword const p = ndim;

    int const n_elem = 1L + p + (p * (p + 1L)) / 2L;
    *out = 1.;
    std::fill(out + 1L, out + n_elem, 0.);


    double * const mean_part_begin = out + 1L;
    /* Multiplying by the inverse matrix is fast but not smart
     * numerically. */
    double const * sigma_chol_inv_ck = sigma_chol_inv;
    for(unsigned c = 0; c < p; ++c){
      double const mult = *(draw + c),
        * const end = mean_part_begin + c + 1L;
      for(double *rhs = mean_part_begin; rhs != end;
          ++rhs, ++sigma_chol_inv_ck)
        *rhs += mult * *sigma_chol_inv_ck;
    }

    double * o = out + 1L + p;
    for(unsigned c = 0; c < p; c++){
      deriv_integrand_inner_loop(o, mean_part_begin, c);
      o += c + 1L;
    }
  }

  inline void univariate(double * out, double const lw, double const ub) {
    constexpr double const sqrt_2_pi_inv = 0.398942280401433;
    auto dnrm = [&](double const x){
      return std::exp(-x * x / 2.) * sqrt_2_pi_inv;
    };
    bool const f_ub = std::isinf(ub),
               f_lb = std::isinf(lw);

    double const p_ub = f_ub ? 1 : pnorm_std(ub, 1L, 0L),
                 p_lb = f_lb ? 0 : pnorm_std(lw, 1L, 0L),
                 d_ub = f_ub ? 0 : dnrm(ub),
                 d_lb = f_lb ? 0 : dnrm(lw),
              d_ub_ub = f_ub ? 0 : ub * d_ub,
              d_lb_lb = f_lb ? 0 : lw * d_lb,
               sd_inv = *sigma_chol_inv;

    out[0L] = p_ub - p_lb;
    out[1L] = -(d_ub - d_lb) * sd_inv;
    out[2L] = -(d_ub_ub - d_lb_lb) * sd_inv * sd_inv + out[0L] * *sig_inv;
  }

  struct out_type {
    /**
     * minvls Actual number of function evaluations used.
     * inform INFORM = 0 for normal exit, when
     *             ABSERR <= MAX(ABSEPS, RELEPS*||finest||)
     *          and
     *             INTVLS <= MAXCLS.
     *        INFORM = 1 If MAXVLS was too small to obtain the required
     *        accuracy. In this case a value finest is returned with
     *        estimated absolute accuracy ABSERR. */
    int minvls, inform;
    /// maximum norm of estimated absolute accuracy of finest
    double abserr;
    /// likelihood approximation
    double likelihood;
    /// the derivative approximation
    arma::vec derivs;
  };

  out_type get_output(double * const res, int const minvls,
                      int const inform, double const abserr,
                      int const *indices){
    // post process
    double const likelihood = *res;
    double *o = res + 1 + ndim;
    double const * sig_inv_ij = sig_inv;
    for(int c = 0; c < ndim; c++){
      double const * const end = sig_inv_ij + c + 1L;
      for(; sig_inv_ij != end; ++sig_inv_ij, ++o){
        *o -= likelihood * *sig_inv_ij;
        *o /= 2.;
      }
    }

    // gather result
    out_type out;
    out.minvls = minvls;
    out.inform = inform;
    out.abserr = abserr;
    out.likelihood = likelihood;

    arma::vec &derivs = out.derivs;
    int const dim_cov = (ndim * (ndim + 1L)) / 2L;
    derivs.set_size(ndim + dim_cov);

    double * const dmu  = res + 1,
           * const dcov = dmu + ndim,
           * const cmu  = derivs.begin(),
           * const ccov = cmu + ndim;

    /* maps to the indices of the upper triangle matrix */
    auto tri_map = [&](int const r, int const c){
      return r + (c * (c + 1L)) / 2L;
    };

    /* permutate back and return */
    for(int c = 0; c < ndim; ++c){
      int const org_c = indices[c];
      cmu[org_c] = dmu[c];

      for(int r = 0; r <= c; ++r){
        double const val = dcov[tri_map(r, c)];
        int const org_r = indices[r];
        if(org_c >= org_r)
          ccov[tri_map(org_r, org_c)] = val;
        else
          ccov[tri_map(org_c, org_r)] = val;
      }
    }

    return out;
  }
};

class imputation {
public:
  class type_base {
  public:
    virtual int n_ele() const noexcept = 0;
    virtual void set_val(double const, double *& __restrict__)
      const noexcept = 0;
    virtual ~type_base() = default;
  };

  class known final : public type_base {
  public:
    inline int n_ele() const noexcept {
      return 0L;
    };
    inline void set_val(double const, double *& __restrict__)
      const noexcept { };
  };

  class contin final : public type_base {
  public:
    inline int n_ele() const noexcept {
      return 1L;
    }

    inline void set_val(double const v, double *& __restrict__ res)
    const noexcept {
      *res++ = v;
    }
  };

  class binary final : public type_base {
  public:
    double const border;

    binary(double const border): border(border) { }

    inline int n_ele() const noexcept {
      return 2L;
    }

    inline void set_val(double const v, double *& __restrict__ res)
      const noexcept {
      if(v < border){
        *res++ = 1.;
        *res++ = 0.;

      } else {
        *res++ = 0.;
        *res++ = 1.;

      }
    }
  };

  class ordinal final : public type_base {
  public:
    int const n_bs;
    std::unique_ptr<double[]> const borders;

    ordinal(double const *br, int const n_ele):
    n_bs(n_ele - 2L),
    borders(([&](){
#ifdef DO_CHECKS
      if(n_ele < 3L)
        throw std::invalid_argument("ordinal: n_ele less than three");
#endif
      std::unique_ptr<double[]> out(new double[n_bs]);
      for(int i = 0; i < n_bs; ++i)
        *(out.get() + i) = *(br + i + 1L);
      return out;
    })()) { }

    ordinal(ordinal const &o):
    n_bs(o.n_bs),
    borders(([&](){
      std::unique_ptr<double[]> out(new double[n_bs]);
      for(int i = 0; i < n_bs; ++i)
        *(out.get() + i) = *(o.borders.get() + i);
      return out;
    })()) { }

    inline int n_ele() const noexcept {
      return n_bs + 1L;
    }

    inline void set_val(double const v, double *& __restrict__ res)
    const noexcept {
      int i = 0L;
      double const *b = borders.get();
      for(; i < n_bs; ++i, ++b)
        if(v < *b){
          *res++ = 1.;
          break;
        } else
          *res++ = 0.;

      if(i == n_bs)
        *res++ = 1.;
      else {
        ++i;
        for(; i <= n_bs; ++i)
          *res++ = 0.;
      }
    }
  };

private:
  std::vector<type_base const*> cur_list;

  /// working memory
  static cache_mem<double> dmem;
  int const n_integrands_val = get_n_integrands(cur_list),
            ndim;
  double * __restrict__ const mu_vec   = dmem.get_mem(),
         * __restrict__ const sig_chol = mu_vec + ndim,
         * __restrict__ const cdf_mem  = sig_chol + (ndim * (ndim + 1)) / 2,
         * __restrict__ const xtr_mem  = cdf_mem + n_integrands_val;

public:
  /**
   * must be be called before calling other member functions to allocate
   * working memory.
   */
  static void alloc_mem(std::vector<type_base const*> const cur_list,
                        int const max_dim, int const max_threads){
    int const n_ints = get_n_integrands(cur_list);
    rand_Korobov<cdf<imputation> >::alloc_mem(max_dim, n_ints, max_threads);
    cdf<imputation>::alloc_mem(max_dim, max_threads);
    dmem.set_n_mem(
      (max_dim * (max_dim + 5)) / 2 + max_dim * max_dim + n_ints,
      max_threads);
  }

  imputation(std::vector<type_base const *> const &cur_list,
             arma::vec const &mu, arma::vec const &Sig):
    cur_list(cur_list), ndim(mu.n_elem) {
    /* set the mean vector and the Cholesky decomposition of the random
     * of the random effects. */
    std::copy(mu.begin(), mu.end(), mu_vec);

    arma::mat tmp_mat(xtr_mem, ndim, ndim, false, true);
    if(!arma::chol(tmp_mat, Sig))
      throw std::runtime_error("imputation::imputation: chol failed");
    copy_upper_tri(tmp_mat, sig_chol);
  }

  void prep_permutated
    (arma::mat const &sigma_permu, int const *indices){
    // permute the mean vector, the covariance matrix, and the cur_list
    arma::mat tmp_mat(xtr_mem, ndim, ndim, false, true);
    if(!arma::chol(tmp_mat, sigma_permu))
      throw std::runtime_error("imputation::prep_permutated: chol failed");
    copy_upper_tri(tmp_mat, sig_chol);

    for(int i = 0; i < ndim; ++i)
      // we do not use this memory now as of this writting
      cdf_mem[i] = mu_vec[indices[i]];
    std::copy(cdf_mem, cdf_mem + ndim, mu_vec);

    // permutate the type_base list. TODO: memory allocation
    int const n_ele = cur_list.size();
    std::vector<type_base const *> new_list;
    new_list.reserve(n_ele);
    for(int i = 0; i < n_ele; ++i)
      new_list.emplace_back(cur_list[indices[i]]);
    cur_list = new_list;
  }

  static inline int
  get_n_integrands(std::vector<type_base const *> const &ls) noexcept {
    int out(1L);
    for(auto &x : ls)
      out += x->n_ele();
    return out;
  }

  inline int get_n_integrands() noexcept {
    return n_integrands_val;
  }

  double * get_wk_mem(){
    return cdf_mem;
  }

  constexpr static bool needs_last_unif() {
    return true;
  }

  inline void operator()
  (double const * __restrict__ draw, double * __restrict__ out,
   int const *indices, bool const is_permutated){
    double *scale_draw = xtr_mem;
    std::copy(mu_vec, mu_vec + ndim, scale_draw);

    {
      // re-scale
      double const *sig_chol_rp = sig_chol;
      double * scale_draw_c = scale_draw;
      for(int c = 0; c < ndim; ++c, ++scale_draw_c){
        double const * draw_r = draw;
        for(int r = 0; r <= c; ++r, ++sig_chol_rp, ++draw_r)
          *scale_draw_c += *sig_chol_rp * *draw_r;
      }
    }

    double * o = out;
    *o++ = 1.;
    double const *scale_draw_i = scale_draw;
    for(int i = 0; i < static_cast<int>(cur_list.size());
        ++i, ++scale_draw_i)
      cur_list[i]->set_val(*scale_draw_i, o);
  }

  inline void univariate(double * out, double const lw, double const ub) {
    throw std::runtime_error("imputation::univariate: not implemented");
  }

  struct out_type {
    /**
     * minvls Actual number of function evaluations used.
     * inform INFORM = 0 for normal exit, when
     *             ABSERR <= MAX(ABSEPS, RELEPS*||finest||)
     *          and
     *             INTVLS <= MAXCLS.
     *        INFORM = 1 If MAXVLS was too small to obtain the required
     *        accuracy. In this case a value finest is returned with
     *        estimated absolute accuracy ABSERR. */
    int minvls, inform;
    /// maximum norm of estimated absolute accuracy of finest
    double abserr;
    /// likelihood approximation
    double likelihood;
    /// the imputed values
    arma::vec imputations;
  };

  out_type get_output(double * res, int const minvls,
                      int const inform, double const abserr,
                      int const *indices){
    // permutate back the result and return
    out_type out;
    out.minvls = minvls;
    out.inform = inform;
    out.abserr = abserr;
    out.likelihood = *res;

    arma::vec &imputations = out.imputations;
    imputations.set_size(get_n_integrands() - 1);
    arma::uvec offset(ndim);

    /* we make two passes. One to figure out where to write the ith type
     * to and one to do the writting */
    for(int c = 0; c < ndim; ++c)
      offset[indices[c]] = cur_list[c]->n_ele();

    {
      int cur = 1L; // one as the first element is the CDF approximation
      for(int i = 0L; i < ndim; ++i){
        cur += offset[i];
        offset[i] = cur - offset[i];
      }
    }

    // do the copying
    {
      double const *v = res + 1L;
      for(int c = 0; c < ndim; ++c){
        int const org_c = indices[c];
        int const n_ele = cur_list[c]->n_ele();
        double *r = imputations.begin() + offset[org_c] - 1;
        for(int i = 0; i < n_ele; ++i, ++r, ++v)
          *r = *v;
      }
    }

    return out;
  }
};
} // namespace restrictcdf

#endif
