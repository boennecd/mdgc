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

inline double safe_qnorm_w(double const x) MDGC_NOEXCEPT {
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

inline double safe_qnorm_aprx(double const x) MDGC_NOEXCEPT {
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
  (arma::mat const &X, double * MDGC_RESTRICT x) MDGC_NOEXCEPT {
  int const p = X.n_cols;
  for(int c = 0; c < p; c++)
    for(int r = 0; r <= c; r++, x++)
      *x = X.at(r, c);
}

inline void copy_lower_tri
  (arma::mat const &X, double * MDGC_RESTRICT x) MDGC_NOEXCEPT {
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

  double * const lower = dmem.get_mem(),
         * const upper = lower + ndim,
    * const sigma_chol = upper + ndim,
    * const draw       = sigma_chol + (ndim * (ndim + 1L)) / 2L,
    * const dtmp_mem   = draw + ndim * n_qmc_seqs();
  // memory that can be used
  int * const itmp_mem = indices.end();

public:
  /**
   * must be called perior to calling the constructor or any member
   * functions.
   */
  static void alloc_mem(int const max_ndim, int const max_threads) {
    int const n_up_tri = (max_ndim * (max_ndim + 1)) / 2;
    imem.set_n_mem(3 * max_ndim                                 ,
                   max_threads);
    dmem.set_n_mem(
      (6 + n_qmc_seqs()) * max_ndim + n_up_tri + max_ndim * max_ndim +
        2 * n_qmc_seqs(),
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

    double * sds = get_dmem(ndim);
    for(int i = 0; i < ndim; ++i){
      sds[i] = std::sqrt(sigma_in.at(i, i));
      lower[i] = (lower_in[i] - mu_in[i]) / sds[i];
      upper[i] = (upper_in[i] - mu_in[i]) / sds[i];
    }

    is_permutated = false;
    {
      int *idx = indices.begin();
      for(int i = 0; i < ndim; ++i, ++idx)
        *idx = i;
    }

    if(do_reorder and ndim > 1L){
      double * const y     = draw,
             * const A     = get_dmem(ndim),
             * const B     = get_dmem(ndim),
             * const DL    = sds,
             * const delta = get_dmem(ndim);
      std::fill(sds, sds + ndim, 0.);

      auto const correl = pmvnorm::get_cor_vec(sigma_in);
      int const pivot = 1L, doscale = 1L;
      int F_inform = 0L,
             nddim = ndim;
      std::fill(delta, delta + ndim, 0.);
      arma::ivec infi(itmp_mem, ndim, false);

      F77_CALL(mvsort)(
        &ndim, lower, upper, delta,
        correl.cor_vec.memptr(), infin.begin(), y, &pivot, &nddim, A, B,
        DL, sigma_chol, infi.memptr(), &F_inform, indices.begin(),
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
          lower[i] = A[i];
          upper[i] = B[i];
          infin[i] = infi[i];
        }

        arma::mat sigma_permu(get_dmem(ndim * ndim), ndim, ndim, false);
        for(int j = 0; j < ndim; ++j)
          for(int i = 0; i < ndim; ++i)
            sigma_permu.at(i, j) = sigma_in.at(indices[i], indices[j]);
        functor.prep_sim(sigma_permu, indices.begin(), true);
        return;

      } else
        for(int i = 0; i < ndim; ++i){
          lower[i] = A[i];
          upper[i] = B[i];
        }

    } else if(!do_reorder and ndim > 1L) {
      arma::mat tmp(get_dmem(ndim * ndim), ndim, ndim, false);

      tmp = sigma_in;
      for(int i = 0; i < ndim; ++i)
        for(int j = 0; j <ndim; ++j)
          tmp.at(i, j) /= sds[i] * sds[j];

      if(!arma::chol(tmp, tmp)) // TODO: memory allocation
        std::fill(sigma_chol, sigma_chol + (ndim * (ndim + 1L)) / 2L,
                  std::numeric_limits<double>::infinity());
      else
        copy_upper_tri(tmp, sigma_chol);

      if(ndim > 1L){
        /* re-scale such that Cholesky decomposition has ones in the diagonal */
        double * sc = sigma_chol;
        for(int i = 0; i < ndim; ++i){
          double const scal = sc[i];
          lower[i] /= scal;
          upper[i] /= scal;
          double * const sc_end = sc + i + 1L;
          for(; sc != sc_end; ++sc)
            *sc /= scal;
        }
      }
    } else
      *sigma_chol = 1.;

    functor.prep_sim(sigma_in, indices.begin(), false);
  }

  /**
   * Function to be called from mvkbrv.
   *
   * @param ndim_in Passed dimension of the integral.
   * @param unifs (0, 1) variables.
   * @param n_integrands_in Passed dimension of the integrand.
   * @param integrand_val Passed integrand approximation to be updated.
   * @param n_draws number of draws
   */
  void operator()(
      int const *ndim_in, double const * unifs, int const *n_integrands_in,
      double * MDGC_RESTRICT integrand_val, int const n_draws) MDGC_NOEXCEPT {
#ifdef DO_CHECKS
    if(*ndim_in         != ndim)
      throw std::invalid_argument("cdf::eval_integrand: invalid 'ndim_in'");
    if(*n_integrands_in != n_integrands)
      throw std::invalid_argument("cdf::eval_integrand: invalid 'n_integrands_in'");
#endif

    double * const MDGC_RESTRICT out = integrand_val,
           * const MDGC_RESTRICT dr  = draw,
           * const MDGC_RESTRICT su  = dr + n_draws * ndim,
           * const MDGC_RESTRICT w   = su + n_draws;

    std::fill(w, w + n_draws, 1);
    double const * MDGC_RESTRICT sc   = sigma_chol,
                 * MDGC_RESTRICT lw   = lower,
                 * MDGC_RESTRICT up   = upper;
    int const *infin_j = infin.begin();
    /* loop over variables and transform them to truncated normal
     * variables */
    for(int j = 0; j < ndim; ++j, ++sc, ++lw, ++up, ++infin_j){
      std::fill(su, su + n_draws, 0);
      {
        double const *d = dr;
        for(int i = 0; i < j; ++i, sc++, d += n_draws)
          for(int k = 0; k < n_draws; ++k)
            su[k] += *sc * d[k];
      }

      auto pnorm_use = [&](double const x){
        return use_aprx ? pnorm_approx(x) : pnorm_std(x, 1L, 0L);
      };

      int const offset = j * n_draws;
      for(int k = 0; k < n_draws; ++k){
        double lim_l(0.),
               lim_u(1.);
        if(*infin_j == 0L)
          lim_u = pnorm_use(*up - su[k]);
        else if(*infin_j == 1L)
          lim_l = pnorm_use(*lw - su[k]);
        else if(*infin_j == 2L){
          lim_l = pnorm_use(*lw - su[k]);
          lim_u = pnorm_use(*up - su[k]);

        }

        if(lim_l < lim_u){
          w[k] *= lim_u - lim_l;
          if(needs_last_unif() or j + 1 < ndim){
            double const quant_val =
              lim_l + unifs[k * ndim + j] * (lim_u - lim_l);
            dr[offset + k] =
              use_aprx ?
              safe_qnorm_aprx(quant_val) :
              safe_qnorm_w   (quant_val);
          }

        } else {
          w[k] = 0.;

        }
      }
    }

    /* evaluate the integrand and weight the result. */
    functor(dr, out, indices.begin(), is_permutated,
            static_cast<double const*>(w), n_draws);
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

    } else if(std::isinf(*sigma_chol))
      throw std::runtime_error("std::isinf(*sigma_chol)");

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
  (double const *, double * out, int const *, bool const,
   double const *w, int const n_draws)
  MDGC_NOEXCEPT {
#ifdef DO_CHECKS
    if(!out)
      throw std::invalid_argument("likelihood::operator(): invalid out");
#endif
    *out = 0;
    for(int k = 0; k < n_draws; ++k)
      *out += w[k];
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
  inline void prep_sim(arma::mat const&, int const *, bool const) { }
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
   * points to upper triangular matrix of the Cholesky decomposition as a full
   * ndim x ndim matrix.
   */
  double * const sigma_chol = dmem.get_mem();
  /// points to the upper triangular part of the inverse.
  double * const sig_inv = sigma_chol + ndim * ndim;
  /// working memory to be used by cdf
  double * const cdf_mem = sig_inv + (ndim * (ndim + 1)) / 2;
  /// working memory for operator()
  double * const internal_mem = cdf_mem + get_n_integrands(ndim);

  inline void deriv_integrand_inner_loop
    (double * MDGC_RESTRICT o, double const * MDGC_RESTRICT scaled_samp,
     int const c, int const n_draws) MDGC_NOEXCEPT {

    double const * mult = scaled_samp + c * n_draws;
    for(int i = 0; i <= c; ++i, scaled_samp += n_draws)
      for(int k = 0; k < n_draws; ++k)
        o[i] += scaled_samp[k] * mult[k];
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
      max_dim * max_dim + (max_dim * (max_dim + 1)) / 2 +
        get_n_integrands(max_dim) +
        n_qmc_seqs() * max_dim + n_qmc_seqs() + 2 * max_dim * max_dim,
      max_threads);
    cdf<deriv>::alloc_mem(max_dim, max_threads);
  }

  /**
   Note that alloc_mem must have been called before calling this method.
   */
  deriv(arma::vec const &mu, arma::mat const &sig):
  ndim(mu.n_elem) { }

  void prep_sim(arma::mat const &sigma_permu, int const *,
                bool const is_permuted) {
    // need to re-compute the inverse of the Cholesky decomposition
    arma::mat t1(internal_mem, ndim, ndim, false);
    if(!arma::chol(t1, sigma_permu, "upper"))
      throw std::runtime_error("deriv::prep_sim: chol failed");
    std::copy(t1.begin(), t1.end(), sigma_chol);

    if(!arma::inv_sympd(t1, sigma_permu))
      throw std::runtime_error("deriv::prep_sim: inv_sympd failed");
    copy_upper_tri(t1, sig_inv);
  }

  inline static int get_n_integrands(int const x) MDGC_NOEXCEPT {
    return 1 + x + (x * (x + 1)) / 2L;
  }

  inline int get_n_integrands() MDGC_NOEXCEPT {
    return get_n_integrands(ndim);
  }

  double * get_wk_mem(){
    return cdf_mem;
  }

  constexpr static bool needs_last_unif() {
    return true;
  }

  inline void operator()
  (double const * MDGC_RESTRICT draw, double * MDGC_RESTRICT out,
   int const *indices, bool const is_permutated, double const *w,
   int const n_draws) {
    arma::uword const p = ndim;

    *out = 0;
    for(int k = 0; k < n_draws; ++k)
      *out += w[k];

    int const n_elem = 1L + p + (p * (p + 1L)) / 2L;
    double * MDGC_RESTRICT const scaled_samp_mult = internal_mem,
           * MDGC_RESTRICT const wk_mem           = scaled_samp_mult + n_draws * ndim;

    std::fill(out + 1L, out + n_elem, 0.);
    // compute the derivative w.r.t. the mean
    {
      double const * d = draw;
      for(int i = 0; i < ndim; ++i, d += n_draws)
        for(int k = 0; k < n_draws; ++k)
          out[1 + i] += w[k] * d[k];
    }

    // compute the square root of weights times sc_samp
    for(int k = 0; k < n_draws; ++k)
      wk_mem[k] = std::sqrt(w[k]);
    {
        double * MDGC_RESTRICT to = scaled_samp_mult;
        double const * from = draw;
        for(int i = 0; i < ndim; ++i, to += n_draws, from += n_draws)
          for(int k = 0; k < n_draws; ++k)
            to[k] = from[k] * wk_mem[k];
    }

    double * o = out + 1L + p;
    for(unsigned c = 0; c < p; c++){
      deriv_integrand_inner_loop(o, scaled_samp_mult, c, n_draws);
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
               sd_inv = 1 / *sigma_chol;

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
    // multiply by the inverse of the Cholesky factorization.
    if(ndim > 1){
      arma::mat c_mat(sigma_chol, ndim, ndim, false);

      // start with the mean part
      {
        double * MDGC_RESTRICT const d_mu_tmp = internal_mem,
               * MDGC_RESTRICT const d_mu     = res + 1;
        std::copy(d_mu, d_mu + ndim, d_mu_tmp);

        arma::mat lhs(d_mu      , ndim,    1, false),
                  rhs(d_mu_tmp  , ndim,    1, false);

        arma::solve(lhs, arma::trimatu(c_mat), rhs);
      }

      // then the covariance matrix part
      {
        arma::mat d_sig_res(internal_mem   , ndim, ndim, false),
                  d_sig_tmp(d_sig_res.end(), ndim, ndim, false);
        {
          double const * MDGC_RESTRICT d_sig  = res + 1 + ndim;
          for(int j = 0; j < ndim; ++j)
            for(int i = 0; i <= j; ++i){
              d_sig_res.at(i, j) = *d_sig;
              d_sig_res.at(j, i) = *d_sig++;
            }
        }

        arma::solve(d_sig_tmp, arma::trimatu(c_mat), d_sig_res);
        arma::inplace_trans(d_sig_tmp);
        arma::solve(d_sig_res, arma::trimatu(c_mat), d_sig_tmp);

        // copy the result
        double * MDGC_RESTRICT d_sig  = res + 1 + ndim;
        for(int j = 0; j < ndim; ++j)
          for(int i = 0; i <= j; ++i)
            *d_sig++ = d_sig_res.at(i, j);
      }
    }

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

    /* permute back and return */
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
    virtual int n_latent() const MDGC_NOEXCEPT = 0;
    virtual int n_ele() const MDGC_NOEXCEPT = 0;
    virtual void set_val(double const *, double * MDGC_RESTRICT,
                         double const)
      const MDGC_NOEXCEPT = 0;
    virtual ~type_base() = default;
  };

  class known final : public type_base {
    int const n_latent_val;
  public:
    inline int n_latent() const MDGC_NOEXCEPT {
      return n_latent_val;
    }
    inline int n_ele() const MDGC_NOEXCEPT {
      return 0L;
    };
    inline void set_val(double const *, double * MDGC_RESTRICT,
                        double const)
      const MDGC_NOEXCEPT { };

    known(int const n_latent_val = 1): n_latent_val(n_latent_val) { }
  };

  class contin final : public type_base {
  public:
    inline int n_latent() const MDGC_NOEXCEPT {
      return 1L;
    }

    inline int n_ele() const MDGC_NOEXCEPT {
      return 1L;
    }

    inline void set_val(double const *v, double * MDGC_RESTRICT res,
                        double const weight)
    const MDGC_NOEXCEPT {
      *res += weight * *v;
    }
  };

  class binary final : public type_base {
  public:
    double const border;

    binary(double const border): border(border) { }

    inline int n_latent() const MDGC_NOEXCEPT {
      return 1L;
    }

    inline int n_ele() const MDGC_NOEXCEPT {
      return 2L;
    }

    inline void set_val(double const *v, double * MDGC_RESTRICT res,
                        double const weight)
      const MDGC_NOEXCEPT {
      if(*v < border)
        res[0] += weight;
      else
        res[1] += weight;
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
      std::copy(br + 1L, br + 1L + n_bs, out.get());
      return out;
    })()) { }

    ordinal(ordinal const &o):
    n_bs(o.n_bs),
    borders(([&](){
      std::unique_ptr<double[]> out(new double[o.n_bs]);
      std::copy(o.borders.get(), o.borders.get() + o.n_bs, out.get());
      return out;
    })()) { }

    inline int n_latent() const MDGC_NOEXCEPT {
      return 1L;
    }

    inline int n_ele() const MDGC_NOEXCEPT {
      return n_bs + 1L;
    }

    inline void set_val(double const *v, double * MDGC_RESTRICT res,
                        double const weight)
    const MDGC_NOEXCEPT {
      int i = 0L;
      double const *b = borders.get();
      for(; i < n_bs; ++i, ++b)
        if(*v < *b){
          *res += weight;
          break;
        } else
          ++res;

      if(i == n_bs)
        *res += weight;
    }
  };

  class multinomial final : public type_base {
  public:
    int const n_lvls;

    multinomial(int const n_lvls): n_lvls(n_lvls) { }

    int n_latent() const MDGC_NOEXCEPT {
      return n_lvls;
    }
    int n_ele() const MDGC_NOEXCEPT {
      return n_lvls;
    };

    void set_val
      (double const *v, double * MDGC_RESTRICT res,
       double const weight) const MDGC_NOEXCEPT {
      int max_lvl(0);
      double max_val = *v;

      for(int i = 1; i < n_lvls; ++i)
        if(v[i] > max_val){
          max_lvl = i;
          max_val = v[i];
        }

      res[max_lvl] += weight;
    }
  };

private:
  std::vector<type_base const*> cur_list;

  /// working memory
  static cache_mem<double> dmem;
  int const n_integrands_val = get_n_integrands(cur_list),
            ndim;
  double * MDGC_RESTRICT const mu_vec   = dmem.get_mem(),
         * MDGC_RESTRICT const sig_chol = mu_vec + ndim,
         * MDGC_RESTRICT const cdf_mem  = sig_chol + (ndim * (ndim + 1)) / 2,
         * MDGC_RESTRICT const xtr_mem  = cdf_mem + n_integrands_val;

public:
  /**
   * must be be called before calling other member functions to allocate
   * working memory.
   */
  static void alloc_mem(std::vector<type_base const*> const cur_list,
                        int const max_threads){
    int const n_ints = get_n_integrands(cur_list);
    int n_latent(0);
    for(auto &x : cur_list)
      n_latent += x->n_latent();
    rand_Korobov<cdf<imputation> >::alloc_mem(n_latent, n_ints, max_threads);
    cdf<imputation>::alloc_mem(n_latent, max_threads);
    dmem.set_n_mem(
      // mu_vec and sig_chol, and an extra n_latent vector
      (n_latent * (n_latent + 5)) / 2 +
        // a full n_latent x n_latent matrix and one n_latent vector
        n_latent * (n_latent + 1) +
        // memory for cdf_mem
        n_ints +
        // two extra matrices of size n_latent * ndraws
        2 * n_latent * n_qmc_seqs(),
      max_threads);
  }

  imputation(std::vector<type_base const *> const &cur_list,
             arma::vec const &mu, arma::mat const &Sig):
    cur_list(cur_list), ndim(mu.n_elem) {
    /* set the mean vector. */
    std::copy(mu.begin(), mu.end(), mu_vec);
  }

  void prep_sim
    (arma::mat const &sigma_permu, int const *indices, bool const is_permuted){
    // permute the mean vector, the covariance matrix, and the cur_list
    arma::mat tmp_mat(xtr_mem, ndim, ndim, false, true);
    if(!arma::chol(tmp_mat, sigma_permu))
      throw std::runtime_error("imputation::prep_sim: chol failed");
    copy_upper_tri(tmp_mat, sig_chol);

    if(!is_permuted)
      return;

    for(int i = 0; i < ndim; ++i)
      // we do not use this memory now as of this writing
      cdf_mem[i] = mu_vec[indices[i]];
    std::copy(cdf_mem, cdf_mem + ndim, mu_vec);
  }

  static inline int
  get_n_integrands(std::vector<type_base const *> const &ls) MDGC_NOEXCEPT {
    int out(1L);
    for(auto &x : ls)
      out += x->n_ele();
    return out;
  }

  inline int get_n_integrands() MDGC_NOEXCEPT {
    return n_integrands_val;
  }

  double * get_wk_mem(){
    return cdf_mem;
  }

  constexpr static bool needs_last_unif() {
    return true;
  }

  inline void operator()
  (double const * MDGC_RESTRICT draw, double * MDGC_RESTRICT out,
   int const *indices, bool const is_permutated, double const *w,
   int const n_draws){
    double * MDGC_RESTRICT scale_draw       = xtr_mem,
           * MDGC_RESTRICT scale_draw_permu = scale_draw + ndim * n_draws;
    std::fill(scale_draw, scale_draw + ndim * n_draws, 0.);

    {
      // re-scale
      double const * sig_chol_rp = sig_chol;
      double * MDGC_RESTRICT sc_samp = scale_draw;
      for(int c = 0; c < ndim; ++c, sc_samp += n_draws){
        double const * mult = draw;
        for(int r = 0; r <= c; ++r, ++sig_chol_rp, mult += n_draws)
          for(int k = 0; k < n_draws; ++k)
            sc_samp[k] += mult[k] * *sig_chol_rp;
      }
    }

    // re-locate
    for(int i = 0; i < ndim; ++i){
      double * MDGC_RESTRICT sc_samp = scale_draw + i * n_draws;
      for(int k = 0; k < n_draws; ++k)
        sc_samp[k] += mu_vec[i];
    }

    // permute and transpose
    for(int i = 0; i < ndim; ++i)
      for(int k = 0; k < n_draws; ++k)
        scale_draw_permu[k * ndim + indices[i]] =
          scale_draw[k + i * n_draws];

    std::fill(out, out + get_n_integrands(), 0);
    for(int k = 0; k < n_draws; ++k)
      *out += w[k];

    double const *scale_draw_i = scale_draw_permu;
    for(int k = 0; k < n_draws; ++k){
      double * o = out + 1;

      for(int i = 0; i < static_cast<int>(cur_list.size()); ++i){
        cur_list[i]->set_val(scale_draw_i, o, w[k]);

        // increment the pointers
        o            += cur_list[i]->n_ele();
        scale_draw_i += cur_list[i]->n_latent();
      }
    }
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
    // copy the result and return
    out_type out;
    out.minvls = minvls;
    out.inform = inform;
    out.abserr = abserr;
    out.likelihood = *res;

    arma::vec &imputations = out.imputations;
    imputations.set_size(std::max(0, get_n_integrands() - 1));
    std::copy(res + 1, res + 1 + imputations.n_elem, imputations.begin());

    return out;
  }
};
} // namespace restrictcdf

#endif
