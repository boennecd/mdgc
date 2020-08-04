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

/**
 * draws from a truncated normal distribution.
 *
 * @tparam i Code for finite upper bound. -1: unbounded, 0: -Inf lower
 * bound, 1: Inf upper bound, and 2: both are unbounded.
 * @param a The lower truncation point.
 * @param b The upper truncation point.
 * @param u Random uniform draw on (0, 1).
 * @param comp_quantile Boolean for whether to compute the CDF.
 */
template<int i>
std::array<double, 2> draw_trunc_mean
(double const a, double const b, double const u,
 bool const comp_quantile){
  throw std::runtime_error("draw_trunc_mean: not implemented");
  return { 0, 0 };
}

template<> inline std::array<double, 2> draw_trunc_mean<-1L>
(double const a, double const b, double const u,
 bool const comp_quantile){
  return { 1., u };
}

template<> inline std::array<double, 2> draw_trunc_mean<0L>
(double const a, double const b, double const u,
 bool const comp_quantile){
  double const qb = pnorm_std(b, 1L, 0L);
  if(comp_quantile)
    return { qb, qnorm_w(qb * u, 0, 1, 1L, 0L) };
  return { qb, std::numeric_limits<double>::quiet_NaN() };
}

template<> inline std::array<double, 2> draw_trunc_mean<1L>
(double const a, double const b, double const u,
 bool const comp_quantile){
  double const qa = pnorm_std(a, 1L, 0L);
  if(comp_quantile)
    return { 1 - qa, qnorm_w(qa + u * (1 - qa), 0, 1, 1L, 0L) };
  return { 1 - qa, std::numeric_limits<double>::quiet_NaN() };
}

template<> inline std::array<double, 2> draw_trunc_mean<2L>
(double const a, double const b, double const u,
 bool const comp_quantile){
  double const qa = pnorm_std(a, 1L, 0L),
               qb = pnorm_std(b, 1L, 0L);
#ifdef DO_CHECKS
  if(qb <= qa)
    throw std::runtime_error("draw_trunc_mean: qb <= qa");
#endif
  if(comp_quantile)
    return { qb - qa, qnorm_w(qa + u * (qb - qa), 0, 1, 1L, 0L) };
  return { qb - qa, std::numeric_limits<double>::quiet_NaN() };
}

/**
 * Holds output of the integral approximation.
 */
struct output {
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
  /// Maximum norm of estimated absolute accuracy of finest
  double abserr;
  /// Estimated NF-vector of values of the integrals.
  arma::vec finest;
};

typedef void (*mvkbrv_ptr)(int const*, double*, int const*, double*);

/**
 * Sets the function to be approximated.
 */
void set_mvkbrv_ptr(mvkbrv_ptr);

/**
 * copies the upper upper triangular matrix.
 *
 * @param X Matrix top copy.
 * @param x Pointer to copy to.
 */
inline void copy_upper_tri(arma::mat const &X, double *x){
  size_t const p = X.n_cols;
  for(unsigned c = 0; c < p; c++)
    for(unsigned r = 0; r <= c; r++, x++)
      *x = X.at(r, c);
}

/**
 * Approximates the function set by mvkbrv_ptr.
 *
 * @param ndim Dimension of the integral.
 * @param n_integrands Dimension of the integrand.
 * @param maxvls Maximum number of integrand evaluations.
 * @param abseps Absolute convergence threshold.
 * @param releps Relative convergence threshold.
 */
output approximate_integral(
    int const ndim, int const n_integrands, int const maxvls,
    double const abseps, double const releps);

/**
 * Approximates the integrals of the multivariate normal cdf over a
 * rectangle where the density is multiplied by some function in the
 * integrand.
 *
 * First, call this constructor. Then call the member functions to perform
 * the approximation.
 *
 * @tparam funcs Class with static member functions which determines the
 * type of integral which is approximated.
 *
 * It also needs the following static member functions.
 *
 * set_child_wk_mem should setup needed intermediaries.
 *
 * get_n_integrands should return the dimension of the integral.
 *
 * integrand should evaluate the integrand given a sample of truncated
 * normal variables.
 *
 * post_process should post process the result from the mvkbrv subroutine.
 * This may reduce the computation time.
 *
 * needs_last_unif returns a boolean for whether the last truncated normal
 * variable is needed.
 *
 * univariate handles the univariate case.
 */
template<class funcs>
class cdf {
  /** maps working memory to vectors */
  class ptr_to_dat {
    double * const wk_mem;
    int const ndim;

  public:
    /// objects for the parent clas
    double * const lower,
           * const upper,
           * const draw,
           * const tmp_vec,
           * const sigma_chol,
    /// objects for the child class
           * const child_mem;

    ptr_to_dat(double * const wk_mem, int const ndim):
      wk_mem(wk_mem), ndim(ndim),
      lower     (wk_mem            ),
      upper     (wk_mem + ndim     ),
      draw      (wk_mem + 2L * ndim),
      tmp_vec   (wk_mem + 3L * ndim),
      sigma_chol(wk_mem + 4L * ndim),
      child_mem (wk_mem + 4L * ndim + (ndim * (ndim + 1L)) / 2L)
      { }
  };

  static int ndim, n_integrands;
  static arma::ivec infin, idx;
  static double * wk_mem;
  static bool is_permutated;
  static constexpr bool const needs_last_unif =
    funcs::needs_last_unif();

#ifdef _OPENMP
#pragma omp threadprivate(ndim, n_integrands, infin, idx, wk_mem, is_permutated)
#endif

  static double * get_working_memory();

public:
  /**
   * Sets the working memory. Must be called prior to calling the other
   * member functions.
   *
   * @param max_dim Maximum number of variables to integrate out.
   * @param n_threads Maximum number of threads.
   */
  static void set_working_memory(size_t max_dim, size_t const n_threads);

  /**
   * Function to be called from mvkbrv.
   *
   * @param ndim_in Passed dimension of the integral.
   * @param unifs (0, 1) variables.
   * @param n_integrands_in Passed dimension of the integrand.
   * @param integrand_val Passed integrand approximation to be updated.
   */
  static void eval_integrand(
      int const *ndim_in, double *unifs, int const *n_integrands_in,
      double *integrand_val){
    size_t const udim = ndim;

#ifdef DO_CHECKS
    if(*ndim_in         != ndim)
      throw std::invalid_argument("cdf::eval_integrand: invalid 'ndim_in'");
    if(*n_integrands_in != n_integrands)
      throw std::invalid_argument("cdf::eval_integrand: invalid 'n_integrands_in'");
#endif

    ptr_to_dat map_obj(wk_mem, udim);
    arma::vec u(unifs        , udim        , false, false),
            out(integrand_val, n_integrands, false, false),
           draw(map_obj.draw , udim        , false, false);
    double const * const lower      = map_obj.lower,
                 * const upper      = map_obj.upper,
                 * const sigma_chol = map_obj.sigma_chol;

    double w(1.);
    double const * sc = sigma_chol,
                 * lw = lower,
                 * up = upper;
    /* loop over variables and transform them to truncated normal
     * variables */
    for(size_t j = 0; j < udim; ++j, ++sc, ++lw, ++up){
      auto const draw_n_p = ([&](){
        bool const needs_q =
          needs_last_unif or j + 1 < udim;
        double const * const draw_end = draw.begin() + j;

       if(infin[j] == 0L){
          double b(*up);
          for(double const *d = draw.begin(); d != draw_end; ){
            double const term = *sc++ * *d++;
            b -= term;
          }
          b /= *sc;

          return draw_trunc_mean<0L>(0, b, u[j], needs_q);

        } else if(infin[j] == 1L){
          double a(*lw);
          for(double const *d = draw.begin(); d != draw_end; ){
            double const term = *sc++ * *d++;
            a -= term;
          }
          a /= *sc;

          return draw_trunc_mean<1L>(a, 0, u[j], needs_q);

        } else if(infin[j] == 2L){
          double a(*lw),
                 b(*up);
          for(double const *d = draw.begin(); d != draw_end; ){
            double const term = *sc++ * *d++;
            a -= term;
            b -= term;
          }
          a /= *sc;
          b /= *sc;

          return draw_trunc_mean<2L>(a, b, u[j], needs_q);
        }
        else if(infin[j] == -1L){
          sc += j;
          return draw_trunc_mean<-1L>(0, 0, u[j], needs_q);

        }

        throw std::runtime_error("draw_trunc_mean: not implemented");
        sc += j;
        return draw_trunc_mean<-1L>(0, 0, u[j], needs_q);
      })();

      w       *= draw_n_p[0];
      draw[j]  = draw_n_p[1];
    }

    /* evaluate the integrand and weigth the result. */
    funcs::integrand(draw, ndim, out, map_obj.child_mem);
    out *= w;
  }

  /**
   * @param lower_in Lower bounds in the CDF.
   * @param upper_in Upper bounds in the CDF.
   * @param mu_in Mean vector.
   * @param sigma_in Covariance matrix.
   * @param do_reorder true if the order of integrations may be reordered.
   */
  cdf(arma::vec const &lower_in, arma::vec const &upper_in,
      arma::vec const &mu_in, arma::mat const &sigma_in,
      bool const do_reorder){
    ndim = mu_in.n_elem;
    n_integrands = funcs::get_n_integrands(mu_in, sigma_in);

    /* checks */
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

    size_t const udim = ndim;
    infin = pmvnorm::get_infin(lower_in, upper_in);
    wk_mem = get_working_memory();
    ptr_to_dat map_obj(wk_mem, udim);
    arma::vec lower     (map_obj.lower     , udim, false, false),
              upper     (map_obj.upper     , udim, false, false),
              sigma_chol(map_obj.sigma_chol, (udim * (udim + 1L)) / 2L, false, false);

    /* re-scale */
    arma::vec sds = arma::sqrt(arma::diagvec(sigma_in)),
               mu = mu_in / sds;

    lower  = lower_in;
    lower /= sds;
    lower -= mu;

    upper  = upper_in;
    upper /= sds;
    upper -= mu;

    is_permutated = false;
    if(do_reorder and udim > 1L){
      double * const y     = map_obj.draw,
             * const A     = map_obj.tmp_vec,
             * const B     = map_obj.child_mem,
             * const DL    = map_obj.child_mem + udim,
             * const delta = mu.begin();
      auto const correl = pmvnorm::get_cor_vec(sigma_in);
      int const pivot = 1L, doscale = 0L;
      int F_inform = 0L, nddim = udim;
      for(size_t i = 0; i < udim; ++i)
        *(delta + i) = 0.;

      idx.set_size(udim);
      for(size_t i = 0; i < udim; ++i)
        idx[i] = i;
      arma::ivec infi(udim);

      F77_CALL(mvsort)(
        &ndim, lower.memptr(), upper.memptr(), delta,
        correl.cor_vec.memptr(), infin.memptr(), y, &pivot, &nddim, A, B,
        DL, sigma_chol.memptr(), infi.memptr(), &F_inform, idx.memptr(),
        &doscale);

      if(F_inform != 0 or ndim != nddim)
        throw std::runtime_error("cdf::cdf: error in mvsort");

      int prev = idx[0];
      for(auto next = idx.begin() + 1L;
          next != idx.end() and !is_permutated; prev = *next, ++next){
        is_permutated = is_permutated or *next != prev + 1L;
      }

      if(is_permutated){
        for(size_t i = 0; i < udim; ++i){
          lower[i] = *(A + i);
          upper[i] = *(B + i);
          infin[i] = infi[i];
        }

        arma::uvec uidx(udim);
        arma::vec mu_permu(map_obj.tmp_vec, udim, false, false);
        for(size_t i = 0; i < udim; ++i){
          uidx[i]     = idx[i];
          mu_permu[i] = mu[uidx[i]];
        }
        arma::mat sigma_permu = sigma_in.submat(uidx, uidx);
        funcs::set_child_wk_mem(mu_permu, sigma_permu, map_obj.child_mem);
        return;
      }

    } else if(!do_reorder and udim > 1L) {
      arma::mat tmp = sigma_in;
      tmp.each_row() /= sds.t();
      tmp.each_col() /= sds;
      if(!arma::chol(tmp, tmp))
        sigma_chol.fill(std::numeric_limits<double>::infinity());
      else
        copy_upper_tri(tmp, sigma_chol.memptr());

    }

    funcs::set_child_wk_mem(mu_in, sigma_in, map_obj.child_mem);
  }

  /**
   * Integrates over (-Inf, 0)^[# of random effects]
   *
   * @param mu_in Mean vector.
   * @param sigma_in Covariance matrix.
   * @param do_reorder true if the order of integrations may be reordered.
   */
  cdf(arma::vec const &mu_in, arma::mat const &sigma_in,
      bool const do_reorder):
  cdf(
    ([&](){
      arma::vec out(mu_in.n_elem);
      out.fill(-std::numeric_limits<double>::infinity());
      return out;
    })(), arma::vec(mu_in.n_elem, arma::fill::zeros), mu_in, sigma_in,
    do_reorder) { }

  /**
   * Performs the approximation.
   *
   * @param maxvls Maximum number of function evaluations allowed.
   * @param abseps Required absolute accuracy.
   * @param releps Required relative accuracy.
   *
   * @return Either a one-dimension or a (1 + p + p(p + 1)/2)-dimensional
   * vector. The latter contains the likelihood approximation, derivative
   * with respect to the mean, and a (p(p + 1) /2)-dimensional vector
   * containing an upper diagonal matrix with derivatives with respect to
   * the covariance matrix. Notice that the latter are not scaled by 2.
   */
  static output approximate
  (int const maxvls, double const abseps, double const releps){
#ifdef DO_CHECKS
    if(abseps <= 0 and releps <= 0)
      throw std::invalid_argument("cdf::approximate: no valid convergence threshold");
    if(maxvls <= 0L)
      throw std::invalid_argument("cdf::approximate: invalid 'maxvls'");
#endif

    ptr_to_dat map_obj(wk_mem, ndim);
    if(ndim == 1L){
      /* handle the one-dimensional case as a special case */
      double const lw = *map_obj.lower,
                   up = *map_obj.upper;
      output out;
      out.finest = funcs::univariate(lw, up, map_obj.child_mem);
      out.inform = 0L;
      out.abserr = 0;
      return out;

    } else if(std::isinf(*map_obj.sigma_chol)){
      output out;
      out.finest.resize(n_integrands);
      out.finest.fill(std::numeric_limits<double>::quiet_NaN());
      out.inform = -1L;
      return out;

    }

    /* set pointer to this class' member function */
    wk_mem = get_working_memory();
    set_mvkbrv_ptr(&cdf<funcs>::eval_integrand);
    output out =
      approximate_integral(ndim, n_integrands, maxvls, abseps, releps);

    funcs::post_process(out.finest, ndim, map_obj.child_mem);
    if(is_permutated)
      funcs::permutate(out.finest, ndim, idx, wk_mem);

    return out;
  }
};

/**
 * func classes used as template argument for cdf used to approximate the
 * likelihood. */
class likelihood {
public:
  static void set_child_wk_mem
  (arma::vec const&, arma::mat const&, double const * const) { }

  static int constexpr get_n_integrands(arma::vec const&, arma::mat const&){
    return 1L;
  }
  static void integrand(arma::vec const&, int const, arma::vec&,
                        double const * const);
  static void post_process(arma::vec&, int const, double const * const) { }
  constexpr static bool needs_last_unif() {
    return false;
  }

  static arma::vec univariate(double const lw, double const ub,
                              double const * const wk_mem){
    arma::vec out(1L);
    double const p_ub = std::isinf(ub) ? 1 : pnorm_std(ub, 1L, 0L),
                 p_lb = std::isinf(lw) ? 0 : pnorm_std(lw, 1L, 0L);
    out[0] = p_ub - p_lb;
    return out;
  }

  static void permutate(arma::vec&, int const, arma::ivec const&,
                        double * const) { }
};

/**
 * func classes used as template argument for cdf used to approximates both
 * the probability and the derivatives w.r.t. the mean and the
 * covariance matrix. */
class deriv {
public:
  static void set_child_wk_mem
  (arma::vec const &mu, arma::mat const &sigma,
   double * const wk_mem){
    size_t const p = mu.n_elem;
    arma::mat signa_inv(arma::inv(sigma)),
    sigma_chol_inv(arma::inv(arma::trimatu(arma::chol(sigma))));

    copy_upper_tri(sigma_chol_inv, wk_mem);
    copy_upper_tri(signa_inv     , wk_mem + (p * (p + 1L)) / 2L);
  }

  static int get_n_integrands(arma::vec const&, arma::mat const&);
  static void integrand(arma::vec const&, int const, arma::vec&,
                        double const * const);
  static void post_process(arma::vec&, int const, double const * const);
  constexpr static bool needs_last_unif() {
    return true;
  }

  static arma::vec univariate(double const lw, double const ub,
                              double const * const wk_mem){
    arma::vec out(3L);
    static double const sqrt_2_pi = std::sqrt(2 * M_PI);
    auto dnrm = [&](double const x){
      return std::exp(-x * x / 2.) / sqrt_2_pi;
    };

    bool const f_ub = std::isinf(ub),
               f_lb = std::isinf(lw);

    double const p_ub = f_ub ? 1 : pnorm_std(ub, 1L, 0L),
                 p_lb = f_lb ? 0 : pnorm_std(lw, 1L, 0L),
               sig_inv = *wk_mem,
                  d_ub = f_ub ? 0 : dnrm(ub),
                  d_lb = f_lb ? 0 : dnrm(lw),
               d_ub_ub = f_ub ? 0 : ub * d_ub,
               d_lb_lb = f_lb ? 0 : lw * d_lb;
    out[0L] = p_ub - p_lb;
    out[1L] = -(d_ub - d_lb) * sig_inv;
    out[2L] = -(d_ub_ub - d_lb_lb) / 2 * sig_inv * sig_inv;
    return out;
  }

  static void permutate(arma::vec &finest, int const ndim,
                        arma::ivec const &idx, double * const work_mem) {
    size_t const dim_cov = (ndim * (ndim + 1L)) / 2L;

    arma::vec dmu(finest.memptr() + 1L       , ndim   , false, true),
             dcov(finest.memptr() + 1L + ndim, dim_cov, false, true),
              cmu(work_mem                   , ndim   , false, true),
             ccov(work_mem + ndim            , dim_cov, false, true);

    /* maps to the indices of the upper triangle matrix */
    auto tri_map = [&](int const r, int const c){
      return r + (c * (c + 1L)) / 2L;
    };

    /* copy to temporaries */
    for(int c = 0; c < ndim; ++c){
      int const org_c = idx[c];
      cmu[org_c] = dmu[c];

      for(int r = 0; r <= c; ++r){
        double const val = dcov[tri_map(r, c)];
        int const org_r = idx[r];
        if(org_c >= org_r)
          ccov[tri_map(org_r, org_c)] = val;
        else
          ccov[tri_map(org_c, org_r)] = val;
      }
    }

    /* copy back */
    for(int c = 0; c < ndim; ++c)
      dmu[c] = cmu[c];
    for(size_t i = 0; i < dim_cov; ++i)
      dcov[i] = ccov[i];
  }
};

/* initialize static members */
template<class funcs>
int cdf<funcs>::ndim = 0L;
template<class funcs>
int cdf<funcs>::n_integrands = 0L;
template<class funcs>
arma::ivec cdf<funcs>::infin = arma::ivec();
template<class funcs>
arma::ivec cdf<funcs>::idx = arma::ivec();
template<class funcs>
double * cdf<funcs>::wk_mem = nullptr;
template<class funcs>
bool cdf<funcs>::is_permutated = false;
}

#endif
