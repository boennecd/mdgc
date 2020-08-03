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
           * const sigma_chol,
    /// objects for the child class
           * const child_mem;


    ptr_to_dat(double * const wk_mem, int const ndim):
      wk_mem(wk_mem), ndim(ndim),
      lower     (wk_mem            ),
      upper     (wk_mem + ndim     ),
      draw      (wk_mem + 2L * ndim),
      sigma_chol(wk_mem + 3L * ndim),
      child_mem (wk_mem + 3L * ndim + (ndim * (ndim + 1L)) / 2L)
      { }
  };

  static int ndim, n_integrands;
  static arma::ivec infin;
  static double * wk_mem;
  static constexpr bool const needs_last_unif =
    funcs::needs_last_unif();

#ifdef _OPENMP
#pragma omp threadprivate(ndim, n_integrands, infin, wk_mem)
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
#ifdef DO_CHECKS
    if(*ndim_in         != ndim)
      throw std::invalid_argument("cdf::eval_integrand: invalid 'ndim_in'");
    if(*n_integrands_in != n_integrands)
      throw std::invalid_argument("cdf::eval_integrand: invalid 'n_integrands_in'");
#endif

    ptr_to_dat map_obj(wk_mem, ndim);
    arma::vec u(unifs        , ndim        , false, false),
            out(integrand_val, n_integrands, false, false),
           draw(map_obj.draw , ndim        , false, false);
    double const * const lower      = map_obj.lower,
                 * const upper      = map_obj.upper,
                 * const sigma_chol = map_obj.sigma_chol;

    double w(1.);
    double const * sc = sigma_chol,
                 * lw = lower,
                 * up = upper;
    /* loop over variables and transform them to truncated normal
     * variables */
    for(size_t j = 0; j < static_cast<size_t>(ndim); ++j, ++sc, ++lw, ++up){
      auto const draw_n_p = ([&](){
        bool const needs_q =
          needs_last_unif or j + 1 < static_cast<size_t>(ndim);
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

    infin = pmvnorm::get_infin(lower_in, upper_in);
    wk_mem = get_working_memory();
    ptr_to_dat map_obj(wk_mem, ndim);
    arma::vec lower     (map_obj.lower     , ndim, false, false),
              upper     (map_obj.upper     , ndim, false, false),
              sigma_chol(map_obj.sigma_chol, (ndim * (ndim + 1L)) / 2L, false, false);

    /* re-scale */
    arma::vec const sds = arma::sqrt(arma::diagvec(sigma_in)),
                    mu = mu_in / sds;

    lower  = lower_in;
    lower /= sds;
    lower -= mu;

    upper  = upper_in;
    upper /= sds;
    upper -= mu;

    {
      arma::mat tmp = sigma_in;
      tmp.each_row() /= sds.t();
      tmp.each_col() /= sds;
      if(!arma::chol(tmp, tmp))
        sigma_chol += std::numeric_limits<double>::infinity();
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
   */
  cdf(arma::vec const &mu_in, arma::mat const &sigma_in):
  cdf(
    ([&](){
      arma::vec out(mu_in.n_elem);
      out.fill(-std::numeric_limits<double>::infinity());
      return out;
    })(), arma::vec(mu_in.n_elem, arma::fill::zeros), mu_in, sigma_in,
    false) { }

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
    if(std::isinf(map_obj.sigma_chol[0])){
      output out;
      out.finest.resize(n_integrands);
      out.finest.fill(std::numeric_limits<double>::quiet_NaN());
      out.inform = -1L;
      return out;
    }

    // TODO: there is no variable reordering at the moment. This may reduce
    // the variance of the estimators.

    /* set pointer to this class' member function */
    wk_mem = get_working_memory();
    set_mvkbrv_ptr(&cdf<funcs>::eval_integrand);
    output out =
      approximate_integral(ndim, n_integrands, maxvls, abseps, releps);

    funcs::post_process(out.finest, ndim, map_obj.child_mem);
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
};

/* initialize static members */
template<class funcs>
int cdf<funcs>::ndim = 0L;
template<class funcs>
int cdf<funcs>::n_integrands = 0L;
template<class funcs>
arma::ivec cdf<funcs>::infin = arma::ivec();
template<class funcs>
double * cdf<funcs>::wk_mem = nullptr;
}

#endif
