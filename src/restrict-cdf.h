#ifndef RESTRICT_CDF_H
#define RESTRICT_CDF_H

#include "arma-wrap.h"
#include <array>
#include <limits>
#include <memory>
#include <cmath>
#include "pnorm.h"
#include "qnorm.h"

namespace restrictcdf {

/**
 * draws from a truncated normal distribution.
 *
 * @param b The upper truncation point.
 * @param u Random uniform draw on (0, 1).
 * @param comp_quantile Boolean for whether to compute the CDF.
 */
inline std::array<double, 2> draw_trunc_mean
  (double const b, const double u, bool comp_quantile){
  double const qb = pnorm_std(b, 1L, 0L);
  if(comp_quantile)
    return { qb, qnorm_w(qb * u, 0, 1, 1L, 0L) };
  return { qb, std::numeric_limits<double>::quiet_NaN() };
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
 * sets the function to be approximated.
 */
void set_mvkbrv_ptr(mvkbrv_ptr);

/**
 * approximates the function set by mvkbrv_ptr.
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
 * rectangle from (-Inf, 0)^[# of random effects] where the PDF is
 * multiplied by some function in the integrand.
 *
 * First, call this constructor. Then call the member functions to perform
 * the approximation.
 *
 * @tparam funcs Class with static member functions which determines the
 * type of integral which is approximated.
 *
 * The funcs class needs a nested class comp_dat to store the data needed
 * to perform the approximation. It also needs the following static member
 * functions:
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
  using comp_dat = typename funcs::comp_dat;

  static int ndim, n_integrands;
  static arma::vec mu;
  static arma::vec sigma_chol;
  /* TODO: use omp threadprivate... */
  thread_local static std::unique_ptr<comp_dat> dat;
  static constexpr bool const needs_last_unif =
    funcs::needs_last_unif();

#ifdef _OPENMP
#pragma omp threadprivate(ndim, n_integrands, mu, sigma_chol)
#endif

public:
  /**
   * function to be called from mvkbrv.
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

    arma::vec u(unifs        , ndim        , false),
            out(integrand_val, n_integrands, false),
           draw(ndim);

    double w(1.), *sc = sigma_chol.memptr();
    /* loop over variables and transform them to truncated normal
     * variables */
    for(unsigned j = 0; j < (unsigned)ndim; ++j){
      double b(-mu[j]);
      for(unsigned k = 0; k < j; ++k)
        b -= *sc++ * draw[k];
      b /= *sc++;

      auto const draw_n_p = draw_trunc_mean
        (b, u[j], needs_last_unif or j + 1 < (unsigned)ndim);
      w       *= draw_n_p[0];
      draw[j]  = draw_n_p[1];
    }

    /* evaluate the integrand and weigth the result. */
    funcs::integrand(draw, *dat, out);
    out *= w;
  }

  /**
   * @param mu_in Mean vector.
   * @param sigma_in Covariance matrix.
   */
  cdf(arma::vec const &mu_in, arma::mat const &sigma_in){
    ndim = mu_in.n_elem;
    n_integrands = funcs::get_n_integrands(mu_in, sigma_in);

    /* checks */
#ifdef DO_CHECKS
    if(sigma_in.n_cols != static_cast<size_t>(ndim) or
         sigma_in.n_rows != static_cast<size_t>(ndim))
      throw std::invalid_argument("cdf::cdf: invalid 'sigma_in'");
    if(n_integrands <= 0L)
      throw std::invalid_argument("cdf::cdf: invalid 'n_integrands'");
#endif

    /* re-scale */
    arma::vec const sds = arma::sqrt(arma::diagvec(sigma_in));
    mu = mu_in / sds;

    sigma_chol = ([&]{
      arma::uword const p = sigma_in.n_cols;
      arma::vec out((p * (p + 1L)) / 2L);


      arma::mat tmp = sigma_in;
      tmp.each_row() /= sds.t();
      tmp.each_col() /= sds;
      if(!arma::chol(tmp, tmp)){
        out += std::numeric_limits<double>::infinity();
        return out;
      }

      double *o = out.memptr();
      for(unsigned c = 0; c < p; c++)
        for(unsigned r = 0; r <= c; r++)
          *o++ = tmp.at(r, c);

      return out;
    })();

    dat.reset(new comp_dat(mu_in, sigma_in, sigma_chol));
  }

  /**
   * Performs the approximation.
   *
   * @param maxvls Maximum number of function evaluations allowed.
   * @param abseps Required absolute accuracy.
   * @param releps Required relative accuracy.
   */
  static output approximate
  (int const maxvls, double const abseps, double const releps){
#ifdef DO_CHECKS
    if(abseps <= 0 and releps <= 0)
      throw std::invalid_argument("cdf::approximate: no valid convergence threshold");
    if(maxvls <= 0L)
      throw std::invalid_argument("cdf::approximate: invalid 'maxvls'");
#endif

    if(std::isinf(sigma_chol[0])){
      output out;
      out.finest.resize(n_integrands);
      out.finest.fill(std::numeric_limits<double>::quiet_NaN());
      out.inform = -1L;
      return out;
    }

    /* set pointer to this class' member function */
    set_mvkbrv_ptr(&cdf<funcs>::eval_integrand);
    output out =
      approximate_integral(ndim, n_integrands, maxvls, abseps, releps);

    funcs::post_process(out.finest, *dat);
    return out;
  }
};

/**
 * func classes used as template argument for cdf used to approximate the
 * likelihood. */
class likelihood {
public:
  class comp_dat {
  public:
    comp_dat(arma::vec const&, arma::mat const&, arma::vec const&) { }
  };

  static int get_n_integrands(arma::vec const&, arma::mat const&);
  static void integrand(arma::vec const&, comp_dat const&,arma::vec&);
  static void post_process(arma::vec&, comp_dat const&);
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
  /**
   * stores data need to perform the approximation.
   */
  class comp_dat {
  public:
    arma::vec const *mu;
    arma::mat const *sigma,
                     signa_inv,
                     sigma_chol_inv;

    comp_dat(arma::vec const &mu_in, arma::mat const &sigma_in,
             arma::vec const &sigma_chol):
      mu(&mu_in), sigma(&sigma_in), signa_inv(arma::inv(sigma_in)),
      sigma_chol_inv(arma::inv(arma::trimatu(
        arma::chol(sigma_in)))) { }
  };

  static int get_n_integrands(arma::vec const&, arma::mat const&);
  static void integrand(arma::vec const&, comp_dat const&, arma::vec&);
  static void post_process(arma::vec&, comp_dat const&);
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
arma::vec cdf<funcs>::mu = arma::vec();
template<class funcs>
arma::vec cdf<funcs>::sigma_chol = arma::vec();
template<class funcs>
thread_local std::unique_ptr<typename cdf<funcs>::comp_dat >
  cdf<funcs>::dat = std::unique_ptr<cdf<funcs>::comp_dat>();
}

#endif
