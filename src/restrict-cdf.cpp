#include "restrict-cdf.h"
#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif

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
   * @param inform non-zero if something went wrong.
   * @param idx N-dimensional vector with indices of applied permutation.
   */
  void F77_NAME(mvsqrt)(
      int const* /* N */, double const* /* lower */,
      double const* /* upper */, double const* /* delta */,
      double const* /* correl */, int const* /* infin */,
      double const* /* y */, int const* /* pivot */,
      int* /* ND */, double* /* A */, double* /* B */, double* /* DL */,
      double* /* cov */, int* /* infi */, int* /* inform */,
      int* /* idx */);
}

using std::invalid_argument;

static restrictcdf::mvkbrv_ptr current_mvkbrv_ptr = nullptr;
#ifdef _OPENMP
#pragma omp threadprivate(current_mvkbrv_ptr)
#endif

void restrictcdf::set_mvkbrv_ptr(mvkbrv_ptr new_ptr){
  current_mvkbrv_ptr = new_ptr;
}

extern "C"
{
  void F77_NAME(mvkbrveval)(
      int const* /* NDIM */, int const* /* MAXVLS */, int const* /* NF */,
      double const* /* ABSEPS */, double const* /* RELEPS */,
      double* /* ABSERR */, double* /* FINEST */, int* /* INFORM */);

  void F77_SUB(mvkbrvintegrand)
    (int const *m, double *unifs, int const *mf, double *out){
#ifdef DO_CHECKS
    if(!current_mvkbrv_ptr)
      throw invalid_argument("mvkbrvintegrand: 'current_mvkbrv_ptr' not set");
    if(!out)
      throw invalid_argument("mvkbrvintegrand: 'out' not set");
#endif
    (*current_mvkbrv_ptr)(m, unifs, mf, out);
  }
}

namespace restrictcdf {
static size_t wk_mem_per_thread = 0L,
              current_wk_size   = 0L;
static std::unique_ptr<double[]> current_wk_mem =
  std::unique_ptr<double[]>();

template<class funcs>
double * cdf<funcs>::get_working_memory(){
#ifdef _OPENMP
  size_t const my_num = omp_get_thread_num();
#else
  size_t const my_num(0L);
#endif

  return current_wk_mem.get() + my_num * wk_mem_per_thread;
}

template<class funcs>
void cdf<funcs>::set_working_memory
  (size_t max_dim, size_t const n_threads){
  // assume thte cacheline is 128 bytes. Then make sure we avoid false
  // sharing by having 2 x 128 bytes per thread.
  constexpr size_t const mult = 128L / 8L,
                     min_size = 2L * mult;
  max_dim = std::max(max_dim, min_size);
  max_dim = (max_dim + mult - 1L) / mult;
  max_dim *= mult;
  wk_mem_per_thread = max_dim;

  size_t const new_size =
    std::max(n_threads, static_cast<size_t>(1L)) * max_dim;
  if(new_size > current_wk_size){
    current_wk_mem.reset(new double[new_size]);
    current_wk_size = new_size;
    return;
  }
}

output approximate_integral(
    int const ndim, int const n_integrands, int const maxvls,
    double const abseps, double const releps){
  output out;
  out.finest.resize(n_integrands);

  F77_CALL(mvkbrveval)(
      &ndim, &maxvls, &n_integrands, &abseps, &releps,
      &out.abserr, out.finest.memptr(), &out.inform);

  return out;
}

void likelihood::integrand
(arma::vec const &draw, likelihood::comp_dat const& dat, arma::vec &out){
#ifdef DO_CHECKS
  if(out.n_elem != 1L)
    throw invalid_argument("likelihood::integrand: invalid out");
#endif
  out[0] = 1;
}

void likelihood::post_process(arma::vec &finest, comp_dat const &dat) { }

int deriv::get_n_integrands
(arma::vec const &mu, arma::mat const &sigma) {
  arma::uword const p = mu.n_elem;
  return 1 + p + (p * (p + 1)) / 2L;
}

void deriv::integrand
(arma::vec const &draw, deriv::comp_dat const& dat, arma::vec &out){
  arma::uword const p = dat.mu->n_elem;

#ifdef DO_CHECKS
  size_t const n_elem = 1L + p + (p * (p + 1L)) / 2L;
  if(!dat.mu)
    throw std::runtime_error("deriv::integrand: mu is not set");
  if(out.n_elem != n_elem)
    throw invalid_argument("deriv::integrand: invalid out");
#endif
  out.zeros();

  out[0L] = 1.;
  double * const mean_part_begin = out.memptr() + 1L;
  /* Multiplying by the inverse matrix is fast but not smart numerically.
   * TODO: much of this computation can be done later */
  for(unsigned c = 0; c < p; ++c){
    double const mult = draw[c],
                  *r1 = dat.sigma_chol_inv.colptr(c),
          * const end = mean_part_begin + c + 1L;
    for(double *rhs = mean_part_begin; rhs != end; ++r1, ++rhs)
      *rhs += mult * *r1;
  }

  {
    double *o = out.memptr() + 1L + p;
    for(unsigned c = 0; c < p; c++){
      double const mult = *(mean_part_begin + c),
            * const end = mean_part_begin + c + 1;
      for(double *lhs = mean_part_begin; lhs != end; ++o, ++lhs)
        *o = mult * *lhs;
    }
  }
}

void deriv::post_process(arma::vec &finest, comp_dat const &dat) {
  arma::uword const p = dat.mu->n_elem;
  auto const &sig_inv = dat.signa_inv;

  double phat = finest[0L];
  double *o = finest.memptr() + 1L + p;
  for(unsigned c = 0; c < p; c++){
    double const * const end = sig_inv.colptr(c) + c + 1L;
    for(auto rhs = sig_inv.colptr(c); rhs != end; ++rhs, ++o){
      *o -= phat * *rhs;
      *o /= 2.;
    }
  }
}

template class cdf<likelihood>;
template class cdf<deriv>;
}
