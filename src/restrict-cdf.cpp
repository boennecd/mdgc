#include "restrict-cdf.h"
#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif

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
  // assume the cacheline is 128 bytes. Then make sure we avoid false
  // sharing by having 2 x 128 bytes per thread.
  constexpr size_t const mult = 128L / 8L,
                     min_size = 2L * mult;
  size_t const upper_tri_size = (max_dim * (max_dim + 1L)) / 2L;
  max_dim = 4L * max_dim + 3L * upper_tri_size;
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
(arma::vec const &draw, int const ndim, arma::vec &out,
 double const * const){
#ifdef DO_CHECKS
  if(out.n_elem != 1L)
    throw invalid_argument("likelihood::integrand: invalid out");
#endif
  out[0] = 1;
}

int deriv::get_n_integrands
(arma::vec const &mu, arma::mat const &sigma) {
  arma::uword const p = mu.n_elem;
  return 1 + p + (p * (p + 1)) / 2L;
}

void deriv::integrand
(arma::vec const &draw, int const ndim, arma::vec &out,
 double const * const wk_mem){
  arma::uword const p = ndim;

#ifdef DO_CHECKS
  size_t const n_elem = 1L + p + (p * (p + 1L)) / 2L;
  if(out.n_elem != n_elem)
    throw invalid_argument("deriv::integrand: invalid out");
#endif
  out.zeros();

  out[0L] = 1.;
  double * const mean_part_begin = out.memptr() + 1L;
  /* Multiplying by the inverse matrix is fast but not smart numerically.
   * TODO: much of this computation can be done later */
  double const * sigma_chol_inv = wk_mem;
  for(unsigned c = 0; c < p; ++c){
    double const mult = draw[c],
          * const end = mean_part_begin + c + 1L;
    for(double *rhs = mean_part_begin; rhs != end; ++rhs, ++sigma_chol_inv)
      *rhs += mult * *sigma_chol_inv;
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

void deriv::post_process(arma::vec &finest, int const ndim,
                         double const * const wk_mem) {
  arma::uword const p = ndim;

  double phat = finest[0L];
  double *o = finest.memptr() + 1L + p;
  double const * sig_inv = wk_mem + (p * (p + 1L)) / 2L;
  for(unsigned c = 0; c < p; c++){
    double const * const end = sig_inv + c + 1L;
    for(; sig_inv != end; ++sig_inv, ++o){
      *o -= phat * *sig_inv;
      *o /= 2.;
    }
  }
}

template class cdf<likelihood>;
template class cdf<deriv>;
}
