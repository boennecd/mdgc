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
      double* /* ABSERR */, double* /* FINEST */, int* /* INFORM */,
      int* /* MINVLS */ );

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
static size_t wk_mem_per_thread  = 0L,
              current_wk_size    = 0L,
              iwk_mem_per_thread = 0L,
              current_iwk_size   = 0L;
static std::unique_ptr<double[]> current_wk_mem =
  std::unique_ptr<double[]>();
static std::unique_ptr<int[]> current_iwk_mem =
  std::unique_ptr<int[]>();

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
int * cdf<funcs>::get_iworking_memory(){
#ifdef _OPENMP
  size_t const my_num = omp_get_thread_num();
#else
  size_t const my_num(0L);
#endif

  return current_iwk_mem.get() + my_num * iwk_mem_per_thread;
}

template<class funcs>
void cdf<funcs>::set_working_memory
  (size_t const max_dim, size_t const n_threads){
  constexpr size_t const cachline_size = 128L;

  {
    // makes sure we avoid false sharing.
    constexpr size_t const mult = cachline_size / sizeof(double),
                       min_size = 2L * mult;
    size_t const upper_tri_size = (max_dim * (max_dim + 1L)) / 2L;
    size_t m_dim =
      4L * max_dim + 2L * max_dim * max_dim + 3L * upper_tri_size;
    m_dim = std::max(m_dim, min_size);
    m_dim = (m_dim + mult - 1L) / mult;
    m_dim *= mult;
    wk_mem_per_thread = m_dim;

    size_t const new_size =
      std::max(n_threads, static_cast<size_t>(1L)) * m_dim;
    if(new_size > current_wk_size){
      current_wk_mem.reset(new double[new_size]);
      current_wk_size = new_size;

    }
  }

  {
    // makes sure we avoid false sharing.
    constexpr size_t const mult = cachline_size / sizeof(int),
                       min_size = 2L * mult;
    size_t m_dim = 2L * max_dim;
    m_dim = std::max(m_dim, min_size);
    m_dim = (m_dim + mult - 1L) / mult;
    m_dim *= mult;
    iwk_mem_per_thread = m_dim;

    size_t const new_size =
      std::max(n_threads, static_cast<size_t>(1L)) * m_dim;
    if(new_size > current_iwk_size){
      current_iwk_mem.reset(new int[new_size]);
      current_iwk_size = new_size;

    }
  }
}

output approximate_integral(
    int const ndim, int const n_integrands, int const maxvls,
    double const abseps, double const releps){
  output out;
  out.finest.resize(n_integrands);
  out.minvls = 0L;

  F77_CALL(mvkbrveval)(
      &ndim, &maxvls, &n_integrands, &abseps, &releps,
      &out.abserr, out.finest.memptr(), &out.inform,
      &out.minvls);

  return out;
}

void likelihood::integrand
(double const * const, int const ndim, double * const out,
 double const * const){
#ifdef DO_CHECKS
  if(!out)
    throw invalid_argument("likelihood::integrand: invalid out");
#endif
  *out = 1;
}

int deriv::get_n_integrands
(arma::vec const &mu, arma::mat const &sigma) {
  arma::uword const p = mu.n_elem;
  return 1 + p + (p * (p + 1)) / 2L;
}

void deriv::integrand
(double const * const draw, int const ndim, double * const out,
 double const * const wk_mem){
  arma::uword const p = ndim;

  size_t const n_elem = 1L + p + (p * (p + 1L)) / 2L;
  double * const out_end = out + n_elem;

#ifdef DO_CHECKS
  for(double * o = out + 1L; o != out_end; ++o)
    if(!o)
      throw invalid_argument("likelihood::integrand: invalid out");
#endif

  for(double * o = out + 1L; o != out_end; ++o)
    *o = 0.;

  *out = 1.;
  double * const mean_part_begin = out + 1L;
  /* Multiplying by the inverse matrix is fast but not smart numerically.
   * TODO: much of this computation can be done later */
  double const * sigma_chol_inv = wk_mem;
  for(unsigned c = 0; c < p; ++c){
    double const mult = *(draw + c),
          * const end = mean_part_begin + c + 1L;
    for(double *rhs = mean_part_begin; rhs != end; ++rhs, ++sigma_chol_inv)
      *rhs += mult * *sigma_chol_inv;
  }

  {
    double *o = out + 1L + p;
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
