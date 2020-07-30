#include "restrict-cdf.h"

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
  {
    for(unsigned c = 0; c < p; ++c)
      for(unsigned r = 0; r <= c; ++r)
        *(mean_part_begin + r) += dat.sigma_chol_inv.at(r, c) * draw[c];
  }

  {
    double *o = out.memptr() + 1L + p;
    for(unsigned c = 0; c < p; c++){
      double const mean_part_c = *(mean_part_begin + c);
      for(unsigned r = 0; r <= c; r++)
        *o++ = mean_part_c * *(mean_part_begin + r);
    }
  }
}

void deriv::post_process(arma::vec &finest, comp_dat const &dat) {
  arma::uword const p = dat.mu->n_elem;

  double phat = finest[0L];
  double *o = finest.memptr() + 1L + p;
  for(unsigned c = 0; c < p; c++)
    for(unsigned r = 0; r <= c; r++){
      *o   -= phat * dat.signa_inv(r, c);
      *o++ /= 2.;
    }
}

template class cdf<likelihood>;
template class cdf<deriv>;
}
