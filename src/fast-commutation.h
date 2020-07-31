#ifndef FAST_COM_H
#define FAST_COM_H
#include "arma-wrap.h"

inline arma::mat get_commutation_unequal
  (unsigned const n, unsigned const m){
  unsigned const nm = n * m,
             nnm_p1 = n * nm + 1L,
              nm_pm = nm + m;
  arma::mat out(nm, nm, arma::fill::zeros);
  double * o = out.begin();
  for(unsigned i = 0; i < n; ++i, o += nm_pm){
    double *o1 = o;
    for(unsigned j = 0; j < m; ++j, o1 += nnm_p1)
      *o1 = 1.;
  }

  return out;
}

inline arma::mat get_commutation_equal(unsigned const m){
  unsigned const mm = m * m,
                mmm = mm * m,
             mmm_p1 = mmm + 1L,
              mm_pm = mm + m;
  arma::mat out(mm, mm, arma::fill::zeros);
  double * const o = out.begin();
  unsigned inc_i(0L);
  for(unsigned i = 0; i < m; ++i, inc_i += m){
    double *o1 = o + inc_i + i * mm,
           *o2 = o + i     + inc_i * mm;
    for(unsigned j = 0; j < i; ++j, o1 += mmm_p1, o2 += mm_pm){
      *o1 = 1.;
      *o2 = 1.;
    }
    *o1 += 1.;
  }
  return out;
}

inline arma::mat get_commutation(unsigned const n, unsigned const m){
  if(n == m)
    return get_commutation_equal(n);

  return get_commutation_unequal(n, m);
}

#endif
