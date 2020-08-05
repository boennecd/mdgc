#ifndef LP_UTILS_H
#define LP_UTILS_H
#include "arma-wrap.h"

/**
 * computes (X (x) X) x where x is p^2 vector and X is a k x p matrix and
 * store it in the passed pointer which can hold k^2 elements.
 *
 * @param wk working memory with up k x p elements.
 */
inline void X_kron_X_dot_x
(arma::mat const &X, arma::vec const &x, double * const out,
 double * const wk){
  size_t const k = X.n_rows,
               p = X.n_cols,
              kk = k * k;
#ifdef DDO_CHECKS
  if(x.n_elem != p * p)
    throw std::invalid_argument("X_kron_X_dot_x: invalid out");
#endif
  for(size_t i = 0; i < kk; ++i)
    *(out + i) = 0;

  double * const w_end = wk + k * p;
  for(size_t c = 0; c < p; ++c){
    size_t const cp = c * p;

    {
      double const * Xp = X.memptr(),
                   * const xp_end = x.memptr() + cp + p;
      double * w = wk;
      for(double const * xp = x.memptr() + cp; xp != xp_end; ++xp){
        double const mult = *xp,
              * const w_end = w + k;
        for(; w != w_end; ++w, ++Xp)
          *w = mult * *Xp;
      }
    }

    for(size_t r = 0; r < k; ++r){
      size_t const kr = k * r;
      double const mult = X.at(r, c);
      double * const o_end = out + k + kr;
      for(auto w = wk; w != w_end; )
        for(double * o = out + kr; o != o_end; ++o, ++w)
          *o += mult * *w;
    }
  }
}

/**
 * computes (x (x) X) y where x is l vector, X is a k x p matrix, 
 * y is a p vector, and stores it in the passed pointer which can 
 * hold l x k elements.
 * 
 * @param wk working memory with up k elements.
 */
inline void x_kron_X_dot_y
  (arma::vec const &x, arma::mat const &X, arma::vec const &y, 
   double * const out, double * const wk){
  size_t const l = x.n_elem, 
               k = X.n_rows,
              lk = k * l;
#ifdef DDO_CHECKS
  if(y.n_elem != X.n_cols)
    throw std::invalid_argument("x_kron_X_dot_y: invalid y");
#endif
  for(size_t i = 0; i < lk; ++i)
    *(out + i) = 0;
  
  double * const w_end = wk + k;
  for(double * w = wk; w != w_end; ++w)
    *w = 0.;
  
  {
    double const * yp = y.begin();
    for(auto Xp = X.begin(); Xp != X.end(); ++yp){
      double const mult = *yp;
      for(double * w = wk; w != w_end; ++w, ++Xp)
        *w += *Xp * mult;
    }
  }
  
  double * o = out;
  for(auto xp = x.begin(); xp != x.end(); ++xp)    
    for(double * w = wk; w != w_end; ++w, ++o)
      *o = *w * *xp;
}

/**
 * computes (X (x) I) where X is a k x p matrix, I is the l-dimensional
 * identity matrix, x is an l * p vector and stores the result  in the
 * passed pointer which can hold k * l elements.
 */
inline void X_kron_I_dot_x
  (arma::mat const &X, size_t const l, arma::vec const &x,
   double * const out, bool const set_zero){
  size_t const k = X.n_rows,
               p = X.n_cols,
              kl = k * l;
#ifdef DDO_CHECKS
  if(x.n_elem != l * p)
    throw std::invalid_argument("X_kron_I_dot_x: invalid out");
#endif
  for(size_t i = 0; i < kl and set_zero; ++i)
    *(out + i) = 0;

  for(size_t c = 0; c < p; ++c){
    for(size_t r = 0; r < k; ++r){
      double const mult = X.at(r, c);
      double const * const x_end = x.memptr() + c * l + l;
      double * o = out + r * l;
      for(double const * xp = x.memptr() + c * l;
          xp != x_end; ++xp, ++o)
        *o += *xp * mult;
    }
  }
}

/**
 * computes (I (x) X) where X is a k x p matrix, I is the l-dimensional
 * identity matrix, x is an l * p vector and stores the result  in the
 * passed pointer which can hold k * l elements.
 */
inline void I_kron_X_dot_x
  (arma::mat const &X, size_t const l, arma::vec const &x,
   double * const out){
  size_t const k = X.n_rows,
               p = X.n_cols,
              kl = k * l;
#ifdef DDO_CHECKS
  if(x.n_elem != l * p)
    throw std::invalid_argument("I_kron_X_dot_x: invalid out");
#endif
  for(size_t i = 0; i < kl; ++i)
    *(out + i) = 0;

  for(size_t c = 0; c < l; ++c){
    double * const o_begin = out + c * k,
           * const o_end   = o_begin + k;
    double const * Xp = X.memptr(),
                 * xp_end = x.memptr() +  c * p + p;
    for(double const * xp = x.memptr() + c * p; xp != xp_end; ++xp){
      double const mult = *xp;
      for(auto o = o_begin; o != o_end; ++o, ++Xp)
        *o += mult * *Xp;
    }
  }
}

#endif
