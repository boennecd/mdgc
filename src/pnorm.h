#ifndef PNORM_H
#define PNORM_H

#include <Rmath.h>

#ifdef __cplusplus
#include <cmath>
using std::isinf;
using std::isnan;
#else
#include <math.h>
#endif

/**
 * evaluates the standard normal CDF after avoiding some checks in the
 * R function.
 */
static inline double pnorm_std(double const x, int lower, int is_log) {
  if(isinf(x) || isnan(x))
    return NAN;

  double p, cp;
  p = x;
  Rf_pnorm_both(x, &p, &cp, lower ? 0 : 1, is_log);
  return lower ? p : cp;
}

/**
 * evaluates the normal CDF after avoiding some checks in the R function.
 */
static inline double pnorm_w(double const x, double const mu, double const sigma,
                             int lower, int is_log) {
  return pnorm_std((x - mu) / sigma, lower, is_log);
}

#endif
