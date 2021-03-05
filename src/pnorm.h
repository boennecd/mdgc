#ifndef PNORM_H
#define PNORM_H

#include <Rmath.h>

#ifdef __cplusplus
#include <cmath>
#include <limits>
#define PNORM_ISNAN std::isnan
#define PNORM_ISINF std::isinf
#define PNORM_NAN std::numeric_limits<double>::quiet_NaN();
#else
#include <math.h>
#define PNORM_ISNAN isnan
#define PNORM_ISINF isinf
#define PNORM_NAN R_NaN
#endif

/**
 * evaluates the standard normal CDF after avoiding some checks in the
 * R function.
 */
static inline double pnorm_std(double const x, int lower, int is_log) {
  if(PNORM_ISINF(x) || PNORM_ISNAN(x))
    return PNORM_NAN;

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

#undef PNORM_ISNAN
#undef PNORM_ISINF
#undef PNORM_NAN

#endif
