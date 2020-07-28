#ifndef PNORM_H
#define PNORM_H

#ifdef __cplusplus
extern "C" {
#endif
#include <Rmath.h>
#include <math.h>

/**
 * evaluates the standard normal CDF after avoiding some checks in the
 * R function.
 */
inline double pnorm_std(double const x, int lower, int is_log){
  if(__builtin_expect((isinf(x) || isnan(x)), 0))
    return NAN;

  double p, cp;
  p = x;
  Rf_pnorm_both(x, &p, &cp, lower ? 0 : 1, is_log);
  return lower ? p : cp;
}

/**
 * evaluates the normal CDF after avoiding some checks in the R function.
 */
inline double pnorm(double const x, double const mu, double const sigma,
                    int lower, int is_log){
  return pnorm_std((x - mu) / sigma, lower, is_log);
}

/*
   log CDF approximation from

   xup <- 0
   dput(off <- pnorm(xup, log = TRUE))
   xs <- seq(-8, xup, length.out = 1000)
   xt <- xs - xup
   fit <- lm(pnorm(xs, log.p = TRUE) ~ poly(xt, raw = TRUE, degree = 3) - 1,
   offset = rep(off, length(xs)))
   xhat <- seq(-9, xup, length.out = 1000)
   matplot(xhat, cbind(pnorm(xhat, log.p = TRUE),
   predict(fit, newdata = data.frame(xt = xhat - xup))),
   type = "l")
   xhat <- seq(xup - 1, xup, length.out = 1000)
   matplot(xhat, cbind(pnorm(xhat, log.p = TRUE),
   predict(fit, newdata = data.frame(xt = xhat - xup))),
   type = "l")
   max(abs(fit$residuals))
   dput(coef(fit))

   xs <- seq(xup, 8, length.out = 1000)
   xt <- xs - xup
   fit <- lm(pnorm(xs, log.p = TRUE) ~ poly(xt, raw = TRUE, degree = 5) - 1,
   offset = rep(off, length(xs)))
   xhat <- seq(xup, 8, length.out = 1000)
   matplot(xhat, cbind(pnorm(xhat, log.p = TRUE),
   predict(fit, newdata = data.frame(xt = xhat - xup))),
   type = "l")
   max(abs(fit$residuals))
   dput(coef(fit))
*/

/**
 * Cheap and crude approximation of the log of the normal CDF.
 */
inline double pnorm_std_aprx(double const x){
  static double const xup = 0,
                    inter = -0.693147180559945;
  if(x > -8 && x < 8){
    if(x >= xup) {
      double const xd = x - xup;
      double out = inter,
              xp = xd;
      out += 0.811401689963717 * xp;
      xp *= xd;
      out -= 0.36464495501633 * xp;
      xp *= xd;
      out += 0.0784688043503829 * xp;
      xp *= xd;
      out -= 0.00809951215678632 * xp;
      xp *= xd;
      out += 0.000321901567031593 * xp;
      return out;

    } else {
      double const xd = x - xup;
      double out = inter,
              xp = xd;
      out += 0.7633873031817 * xp;
      xp *= xd;
      out -= 0.365089939925538 * xp;
      xp *= xd;
      out += 0.0146225647154304 * xp;
      xp *= xd;
      out += 0.000646820246088735 * xp;
      return out;

    }
  }

  return pnorm_std(x, 1L, 1L);
}

#ifdef __cplusplus
}
#endif

#endif
