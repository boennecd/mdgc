#ifndef QNORM_H
#define QNORM_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * slightly modified version of R's qnorm to avoid some checks.
 */
double qnorm_w(double const p, double const mu, double const sigma,
               int const lower_tail, int const log_p);
#ifdef __cplusplus
}

#include "config.h"
#include <limits>

/***
 * 7 digits precision approximation from:
 *    ALGORITHM AS241 APPL. STATIST. (1988) VOL. 37, NO. 3
 */
inline double qnorm_aprx(double const p) MDGC_NOEXCEPT {
  constexpr double       split1 = .425,
                         split2 = 5,
                         const1 = .180625,
                         const2 = 1.6,

                             AO = 3.3871327179,
                             A1 = 5.0434271938e1,
                             A2 = 1.5929113202e2,
                             A3 = 5.9109374720e1,
                             B1 = 1.7895169469e1,
                             B2 = 7.8757757664e1,
                             B3 = 6.7187563600e1,

                             C0 = 1.4234372777,
                             C1 = 2.7568153900,
                             C2 = 1.3067284816,
                             C3 = 1.7023821103e-1,
                             D1 = 7.3700164250e-1,
                             D2 = 1.2021132975e-1,

                             E0 = 6.6579051150,
                             E1 = 3.0812263860,
                             E2 = 4.2868294337e-1,
                             E3 = 1.7337203997e-2,
                             F1 = 2.4197894225e-1,
                             F2 = 1.2258202635e-2;

  // if(p > 0 and p < 1){
    // ok branch
    double const q = p - .5;
    if(std::fabs(q) < split1){
      double const r = const1 - q * q;
      return q * (((A3 * r + A2) * r + A1) * r + AO) /
        (((B3 * r + B2) * r + B1) * r + 1);
    }

    double r = q < 0 ? p : 0.5 - p + 0.5;
    r = std::sqrt(-std::log(r));

    double out;
    if(r <= split2){
      r -= const2;
      out = (((C3 * r + C2) * r + C1) * r + C0) / ((D2 * r + D1) * r + 1);
    } else {
      r -= split2;
      out =  (((E3 * r + E2) * r + E1) * r + E0) / ((F2 * r + F1) * r + 1);

    }

    return q < 0 ? -out : out;
  // }

  // // issues. p <= 0 or p >= 1!
  // if(p < 0 || p > 1)
  //   return std::numeric_limits<double>::quiet_NaN();
  // else if(p == 0)
  //   return -std::numeric_limits<double>::infinity();
  //
  // // p == 1
  // return std::numeric_limits<double>::infinity();
}

#endif

#endif
