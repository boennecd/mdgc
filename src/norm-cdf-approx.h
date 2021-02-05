#ifndef NORM_CDF_APPROX_H
#define NORM_CDF_APPROX_H

#include "config.h"

/** normal CDF approximations using monotone cubic interpolation. */
double pnorm_approx(double const x) MDGC_NOEXCEPT;

#endif
