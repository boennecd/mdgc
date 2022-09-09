#ifndef ARMA_WRAP_H
#define ARMA_WRAP_H

#define _USE_MATH_DEFINES
#include <cmath>

#ifndef ARMA_DONT_USE_OPENMP
#define ARMA_DONT_USE_OPENMP 1
#endif

#ifndef DO_CHECKS
#define ARMA_NO_DEBUG
#endif

#ifdef ARMA_WARN_LEVEL
#undef ARMA_WARN_LEVEL
#endif
#define ARMA_WARN_LEVEL 1

#include <RcppArmadillo.h>

#endif
