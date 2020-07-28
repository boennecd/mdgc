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

#include <RcppArmadillo.h>

#endif
