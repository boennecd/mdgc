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
#endif

#endif
