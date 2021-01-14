#ifndef MULTINOMIAL_PROBIT_H
#define MULTINOMIAL_PROBIT_H

namespace multinomial {

/**
 * approximates
 *     \int \phi(x; \mu_k, \sigma_k)\prod_{i \neq k}\Phi(x;\mu_i, \sigma_i) dx
 *
 * with \sigma_i being one if k \neq 1 or -> 0^+ otherwise.
 *
 * @param mu pointer to the mean vector of the nvars - 1 last elements. The
 * first mean element is assumed to be zero.
 * @param icase integer in 0,...,nvars - 1 with the index k.
 * @param nvars number of categories.
 */
double eval(double const *mu, int const icase, int const nvars);

/**
 * approximates the derivative of eval with pre-allocated working memory.
 *
 * @param mu pointer to the mean vector of the nvars - 1 last elements. The
 * first mean element is assumed to be zero.
 * @param gr_val gradient wrt mu on return.
 * @param icase integer in 0,...,nvars - 1 with the index k.
 * @param nvars number of categories.
 * @param wk working memory of size nvars - 1.
 */
double eval_gr(double const *mu, double *gr_val, int const icase,
               int const nvars, double *wk);

/** finds the mean values corresponding to the passed probabilities.
 *  The methods minimizes the KL diverence using the psqn method. The
 *  returned integer is the info code from the psqn package.
 *
 * @param probs pointer to the nvars probabilities summing to one.
 * @param means vector with space for the nvars - 1 potentially non-zero
 * means.
 * @param nvars number of categories.
 * @param rel_eps relative convergence threshold.
 * @param max_it maximum number of iterations.
 * @param c1,c2 thresholds in the Wolfe conditions.
 */
int find_means(double const *probs, double *means, int const nvars,
               double const rel_eps = 3.000214e-13,
               int const max_it = 100, double const c1 = .0001,
               double const c2 = .9);

} // namespace multinomial

#endif
