/* Given a multinomial variable with K levels with probabilities
 * p_1,\dots,p_K we assume that we observe outcome k if and only if
 * A_k > A_l \forall l \neq k and we assume that
 * (A_1,\dots,A_K)^\top \sim N(\vec\mu, I) where \mu_1 = 0.
 * Below, we provide function to recover the \vec\mu vector given the
 * probabilities of each category. */

#include "multinomial-probit.h"
#include <cmath>
#include <algorithm>
#include <memory>
#include <array>
#include "pnorm.h"
#include "qnorm.h"
#include "psqn-bfgs.h"

namespace multinomial {
/** util function to get the value of the mean vector. */
inline double get_mu(int const idx, double const *mu) noexcept {
  if(idx < 1L)
    return 0;
  return mu[idx - 1L];
}

/** evalutes the log of the normal distribution density. */
inline double eval_log_dnrm(double const x) noexcept {
  // log(2 * pi) / 2
  constexpr double const log2pi_half = 0.918938533204673;
  return -log2pi_half - x * x / 2.;
}

/** evalutes:
 *   \phi(x; \mu_k, 1)\prod_{i \neq k}\Phi(x;\mu_i; 1)
 *  assuming that \mu_1 = 0. */
inline double fn(double const x, double const *mu, int const nvars,
                 int const icase, bool const use_log = false){
  double out(0.);
  // add the log phi term
  out += eval_log_dnrm(x - get_mu(icase, mu));
  // add the log Phi terms
  for(int i = 0; i < nvars; ++i){
    if(i == icase)
      continue;
    out += pnorm_std(x - get_mu(i, mu), 1, 1);
  }

  return use_log ? out : std::exp(out);
}

/** computes the derivative of log fn as a function of x. Returns a 2D
 *  array with the log integrand value and the gradient. */
inline std::array<double, 2L>
log_gr_x(double const x, double const *mu, int const nvars,
         int const icase) noexcept {
  double log_int(0.),
         deriv  (0.);
  // add the log phi term
  {
    double const diff = x - get_mu(icase, mu);
    log_int += eval_log_dnrm(diff);
    deriv -= diff;
  }
  // add the log Phi terms
  for(int i = 0; i < nvars; ++i){
    if(i == icase)
      continue;
    double const diff = x - get_mu(i, mu),
             log_pnrm = pnorm_std(diff, 1, 1);
    log_int += log_pnrm;

    deriv += std::exp(eval_log_dnrm(diff) - log_pnrm);
  }

  return { log_int, deriv };
}

/** computes the second order derivative of log fn as a function of x. */
inline double he_x(double const x, double const *mu, int const nvars,
                   int const icase) noexcept {
  double out(0);
  out -= 1;
  for(int i = 0; i < nvars; ++i){
    if(i == icase)
      continue;
    double const diff = x - get_mu(i, mu),
                    f = pnorm_std(diff, 1, 0);
    double const fp = std::exp(eval_log_dnrm(diff)),
                fpp = -diff * fp;

    out += (f * fpp - fp * fp) / f / f;
  }

  return out;
}

/** computes the derivative of fn as a function of the means (mu)
 *  disregarding the first variable. */
inline double gr(double const x, double const *mu, int const nvars,
                 int const icase, double *gr_val){
  double mult(0.);
  // add the phi term
  {
    double const diff = x - get_mu(icase, mu);
    mult += eval_log_dnrm(diff);
    if(icase > 0L)
      gr_val[icase - 1L] = diff;
  }
  // add the log phi terms
  for(int i = 0; i < nvars; ++i){
    if(i == icase)
      continue;

    double const diff = x - get_mu(i, mu),
             log_pnrm = pnorm_std(diff, 1, 1);
    mult += log_pnrm;

    if(i > 0){
      double const log_dnrm = eval_log_dnrm(diff);
      gr_val[i - 1L] = -std::exp(log_dnrm - log_pnrm);
    }
  }

  // handle other factors
  mult = std::exp(mult);
  for(int i = 0; i < nvars - 1L; ++i)
    gr_val[i] *= mult;

  return mult;
}

/** class used to find the mode for adaptive Gauss-Hermite quadrature. */
class inner_problem final : public PSQN::problem {
  double const * const mu;
  int const icase, nvars;
public:
  inner_problem(double const *mu, int const icase, int const nvars):
  mu(mu), icase(icase), nvars(nvars) { }

  size_t size() const {
    return 1L;
  }

  double func(double const *val){
    return -fn(*val, mu, nvars, icase, true);
  }

  double grad(double const * __restrict__ val,
              double       * __restrict__ gr){
    auto const out = log_gr_x(*val, mu, nvars, icase);
    *gr = -out[1];
    return -out[0];
  }
};

/** finds the mode for adaptive Gauss-Hermite quadrature. */
double solve_inner_problem
(double const *mu, int const icase, int const nvars){
  double out(0);
  inner_problem prob(mu, icase, nvars);
  auto opt_inf = bfgs(prob, &out, 1e-6, 100, 1e-4, .9);

  return opt_inf.info == PSQN::info_code::converged ? out : 0.;
}

/// fixed quadrature nodes and weights
constexpr int const n_nodes = 50L;
constexpr double const ghq_nodes[n_nodes] = { -9.18240695812931, -8.52277103091781, -7.97562236820563, -7.48640942986419, -7.03432350977061, -6.60864797385536, -6.20295251927468, -5.81299467542041, -5.43578608722495, -5.06911758491724, -4.71129366616904, -4.36097316045458, -4.01706817285813, -3.67867706251527, -3.34503831393789, -3.01549776957452, -2.68948470226775, -2.36649390429866, -2.04607196868641, -1.7278065475159, -1.4113177548983, -1.09625112895768, -0.782271729554607, -0.469059056678236, -0.156302546889469, 0.156302546889469, 0.469059056678236, 0.782271729554607, 1.09625112895768, 1.4113177548983, 1.7278065475159, 2.04607196868641, 2.36649390429866, 2.68948470226775, 3.01549776957452, 3.34503831393789, 3.67867706251527, 4.01706817285814, 4.36097316045458, 4.71129366616904, 5.06911758491724, 5.43578608722495, 5.81299467542041, 6.20295251927468, 6.60864797385536, 7.03432350977061, 7.48640942986419, 7.97562236820563, 8.52277103091781, 9.18240695812931 },
                     ghq_weights[n_nodes] = { 1.83379404857338e-37, 1.67380166790777e-32, 1.21524412340452e-28, 2.13765830836018e-25, 1.41709359957341e-22, 4.47098436540802e-20, 7.74238295704335e-18, 8.09426189346537e-16, 5.46594403181544e-14, 2.50665552389957e-12, 8.11187736493018e-11, 1.90904054381192e-09, 3.3467934040214e-08, 4.4570299668178e-07, 4.58168270795559e-06, 3.68401905378071e-05, 0.000234269892109255, 0.00118901178174964, 0.00485326382617194, 0.0160319410684121, 0.0430791591567653, 0.0945489354770862, 0.170032455677164, 0.251130856332002, 0.305085129204399, 0.305085129204399, 0.251130856332002, 0.170032455677164, 0.0945489354770865, 0.0430791591567657, 0.0160319410684121, 0.00485326382617195, 0.00118901178174965, 0.000234269892109255, 3.68401905378071e-05, 4.58168270795546e-06, 4.45702996681779e-07, 3.34679340402137e-08, 1.90904054381188e-09, 8.11187736493008e-11, 2.50665552389963e-12, 5.46594403181558e-14, 8.09426189346515e-16, 7.74238295704327e-18, 4.47098436540803e-20, 1.41709359957345e-22, 2.13765830836021e-25, 1.21524412340452e-28, 1.67380166790769e-32, 1.83379404857343e-37 };
constexpr double const sqrt2 = 1.4142135623731; // sqrt(2)

double eval(double const *mu, int const icase, int const nvars){
  double out(0.);
  double const mode = solve_inner_problem(mu, icase, nvars),
               nhes = -he_x(mode, mu, nvars, icase),
               scal = nhes > 0 ? sqrt2 * std::sqrt(1. / nhes) : sqrt2;

  for(int i = 0; i < n_nodes; ++i){
    int const idx = i % 2 == 0 ? i / 2L : n_nodes - 1L - i / 2L;
    double const x = ghq_nodes[idx],
                xa = scal * x + mode;
    out += ghq_weights[idx] * std::exp(x * x) * fn(xa, mu, nvars, icase);
  }

  return out * scal;
}

/** approximates the derivative of eval with pre-allocated working
 *  memory. */
double eval_gr(double const *mu, double *gr_val, int const icase,
               int const nvars, double *wk){
  double out(0.);
  double const mode = solve_inner_problem(mu, icase, nvars),
               nhes = -he_x(mode, mu, nvars, icase),
               scal = nhes > 0 ? sqrt2 * std::sqrt(1 / nhes) : sqrt2;

  int const nvars_m1 = nvars - 1;
  std::fill(gr_val, gr_val + nvars_m1, 0.);
  std::fill(wk    , wk     + nvars_m1, 0.);

  for(int i = 0; i < n_nodes; ++i){
    int const idx = i % 2 == 0 ? i / 2L : n_nodes - 1L - i / 2L;
    double const x = ghq_nodes[idx],
                xa = scal * x + mode,
                 w = ghq_weights[idx] * std::exp(x * x);
    out += w * gr(xa, mu, nvars, icase, wk);
    for(int j = 0; j < nvars_m1; ++j)
      gr_val[j] += w * wk[j];
  }

  for(int j = 0; j < nvars_m1; ++j)
    gr_val[j] *= scal;
  return out * scal;
}

/** approximates the derivative of eval without pre-allocated working
 *  memory. */
double eval_gr(double const *mu, double *gr_val, int const icase,
               int const nvars){
  std::unique_ptr<double[]> wk_mem(new double[nvars - 1]);
  return eval_gr(mu, gr_val, icase, nvars, wk_mem.get());
}

/** class to minimize the KL divergence. */
class mult_problem final : public PSQN::problem {
  int const nvars;
  double const * const probs;
  std::unique_ptr<double[]> wk_mem =
    std::unique_ptr<double[]>(new double[2 * (nvars - 1)]);

public:
  mult_problem(double const *probs, int const nvars):
  nvars(nvars), probs(probs) { }

  size_t size() const {
    return nvars - 1L;
  }

  double func(double const *val){
    double out(0.);
    for(int i = 0; i < nvars; ++i)
      out -= probs[i] * std::log(eval(val, i, nvars));
    return out;
  }

  double grad(double const * __restrict__ val,
              double       * __restrict__ gr){
    double out(0.);
    std::fill(gr, gr + nvars - 1, 0.);
    double * const g = wk_mem.get(),
           * const w = g + nvars - 1;
    for(int i = 0; i < nvars; ++i){
      double const likelihood = eval_gr(val, g, i, nvars, w);
      out -= probs[i] * std::log(likelihood);
      for(int j = 0; j < nvars - 1; ++j)
        gr[j] -= probs[i] * g[j] / likelihood;
    }
    return out;
  }
};

int find_means(double const *probs, double *means, int const nvars,
               double const rel_eps, int const max_it,
               double const c1, double const c2){
  mult_problem prob(probs, nvars);

  // find the starting values
  {
    std::unique_ptr<double[]> tmp_mem(new double[nvars]);
    double * const w = tmp_mem.get();
    for(int i = 0; i < nvars; ++i)
      w[i] = qnorm_w(probs[i], 0, 1., 1L, 0L);
    for(int i = 1; i < nvars; ++i)
      means[i - 1] = w[i] - w[0];
  }

  auto out = bfgs(prob, means, rel_eps, max_it, c1, c2);
  return static_cast<int>(out.info);
}

} // namespace multinomial
