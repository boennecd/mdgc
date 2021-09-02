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
/** evaluates the log of the normal distribution density. */
inline double eval_log_dnrm(double const x) MDGC_NOEXCEPT {
  // log(2 * pi) / 2
  constexpr double log2pi_half = 0.918938533204673;
  return -log2pi_half - x * x / 2.;
}

/// fixed quadrature nodes and weights
/*
 gw <- pracma::gaussLaguerre(100L)
 keep <- gw$w > 1e-100
 dput(gw$x[keep])
 dput(gw$w[keep])
 sum(keep)
*/

constexpr int n_nodes = 56L;
constexpr double glq_nodes[n_nodes] = { 0.0143861469954186, 0.0758036120233572, 0.186314102057189, 0.345969180991432, 0.554810937580917, 0.812891284115669, 1.12027383500754, 1.47703432992383, 1.8832608263424, 2.33905384964603, 2.84452654275536, 3.39980482744571, 4.00502758175865, 4.66034683556891, 5.36592798558512, 6.12195003080402, 6.92860582937617, 7.78610237786251, 8.69466111392217, 9.65451824355508, 10.6659250941217, 11.7291484944722, 12.844471183641, 14.0121922496943, 15.2326276004667, 16.506110468082, 17.8329919493264, 19.2136415841361, 20.6484479746683, 22.1378194476567, 23.6821847630024, 25.281993871834, 26.9377187275743, 28.6498541538913, 30.4189187737909, 32.2454560045207, 34.1300351234215, 36.07325241038, 38.0757323731071, 40.1381290621155, 42.2611274829848, 44.4454451143118, 46.6918335406515, 49.0010802107724, 51.374010332704, 53.8114889183557, 56.3144229919617, 58.8837639782821, 61.5205102883961, 64.2257101231016, 67.0004645164193, 69.8459306445584, 72.7633254289746, 75.7539294659399, 78.8190913194115, 81.9602322190601 },
                     glq_weights[n_nodes] = { 0.0363926058834001, 0.0796767462129539, 0.112115103342487, 0.130356612975145, 0.134043339728462, 0.125407090780663, 0.108314112097261, 0.0870966384699589, 0.0655510093123107, 0.0463401335826442, 0.0308463086276817, 0.0193678281139787, 0.0114854423601796, 0.00643895100161029, 0.00341497998969265, 0.0017143197401822, 0.000814871591587831, 0.000366854836599486, 0.000156452074178106, 6.32108705288851e-05, 2.41957522651891e-05, 8.77430976375552e-06, 3.01426748600091e-06, 9.80833589934502e-07, 3.0226387435323e-07, 8.82005839529586e-08, 2.43642585620072e-08, 6.3697113739018e-09, 1.57560032045964e-09, 3.68632920134617e-10, 8.15479892422471e-11, 1.70506255682607e-11, 3.36821417086667e-12, 6.28352495536662e-13, 1.1064980159833e-13, 1.83835017545486e-14, 2.8801150572306e-15, 4.25261289733722e-16, 5.91443246725226e-17, 7.74308910250001e-18, 9.53624944907477e-19, 1.10409451698693e-19, 1.20084697043049e-20, 1.22600700749045e-21, 1.17402171464669e-22, 1.05359406460872e-23, 8.85323176847421e-25, 6.95919877548308e-26, 5.1123698154656e-27, 3.50627214881706e-28, 2.24264627269388e-29, 1.33620942344572e-30, 7.40742895757554e-32, 3.81586013713635e-33, 1.82420019768289e-34, 8.08163233606434e-36 };

double eval(double const *mu, int const icase, int const nvars){
  if(icase < 1L){
    double out(0.);
    for(int i = 0; i < nvars - 1L; ++i)
      out += pnorm_std(-mu[i], 1, 1);

    return std::exp(out);
  }

  double out(0.);
  for(int i = 0; i < n_nodes; ++i){
    double const x = glq_nodes[i];
    double new_term = x + eval_log_dnrm(x - mu[icase - 1]);
    for(int j = 1; j < nvars; ++j){
      if(j == icase)
        continue;
      new_term += pnorm_std(x - mu[j - 1], 1, 1);
    }

    out += std::exp(new_term) * glq_weights[i];
  }

  return out;
}

double eval_gr(double const *mu, double *gr_val, int const icase,
               int const nvars, double *wk){
  int const nvars_m1 = nvars - 1;
  if(icase < 1L){
    double out(0.);
    for(int i = 0; i < nvars_m1; ++i){
      double const log_pnorm = pnorm_std(-mu[i], 1, 1);
      out += log_pnorm;
      gr_val[i] = eval_log_dnrm(-mu[i]) - log_pnorm;
    }

    for(int i = 0; i < nvars_m1; ++i)
      gr_val[i] = -std::exp(gr_val[i] + out);

    return std::exp(out);
  }

  double out(0.);
  std::fill(gr_val, gr_val + nvars_m1, 0.);
  for(int i = 0; i < n_nodes; ++i){
    double const x = glq_nodes[i];

    double diff_icase = x - mu[icase - 1],
           dnrm_icase = eval_log_dnrm(diff_icase),
             new_term = x + dnrm_icase;
    wk[icase - 1] = diff_icase;

    for(int j = 1; j < nvars; ++j){
      if(j == icase)
        continue;
      double const muj = mu[j - 1],
                  diff = x - muj,
             log_pnorm = pnorm_std(diff, 1, 1);
      new_term += log_pnorm;
      wk[j - 1] = -std::exp(eval_log_dnrm(diff) - log_pnorm);
    }

    new_term = std::exp(new_term) * glq_weights[i];
    for(int i = 0; i < nvars_m1; ++i)
      gr_val[i] += wk[i] * new_term;

    out += new_term;
  }

  return out;
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

  double grad(double const * MDGC_RESTRICT val,
              double       * MDGC_RESTRICT gr){
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
