#include "threat-safe-random.h"
#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <Rmath.h>

#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/beta_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/uniform_real.hpp>

namespace {
using rng_type = boost::mt19937;

static std::unique_ptr<rng_type[]> generators;

rng_type& get_generator(){
#ifdef _OPENMP
  int const me = omp_get_thread_num();
#else
  int const me = 0L;
#endif

#ifdef DO_CHECKS
  if(!(generators.get() + me))
    throw std::runtime_error("get_generator: no generator for this thread");
#endif
  return generators[me];
}
}

namespace parallelrng {
void set_rng_seeds(std::vector<unsigned> const &seeds){
  std::size_t const n = seeds.size();
  generators.reset(new rng_type[n]);

  unsigned i = 0L;
  for(auto s : seeds)
    generators[i++] = rng_type(s);
}

void set_rng_seeds(unsigned const n_threads){
#ifdef DO_CHECKS
  if(n_threads <= 0L)
    throw std::invalid_argument("set_rng_seeds: invalid 'n_threads'");
#endif

  std::vector<unsigned> seeds;
  seeds.reserve(n_threads);

  double const max = 10000000.;
  for(unsigned i = 0; i < n_threads; ++i)
    seeds.emplace_back(static_cast<unsigned>(unif_rand() * max + .5));

  set_rng_seeds(seeds);
}
}

extern "C"
{
double rngnorm_wrapper(){
  typedef boost::normal_distribution<> normal_generic;
  static normal_generic const n01(0.0, 1.0);

  boost::variate_generator<rng_type&, normal_generic>
    rng_normal(get_generator(), n01);

  return rng_normal();
}

double rngbeta_wrapper(double const a, double const b){
  typedef boost::random::beta_distribution<> beta_genric;
  beta_genric betaD(a, b);
  boost::variate_generator<rng_type&, beta_genric>
    rng_beta(get_generator(), betaD);

  return rng_beta();
}

double rnggamma_wrapper(double const a){
  typedef boost::random::gamma_distribution<> gamma_generic;
  gamma_generic gammaD(a);

  boost::variate_generator<rng_type&, gamma_generic>
    rng_gamma(get_generator(), gammaD);

  return rng_gamma();
}

double rngunif_wrapper(){
  typedef boost::uniform_real<> unif_generic;

  unif_generic unifD(0., 1.);
  boost::variate_generator<rng_type&, unif_generic>
    rng_unif(get_generator(), unifD);

  return rng_unif();
}
}
