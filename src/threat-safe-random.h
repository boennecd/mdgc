#ifndef THREAT_SAFE_RNG_H
#define THREAT_SAFE_RNG_H

#ifdef __cplusplus
#include <vector>
#include <stdint.h>

/// forward declare class
namespace boost {
namespace random {
template<typename UIntType, std::size_t w, std::size_t n, std::size_t m,
         std::size_t r, UIntType a, std::size_t u, UIntType d, std::size_t s,
         UIntType b, std::size_t t, UIntType c, std::size_t l, UIntType f>
class mersenne_twister_engine;
} // namespace random
} // namespace boost

namespace parallelrng {
/**
 * set the seeds for up the number of element of the vector. Must be called
 * before calling any of the subsequent methods in this header.
 *
 * @param seeds Vector with seeds to use.
 */
void set_rng_seeds(std::vector<unsigned> const &seeds);
/**
 * set a given number of random seeds.
 *
 * @param n_threads Number of seeds to set.
 */
void set_rng_seeds(unsigned const n_threads);

typedef boost::random::mersenne_twister_engine
  <uint32_t, 32, 624, 397, 31, 0x9908b0df, 11, 0xffffffff, 7, 0x9d2c5680, 15, 0xefc60000, 18, 1812433253> rng_type;

/** class to draw uniform (0, 1) variables */
class unif_drawer {
  rng_type &gen;

public:
  unif_drawer(rng_type &gen): gen(gen) { }

  /// draws a uniform (0, 1) variable
  double operator()();
};

// get thread specific generator
unif_drawer get_unif_drawer();
}

extern "C"
{
#endif

/**
 * samples from a standard normal distribution.
 */
double rngnorm_wrapper (void);
/**
 * samples from a beta distribution.
 */
double rngbeta_wrapper (double const, double const);
/**
 * samples from a gamma distribution.
 */
double rnggamma_wrapper(double const);
/**
 * samples from a uniform distribution.
 */
double rngunif_wrapper (void);

#ifdef __cplusplus
}
#endif

#endif
