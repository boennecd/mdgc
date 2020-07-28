#ifdef __cplusplus
#include <vector>
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
}

extern "C"
{
#endif

/**
 * samples from a standard normal distribution.
 */
double rngnorm_wrapper ();
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
double rngunif_wrapper ();

#ifdef __cplusplus
}
#endif
