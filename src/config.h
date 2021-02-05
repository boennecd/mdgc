#ifndef MDGC_CONFIG_H
#define MDGC_CONFIG_H

#include <cstddef>

#ifdef DO_CHECKS
#define MDGC_NOEXCEPT
#else
#define MDGC_NOEXCEPT noexcept
#endif

#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#define MDGC_RESTRICT
#else
#define MDGC_RESTRICT __restrict__
#endif

namespace mdgc {
inline constexpr size_t cacheline_size(){
  return 128L;
}
} // namespace mdgc

#endif
