#ifndef MDGC_CONFIG_H
#define MDGC_CONFIG_H

#ifdef DO_CHECKS
#define MDGC_NOEXCEPT
#else
#define MDGC_NOEXCEPT noexcept
#endif

inline constexpr size_t cacheline_size(){
  return 128L;
}

#endif
