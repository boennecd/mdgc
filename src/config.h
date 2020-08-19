#ifndef MDGC_CONFIG_H
#define MDGC_CONFIG_H

#ifdef DO_CHECKS
#define MDGC_NOEXCEPT
#else
#define MDGC_NOEXCEPT noexcept
#endif

#endif
