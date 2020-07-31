#ifndef PROF_H
#define PROF_H
#include <string>
#include <atomic>

class profiler {
#ifdef DO_PROF
  static std::atomic<bool> running_profiler;
#endif

public:
  profiler(const std::string&);

#ifdef DO_PROF
  ~profiler();
#endif
};

#endif
