## Test environments
* Ubuntu 18.04 LTS with gcc 10.1.0
  R version 3.6.3
* Ubuntu 18.04 LTS with gcc 10.1.0
  R version 3.6.3 with valgrind
* Ubuntu 18.04 LTS with gcc 10.1.0
  R devel 2020-11-30 r79529 with LTO checks
* Ubuntu 18.04 LTS with clang-6.0.0
  R devel 2021-02-03 r79933 with ASAN and UBSAN
* Ubuntu 18.04 LTS with gcc 10.1.0
  R devel 2021-02-05 r79941 with ASAN and UBSAN
* Github actions on windows-latest (release), macOS-latest (release), 
  ubuntu-20.04 (release), and ubuntu-20.04 (devel)
* win-builder (devel and release)
* `rhub::check_for_cran()`
* `rhub::check(platform = c("solaris-x86-patched", "macos-highsierra-release-cran"))`

## R CMD check results
There were no WARNINGs or ERRORs.

There is a NOTE about the package size in some cases.
