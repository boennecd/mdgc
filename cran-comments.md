## Test environments
* Ubuntu 20.04 LTS with gcc 10.1.0
  R version 4.2.1
* Ubuntu 20.04 LTS with gcc 10.1.0
  R version 4.2.1 using Valgrind
* GitHub actions on windows-latest (release), macOS-latest (release), 
  ubuntu-20.04 (release), ubuntu-20.04 (old-release), and ubuntu-20.04 (devel)
* win-builder (devel, oldrelease, and release)
* `rhub::check_for_cran()`
* `rhub::check(platform = c("macos-highsierra-release", "macos-highsierra-release-cran"))`
  
## R CMD check results
There were no WARNINGs or ERRORs.

There is a NOTE about the package size in some cases.
