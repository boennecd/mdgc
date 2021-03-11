## Test environments
* Ubuntu 18.04 LTS with gcc 10.1.0
  R version 3.6.3
* Ubuntu 18.04 LTS with gcc 10.1.0
  R version 3.6.3 with LTO checks
* Ubuntu 18.04 LTS with gcc 10.1.0
  R version 3.6.3 with valgrind
* Ubuntu 18.04 LTS with clang-6.0.0
  R devel 2021-02-03 r79933 with ASAN and UBSAN
* Ubuntu 18.04 LTS with gcc 10.1.0
  R devel 2021-02-05 r79941 with ASAN and UBSAN
* Github actions on windows-latest (release), macOS-latest (release), 
  ubuntu-20.04 (release), and ubuntu-20.04 (devel)
* win-builder (devel, oldrelease, and release)
* `rhub::check_for_cran()`
* `rhub::check(platform = c("fedora-clang-devel", "solaris-x86-patched", "macos-highsierra-release-cran"))`

## R CMD check results
There were no WARNINGs or ERRORs except with the ASAN and UBSAN checks with 
clang-6.0.0 which I think are false positives (see below).

There is a NOTE about the package size in some cases.

I have attempted to closely follow the description of the Solaris setup on CRAN. 
In particular, I build GCC 5.2.0 as the the archived version of OpenCSW do not 
work with the math.h file that is updated on 27/02/16 (see 
http://gcc.1065356.n8.nabble.com/the-mystery-of-math-h-in-lib-gcc-triple-name-gcc-version-include-fixed-td1524637.html). However, I still fail to reproduce any of the ERRORs 
on CRAN. 

My configurations for R are: 

./configure \
  CC="/opt/gcc-5.2.0/bin/gcc" \
  CXX="/opt/gcc-5.2.0/bin/g++" \
  CPPFLAGS="-I/opt/gcc-5.2.0/include -I/opt/csw/include -I/usr/local/include" \
  FC="/opt/gcc-5.2.0/bin/gfortran" \
  CFLAGS="-O2" \
  FFLAGS="-O2" \
  CXXFLAGS="-O2" \
  LDFLAGS="-L/opt/gcc-5.2.0/lib -L/usr/local/lib -L/opt/csw/lib" \
  --with-internal-tzcode \
  R_LD_LIBRARY_PATH="/opt/gcc-5.2.0/lib:/usr/local/lib:/opt/csw/lib:/opt/developerstudio12.6/lib:/usr/openwin/lib" \
  --prefix=/opt/R-gcc-5-2
  
where /opt/gcc-5.2.0 contains the GCC-5.2.0 I build. As I have failed again to 
reproduce the errors on Solaris, I have removed the examples and tests on 
Solaris.

The ASAN and UBSAN checks with clang-6.0.0 yields a false positive I think. I 
get the following:	

> ==29904==ERROR: AddressSanitizer: stack-use-after-scope on address 0x7ffe86462b40 at pc 0x7fc45c2ea579 bp 0x7ffe86462950 sp 0x7ffe86462948	
> WRITE of size 8 at 0x7ffe86462b40 thread T0	
>     #0 0x7fc45c2ea578 in .omp_outlined._debug__.22 /tmp/ci-pyMjb0sjid/mdgc.Rcheck/00_pkg_src/mdgc/src/cpp-to-R.cpp:696:32	
>     #1 0x7fc45c2ea578 in .omp_outlined..23 /tmp/ci-pyMjb0sjid/mdgc.Rcheck/00_pkg_src/mdgc/src/cpp-to-R.cpp:696	
>     #2 0x7fc45bd98452 in __kmp_invoke_microtask (/usr/lib/x86_64-linux-gnu/libomp.so.5+0x7c452)	
>     #3 0x7fc45bd52c4a in __kmp_fork_call (/usr/lib/x86_64-linux-gnu/libomp.so.5+0x36c4a)	
>     #4 0x7fc45bd467be in __kmpc_fork_call (/usr/lib/x86_64-linux-gnu/libomp.so.5+0x2a7be)	
>     #5 0x7fc45c2e5533 in impute(arma::Mat<double> const&, arma::Mat<double> const&, arma::Mat<int> const&, arma::Mat<double> const&, arma::Col<double> const&, arma::Mat<double> const&, Rcpp::Vector<19, Rcpp::PreserveStorage>, Rcpp::Vector<19, Rcpp::PreserveStorage>, double, double, unsigned int, Rcpp::Vector<19, Rcpp::PreserveStorage>, Rcpp::Vector<16, Rcpp::PreserveStorage>, int, bool, int, bool) /tmp/ci-pyMjb0sjid/mdgc.Rcheck/00_pkg_src/mdgc/src/cpp-to-R.cpp:690:9	
>     #6 0x7fc45c2afb9a in _mdgc_impute /tmp/ci-pyMjb0sjid/mdgc.Rcheck/00_pkg_src/mdgc/src/RcppExports.cpp:137:34	
>     ...	
> 	
> Address 0x7ffe86462b40 is located in stack of thread T0 at offset 480 in frame	
>     #0 0x7fc45c2e829f in .omp_outlined..23 /tmp/ci-pyMjb0sjid/mdgc.Rcheck/00_pkg_src/mdgc/src/cpp-to-R.cpp:696	
> 	
>   This frame has 28 object(s):	
>     [32, 56) '.kmpc_loc.addr.i.i'	
>     [96, 104) 'ref.tmp.i.i'	
>     [128, 408) 'agg.tmp56.i'	
>     [480, 488) 'i.i' (line 696) <== Memory access at offset 480 is inside this variable	
>     [512, 520) '.omp.lb.i' (line 690)	
>     [544, 552) '.omp.ub.i' (line 690)	
>     [576, 584) '.omp.stride.i' (line 690)	
>     [608, 612) '.omp.is_last.i' (line 690)	
>     [624, 736) 'w_idx_int4.i' (line 690)	
>     [768, 880) 'w_idx_obs5.i' (line 690)	
>     [912, 1088) 'w_obs_val6.i' (line 690)	
>     [1152, 1328) 'w_upper7.i' (line 690)	
>     [1392, 1568) 'w_lower8.i' (line 690)	
>     [1632, 1808) 'mu9.i' (line 690)	
>     [1872, 1984) 'w_idx_cat_obs10.i' (line 690)	
>     [2016, 2128) 'w_idx_cat_not_obs11.i' (line 690)	
>     [2160, 2208) 'known_objs12.i' (line 690)	
>     [2240, 2352) 'a_idx_int.i' (line 690)	
>     [2384, 2496) 'a_idx_obs.i' (line 690)	
>     [2528, 2704) 'a_obs_val.i' (line 690)	
>     [2768, 2944) 'a_upper.i' (line 690)	
>     [3008, 3184) 'a_lower.i' (line 690)	
>     [3248, 3424) 'mu_use.i' (line 690)	
>     [3488, 3512) 'type_i.i' (line 690)	
>     [3552, 3664) 'a_idx_cat_obs.i' (line 690)	
>     [3696, 3808) 'a_idx_cat_not_obs.i' (line 690)	
>     [3840, 4016) 'D.i' (line 690)	
>     [4080, 4104) '.kmpc_loc.addr.i' (line 924)	
> HINT: this may be a false positive if your program uses some custom stack unwind mechanism or swapcontext	
>       (longjmp and C++ exceptions *are* supported)	

This seems very similar to the bug reported here: https://github.com/google/sanitizers/issues/994#issue-356441940.	
In particular, notice the stack trace in the example in the reported issue:	

>     #0 0x5174eb in .omp_outlined._debug__ /home/mgigg/Code/git/martyngigg/sandbox/cpp/asan/omp/omp.cpp:5:36	
>     #1 0x7f3bdf780452 in __kmp_invoke_microtask (/usr/lib/x86_64-linux-gnu/libomp.so.5+0x7c452)	
>     #2 0x7f3bdf73a1b6  (/usr/lib/x86_64-linux-gnu/libomp.so.5+0x361b6)	
>     #3 0x7f3bdf73b2b5 in __kmp_fork_call (/usr/lib/x86_64-linux-gnu/libomp.so.5+0x372b5)	
>     #4 0x7f3bdf72e7be in __kmpc_fork_call (/usr/lib/x86_64-linux-gnu/libomp.so.5+0x2a7be)	
>     #5 0x516f88 in main /home/mgigg/Code/git/martyngigg/sandbox/cpp/asan/omp/omp.cpp:4:11	
>     #6 0x7f3bdf11cb96 in __libc_start_main /build/glibc-OTsEL5/glibc-2.27/csu/../csu/libc-start.c:310	

It is also clang-6.0.0 and the looping variable is blamed: 	

>  This frame has 6 object(s):	
>     [32, 36) 'i' (line 5) <== Memory access at offset 32 is inside this variable	
>     [48, 52) '.omp.lb' (line 4)	
>     [64, 68) '.omp.ub' (line 4)	
>     [80, 84) '.omp.stride' (line 4)	
>     [96, 100) '.omp.is_last' (line 4)	
>     [112, 136) '.kmpc_loc.addr'
