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

## Resubmission
This is a resubmission. I have addressed the issues stated below.

> Please add \value to .Rd files regarding exported methods and explain
> the functions results in the documentation. Please write about the
> structure of the output (class) and also what the output means. (If a
> function does not return a value, please document that too, e.g.
> \value{No return value, called for side effects} or similar)
> Missing Rd-tags:
>       get_mdgc_log_ml.Rd: \value
>       get_mdgc.Rd: \value
>       mdgc_fit.Rd: \value
>       mdgc_start_value.Rd: \value
>       mdgc.Rd: \value

\value has been added to all the listed .Rd files.

> \dontrun{} should only be used if the example really cannot be executed
> (e.g. because of missing additional software, missing API keys, ...) by
> the user. That's why wrapping examples in \dontrun{} adds the comment
> ("# Not run:") as a warning for the user.
> Does not seem necessary.
> 
> Please unwrap the examples if they are executable in < 5 sec, or replace
> \dontrun{} with \donttest{}.

I have replaced \dontrun{} with \donttest{}.

> Please always add all authors, contributors and copyright holders in the
> Authors@R field with the appropriate roles.
> e.g.: Ross Ihaka, Genz and Bretz
> Please explain in the submission comments what you did about this issue.

All the third party code is from the mvtnorm package or the stats package. 
src/qnorm.c is from src/nmath/qnorm.c in R. src/new-mvt.h and mvt.f is from the 
mvtnorm package and parts of it has been re-written in C++ along with some 
changes. The code that I am using is originally by Alan Genz and Frank Bretz.
Like the mnormt package that also uses code by these authors, I have added
them to the Authors field.

I do not think I should use the copyright field in the description file as all 
copyright holders are authors of the code I have used or changed. Writing R 
Extensions states that:

> An optional ‘Copyright’ field can be used where the copyright holder(s) are 
> not the authors. If necessary, this can refer to an installed file: the 
> convention is to use file inst/COPYRIGHTS.

I apologize if I am wrong.
