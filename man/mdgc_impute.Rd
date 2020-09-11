% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mdgc.R
\name{mdgc_impute}
\alias{mdgc_impute}
\title{Impute Missing Values Given Correlation Matrix}
\usage{
mdgc_impute(
  object,
  vcov,
  rel_eps = 0.001,
  maxit = 10000L,
  abs_eps = -1,
  n_threads = 1L,
  do_reorder = TRUE,
  minvls = 1000L
)
}
\arguments{
\item{object}{returned object from \code{\link{get_mdgc}}.}

\item{vcov}{correlation matrix to condition on in the imputation.}

\item{rel_eps}{relative error for each term.}

\item{maxit}{maximum number of iteration.}

\item{abs_eps}{absolute convergence threshold for each term.}

\item{n_threads}{number of threads to use.}

\item{do_reorder}{logical for whether to use heuristic variable reordering.}

\item{minvls}{minimum number of samples.}
}
\description{
Impute Missing Values Given Correlation Matrix
}