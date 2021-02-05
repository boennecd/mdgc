#' mdgc: Missing Data imputation using Gaussian Copulas
#'
#' @description
#'
#' The mdgc package is used to estimate Gaussian Copula models for mixed data
#' types (continuous, binary, ordinal, and multinomial) that can be used for
#' imputations. The main function is the \code{\link{mdgc}} function. The rest
#' of the functions in the package give the user access to lower level
#' functions.
#'
#' Examples are provided at \url{https://github.com/boennecd/mdgc}.
#'
#' @references
#' Christoffersen, B., Clements, M., Humphreys, K., & Kjellström, H. (2021).
#' \emph{Asymptotically Exact and Fast Gaussian Copula Models for Imputation of Mixed Data Types}.
#' \url{https://arxiv.org/abs/2102.02642}.
#'
#' @aliases mdgc-package
"_PACKAGE"

#' @importFrom Rcpp sourceCpp
#' @useDynLib mdgc, .registration = TRUE
NULL
