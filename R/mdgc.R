#' Get mdgc Object
#' @description
#' Creates a mdgc object which is needed for estimation of the
#' covariance matrix and the mean vector and to perform imputation.
#'
#' @details
#' It is important to use appropriate classes for the \code{\link{data.frame}}
#' columns:
#'
#' \itemize{
#'   \item{Continuous variables: }{should be \code{\link{numeric}}s.}
#'   \item{Binary variables: }{should be \code{\link{logical}}s.}
#'   \item{Multinomial variables: }{should be \code{\link{factor}}s.}
#'   \item{Ordinal variables: }{should be \code{\link{ordered}}.}
#' }
#'
#' @return
#' An object of class \code{mdgc}. It has the following elements:
#' \item{lower,upper,code,multinomial,idx_non_zero_mean}{arguments to pass to
#' \code{\link{get_mdgc_log_ml}}.}
#' \item{margs}{functions to get \code{lower} and \code{upper} bounds for each
#' column of \code{dat}.}
#' \item{reals,bins,ords}{indices of continuous, binary, and ordinal variables,
#' respectively.}
#' \item{truth}{the numeric version of \code{dat}.}
#' \item{means}{starting values for the non-zero mean terms
#' (see e.g. \code{\link{mdgc_fit}}).}
#'
#' @param dat \code{\link{data.frame}} with continuous, multinomial, ordinal, and binary
#' variables.
#' @importFrom stats na.omit
#' @importFrom utils head
#' @export
#' @importFrom stats qnorm ecdf quantile
#'
#' @examples
#' # randomly mask data
#' set.seed(11)
#' masked_data <- iris
#' masked_data[matrix(runif(prod(dim(iris))) < .10, NROW(iris))] <- NA
#'
#' # use the functions in the package
#' library(mdgc)
#' obj <- get_mdgc(masked_data)
#' class(obj)
#'
#' @seealso
#' \code{\link{get_mdgc_log_ml}}, \code{\link{mdgc_start_value}}
get_mdgc <- function(dat){
  # checks
  stopifnot(is.data.frame(dat), NROW(dat) > 0L)

  # get the types
  reals <- which(sapply(dat, is.numeric))
  bins  <- which(sapply(dat, is.logical))
  ords  <- which(sapply(dat, is.ordered))
  mults <- which(sapply(dat, is.factor))
  mults <- setdiff(mults, ords)

  nc <- NCOL(dat)
  stopifnot(
    length(reals) + length(bins) + length(ords) + length(mults) == nc)

  # estimate the marginals
  margs <- lapply(1:nc, function(j){
    x <- na.omit(dat[[j]])
    if(length(x) < 1L)
      stop(sprintf("No data for column %d", j))

    if(is.numeric(x)){
      cdf_func <- ecdf(x)
      n <- length(x)
      scal <- n / (n + 1)
      mi <- min(x)
      return(function(x)
        cbind(NA_real_, qnorm(scal * cdf_func(pmax(x, mi)))))

    } else if(is.logical(x) || is.ordered(x)){
      phat <- table(x) / length(x)
      stopifnot(all(phat > 0), length(phat) > 1)
      lb <- c(-Inf, qnorm(cumsum(head(phat, -1))))
      ub <- c(qnorm(cumsum(head(phat, -1))), Inf)

      is_bin <- length(lb) < 3
      if(is_bin){
        mu <- qnorm(phat[2])
        lb[2] <- ub[1] <- 0
      }

      lvls <- if(is.logical(x))
        c(FALSE, TRUE) else levels(x)

      out <- if(is.logical(x))
        function(x){
          idx <- as.integer(x) + 1L
          cbind(lb[idx], ub[idx])
        }
      else
        function(x){
          idx <- as.integer(x)
          cbind(lb[idx], ub[idx])
        }
      attr(out, "levels") <- lvls
      attr(out, "borders") <- c(-Inf, ub)
      if(is_bin)
        attr(out, "mean") <- mu
      return(out)

    } else if(is.factor(x)){
      phat <- table(x) / length(x)
      mu <- multinomial_find_means(phat)
      if(attr(mu, "info-code") != 0)
        warning("multinomial_find_means had non-zero convergence code")

      out <- function(x){
        lv <- c(vapply(!is.na(x), function(x)
          rep(if(x) -Inf else NA_real_, length(phat)),
          numeric(length(phat))))
        ub <- c(vapply(!is.na(x), function(x)
          rep(if(x) 0 else NA_real_, length(phat)),
          numeric(length(phat))))

        n <- length(x)
        list(lower = matrix(lv, n, byrow = TRUE),
             upper = matrix(ub, n, byrow = TRUE))
      }
      attr(out, "levels") <- levels(x)
      attr(out, "means") <- mu

      return(out)
    }

    stop("type not implemented")
  })

  # get upper and lower bounds or points for the latent variables
  Z <- mapply(function(x, func) func(x), x = dat, func = margs,
              SIMPLIFY = FALSE)
  lower <- t(do.call(
    cbind, lapply(Z, function(x){
      if(is.matrix(x)) x[, 1] else x$lower
    })))
  upper <- t(do.call(
    cbind, lapply(Z, function(x){
      if(is.matrix(x)) x[, 2] else x$upper
    })))

  # assign codes for missingness pattern
  code <- matrix(0L, NROW(lower), NCOL(lower))
  is_na <- is.na(lower) + is.na(upper)
  code[is_na == 0L] <- 2L
  code[is_na == 2L] <- 1L

  # construct matrix for multinomial outcomes
  if(length(mults)){
    n_lvls <- sapply(dat[mults], function(x) length(levels(x)))
    n_vars <- sapply(dat, function(x)
      if(is.factor(x) && !is.ordered(x)) length(levels(x)) else 1L)
    idx <- c(0L, cumsum(head(n_vars, -1L)))

    multinomial <- lapply(1:NROW(dat), function(i)
      rbind(as.integer(dat[i, mults]),
            n_lvls,
            idx[mults]))

  } else {
    idx <- seq_along(margs) - 1
    multinomial <- replicate(NROW(dat), matrix(0L, 0L, 0L))
  }

  # get the true values (needed for imputation)
  truth <- t(sapply(dat, as.numeric))
  dimnames(truth) <- rev(dimnames(dat))

  # get the non-zero means
  means <- sapply(
    margs, function(x){
      out <- attr(x, "mean")
      if(!is.null(out))
        return(out)
      attr(x, "means")
  })
  idx_non_zero_mean <- mapply(function(meas, idx){
    if(length(meas) == 0L)
      return(NULL)
    else if(length(meas) == 1L)
      # a binary variable
      return(idx)

    # a multinomial variable
    idx + seq_along(meas)
  }, meas = means, idx = idx, SIMPLIFY = FALSE)

  means <- unlist(means)
  if(is.null(means))
    means <- numeric()
  idx_non_zero_mean <- unlist(idx_non_zero_mean)
  if(is.null(idx_non_zero_mean))
    idx_non_zero_mean <- integer()

  structure(list(
    lower = lower, upper = upper, code = code, margs = margs,
    reals = reals, bins = bins, ords = ords, truth = truth,
    multinomial = multinomial, means = means,
    idx_non_zero_mean = idx_non_zero_mean), class = "mdgc")
}

#' Get Pointer to C++ Object to Approximate the Log Marginal Likelihood
#' @description
#' Creates a C++ object which is needed to approximate the log marginal
#' likelihood. The object cannot be saved.
#'
#' @param object mdgc object from \code{\link{get_mdgc}} or a
#' \code{\link{data.frame}} to pass to \code{\link{get_mdgc}}. Ignored by
#' the default method.
#' @param lower [# variables]x[# observations] matrix with lower bounds
#' for each variable on the normal scale.
#' @param upper [# variables]x[# observations] matrix with upper bounds
#' for each variable on the normal scale.
#' @param code [# variables]x[# observations] matrix integer code for the
#' each variable on the normal scale. Zero implies an observed value (the
#' value in \code{upper}), one implies a missing value, and two implies an
#' interval.
#' @param multinomial \code{\link{list}} where each element is
#' 3x[# multinomial variables] \code{\link{matrix}} with
#' multinomial outcomes. The first index is the outcome as an integer code,
#' the second index is the number of categories, and the third index is the
#' index of each multinomial variable (this is zero-based).
#' @param idx_non_zero_mean indices for non-zero mean variables. Indices
#' should be sorted.
#' @param ... used to pass arguments to S3 methods.
#'
#' @details
#' Indices are zero-based except the outcome index for multinomial
#' variables.
#'
#' \code{idx_non_zero_mean} indices with terms with \code{code} equal to zero
#' (observed values) are ignored.
#'
#' @examples
#' # randomly mask data
#' set.seed(11)
#' masked_data <- iris
#' masked_data[matrix(runif(prod(dim(iris))) < .10, NROW(iris))] <- NA
#'
#' # use the functions in the package
#' library(mdgc)
#' obj <- get_mdgc(masked_data)
#' ptr <- get_mdgc_log_ml(obj)
#'
#' @return
#' A \code{Rcpp::XPtr} to pass to e.g. \code{\link{mdgc_log_ml}}.
#'
#' @seealso
#' \code{\link{mdgc_fit}}, \code{\link{mdgc_log_ml}}
#'
#' @export
get_mdgc_log_ml <- function(object, ...)
  UseMethod("get_mdgc_log_ml")

#' @rdname get_mdgc_log_ml
#' @export
get_mdgc_log_ml.mdgc <- function(object, ...)
  get_mdgc_log_ml.default(lower = object$lower, upper = object$upper,
                          code = object$code,
                          multinomial = object$multinomial,
                          idx_non_zero_mean = object$idx_non_zero_mean)

#' @rdname get_mdgc_log_ml
#' @export
get_mdgc_log_ml.data.frame <- function(object, ...)
  get_mdgc_log_ml(get_mdgc(object))

#' @rdname get_mdgc_log_ml
#' @export
get_mdgc_log_ml.default <- function(object, lower, upper, code, multinomial,
                                    idx_non_zero_mean, ...){
  # checks
  di <- dim(lower)
  stopifnot(
    di[1] > 0L, di[2] > 0L, length(di) == 2L,
    all(di == dim(upper)), all(di == dim(code)),
    is.numeric(lower), is.numeric(upper), is.integer(code),
    all(is.finite(code)), all(range(code) %in% 0:2),
    is.list(multinomial),
    all(sapply(multinomial, function(x) is.matrix(x) && is.integer(x))),
    all(dim(multinomial[[1L]]) == sapply(multinomial, dim)))

  keep <- colSums(is.na(lower) & is.na(upper)) < NROW(upper)
  out <- get_log_lm_terms_cpp(lower = lower[, keep, drop = FALSE],
                              upper = upper[, keep, drop = FALSE],
                              code =  code [, keep, drop = FALSE],
                              multinomial = multinomial[keep],
                              idx_non_zero_mean = idx_non_zero_mean)
  attr(out, "nobs") <- NCOL(upper)
  attr(out, "nvars") <- NROW(lower)
  attr(out, "multinomial") <- multinomial
  attr(out, "idx_non_zero_mean") <- idx_non_zero_mean
  out
}

#' Evaluate the Log Marginal Likelihood and Its Derivatives
#'
#' @description
#' Approximates the log marginal likelihood and the derivatives using
#' randomized quasi-Monte Carlo. The method uses a generalization of the Fortran
#' code by Genz and Bretz (2002).
#'
#' Mean terms for observed continuous variables are always assumed to be
#' zero.
#'
#' The returned log marginal likelihood is not a proper log marginal likelihood
#' if the \code{ptr} object is constructed from a mdgc object from
#' \code{\link{get_mdgc}} as it does not include the log of the determinants of
#' the Jacobians for the transformation of the continuous variables.
#'
#' @return
#' A numeric vector with a single element with the log marginal likelihood
#' approximation. Two attributes are added if \code{comp_derivs} is
#' \code{TRUE}: \code{"grad_vcov"} for the derivative approximation with
#' respect to \code{vcov} and \code{"grad_mea"} for the derivative
#' approximation with respect to \code{mea}.
#'
#' @seealso
#' \code{\link{mdgc_fit}}
#'
#' @param ptr object returned by \code{\link{get_mdgc_log_ml}}.
#' @param vcov covariance matrix.
#' @param mea vector with non-zero mean entries.
#' @param rel_eps relative error for each marginal likelihood factor.
#' @param n_threads number of threads to use.
#' @param comp_derivs logical for whether to approximate the gradient.
#' @param indices integer vector with which terms (observations) to include.
#' Must be zero-based. \code{NULL} yields all observations.
#' @param do_reorder logical for whether to use a heuristic variable
#' reordering. \code{TRUE} is likely the best option.
#' @param maxpts maximum number of samples to draw for each marginal likelihood
#' term.
#' @param abs_eps absolute convergence threshold for each marginal likelihood
#' factor.
#' @param minvls minimum number of samples.
#' @param use_aprx logical for whether to use an approximation of
#' \code{\link{pnorm}} and \code{\link{qnorm}}. This may yield a
#' noticeable reduction in the computation time.
#'
#' @references
#' Genz, A., & Bretz, F. (2002). \emph{Comparison of Methods for the
#' Computation of Multivariate t Probabilities}.
#' Journal of Computational and Graphical Statistics.
#'
#' Genz, A., & Bretz, F. (2008).
#' \emph{Computation of Multivariate Normal and t Probabilities}.
#' Springer-Verlag, Heidelberg.
#'
#' @examples
#' # randomly mask data
#' set.seed(11)
#' masked_data <- iris
#' masked_data[matrix(runif(prod(dim(iris))) < .10, NROW(iris))] <- NA
#'
#' # use the functions in the package
#' library(mdgc)
#' obj <- get_mdgc(masked_data)
#' ptr <- get_mdgc_log_ml(obj)
#' start_vals <- mdgc_start_value(obj)
#' mdgc_log_ml(ptr, start_vals, obj$means)
#' mdgc_log_ml(ptr, start_vals, obj$means, use_aprx = TRUE)
#' mdgc_log_ml(ptr, start_vals, obj$means, use_aprx = TRUE, comp_derivs = TRUE)
#'
#' @export
mdgc_log_ml <- function(ptr, vcov, mea, rel_eps = 1e-2, n_threads = 1L,
                        comp_derivs = FALSE, indices = NULL,
                        do_reorder = TRUE, maxpts = 100000L,
                        abs_eps = -1., minvls = 100L,
                        use_aprx = FALSE){
  nvars <- attr(ptr, "nvars")
  nobs <- attr(ptr, "nobs")
  idx_non_zero_mean <- attr(ptr, "idx_non_zero_mean")
  stopifnot(!is.null(nvars), !is.null(nobs), !is.null(idx_non_zero_mean))
  if(is.null(indices))
    indices <- 0:(nobs - 1L)

  stopifnot(
    all(indices >= 0L & indices < nobs),
    all(dim(vcov) == c(nvars, nvars)),
    is.numeric(rel_eps), length(rel_eps) == 1L, is.finite(rel_eps),
    is.integer(n_threads), length(n_threads) == 1L, n_threads > 0L,
    is.logical(comp_derivs), length(comp_derivs) == 1L, !is.na(comp_derivs),
    is.logical(do_reorder), length(do_reorder) == 1L, !is.na(do_reorder),
    is.integer(maxpts), length(maxpts) == 1L, maxpts > 0L,
    is.numeric(abs_eps), length(abs_eps) == 1L, is.finite(abs_eps),
    is.integer(minvls), length(minvls) == 1L, minvls <= maxpts,
    minvls >= 0L,
    is.logical(use_aprx), length(use_aprx) == 1L,
    is.finite(use_aprx),
    length(mea) == length(idx_non_zero_mean), is.numeric(mea),
    all(is.finite(mea)))

  .mdgc_log_ml(
    ptr = ptr, vcov = vcov, mu = mea, indices = indices, maxpts = maxpts,
    abs_eps = abs_eps, rel_eps = rel_eps, n_threads = n_threads,
    comp_derivs = comp_derivs, do_reorder = do_reorder, minvls = minvls,
    use_aprx = use_aprx)
}

.mdgc_log_ml <- function(
  ptr, vcov, mu, rel_eps, n_threads, comp_derivs, indices, do_reorder,
  maxpts, abs_eps, minvls, use_aprx)
  eval_log_lm_terms(
    ptr = ptr, vcov = vcov, mu = mu, indices = indices, maxpts = maxpts,
    abs_eps = abs_eps, rel_eps = rel_eps, n_threads = n_threads,
    comp_derivs = comp_derivs, do_reorder = do_reorder, minvls = minvls,
    use_aprx = use_aprx)

#' Get Starting Value for the Covariance Matrix Using a Heuristic
#'
#' @description
#' Uses a heuristic to get starting values for the covariance matrix. These
#' can be passed e.g. to \code{\link{mdgc_fit}}.
#'
#' @inheritParams get_mdgc_log_ml
#' @inheritParams mdgc_fit
#' @param mea vector with non-zero mean entries.
#' @param n_threads number of threads to use.
#' @param object mdgc object from \code{\link{get_mdgc}}. Ignored by
#' the default method.
#'
#' @return
#' The starting value for the covariance matrix.
#'
#' @examples
#' # randomly mask data
#' set.seed(11)
#' masked_data <- iris
#' masked_data[matrix(runif(prod(dim(iris))) < .10, NROW(iris))] <- NA
#'
#' # use the functions in the package
#' library(mdgc)
#' obj <- get_mdgc(masked_data)
#' ptr <- get_mdgc_log_ml(obj)
#' start_vals <- mdgc_start_value(obj)
#' start_vals # starting value for the covariance matrix
#'
#' @importFrom stats cov cov2cor
#' @export
mdgc_start_value <- function(object, ...)
  UseMethod("mdgc_start_value")

#' @rdname mdgc_start_value
#' @method mdgc_start_value mdgc
#' @export
mdgc_start_value.mdgc <- function(object, ...)
  mdgc_start_value.default(
    lower = object$lower, upper = object$upper, code = object$code,
    multinomial = object$multinomial,
    idx_non_zero_mean = object$idx_non_zero_mean, mea = object$means, ...)

#' @rdname mdgc_start_value
#' @method mdgc_start_value default
#' @export
mdgc_start_value.default <- function(object, lower, upper, code,
                                     multinomial, idx_non_zero_mean, mea,
                                     n_threads = 1L, ...){
  stopifnot(length(idx_non_zero_mean) == length(mea),
            all(is.finite(mea)), is.numeric(mea),
            !anyDuplicated(idx_non_zero_mean),
            all(idx_non_zero_mean < NROW(lower)))

  # get the mean vector to use
  mea_use <- numeric(NROW(lower))
  if(length(idx_non_zero_mean) > 0)
    mea_use[idx_non_zero_mean + 1L] <- mea

  Z <- get_z_hat(lower - mea_use, upper - mea_use,
                 code, multinomial = multinomial, n_threads = n_threads)
  out <- cov(t(Z), use = "pairwise.complete.obs")

  first_mult <- multinomial[[1L]]
  any_cate <- length(first_mult) > 0
  if(any_cate)
    for(i in first_mult[3, ] + 1)
      out[i, i] <- 1

  # handle non-postive definite estimates
  eg <- eigen(out, symmetric = TRUE)
  eps <- 1e-2 * max(eg$values)
  if(any(is_small <- eg$values < eps)){
    eg$values[is_small] <- eps
    out <- tcrossprod(eg$vectors * rep(eg$values, each = NROW(Z)), eg$vectors)
  }

  # rescale
  out <- cov2cor(out)
  if(any_cate)
    for(i in first_mult[3, ] + 1)
      out[i, i] <- 0
  out
}

#' Estimate the Model Parameters
#'
#' @description
#' Estimates the covariance matrix and the non-zero mean terms.
#' The \code{lr} parameter and the \code{batch_size} parameter are likely
#' data dependent.
#' Convergence should be monitored e.g. by using \code{verbose = TRUE}
#' with \code{method = "svrg"}.
#'
#' See the README at \url{https://github.com/boennecd/mdgc} for examples.
#'
#' @seealso
#' \code{\link{mdgc_log_ml}}, \code{\link{mdgc_start_value}},
#' \code{\link{mdgc_impute}}.
#'
#' @param ptr returned object from \code{\link{get_mdgc_log_ml}}.
#' @param vcov,mea starting value for the covariance matrix and the
#' non-zero mean entries.
#' @param batch_size number of observations in each batch.
#' @param lr learning rate.
#' @param decay the learning rate used by SVRG is given by \code{lr * decay^iteration_number}.
#' @param method estimation method to use. Can be \code{"svrg"}, \code{"adam"},
#' or \code{"aug_Lagran"}.
#' @param maxit maximum number of iteration.
#' @param seed fixed seed to use. Use \code{NULL} if the seed should not be
#' fixed.
#' @param epsilon,beta_1,beta_2 ADAM parameters.
#' @param conv_crit relative convergence threshold.
#' @param verbose logical for whether to print output during the estimation.
#' @param mu starting value for the penalty in the augmented Lagrangian
#' method.
#' @param lambda starting values for the Lagrange multiplier estimates.
#' \code{NULL} yields a default.
#' @inheritParams mdgc_log_ml
#'
#' @return
#' An \code{\link{list}} with the following elements:
#' \item{result}{\code{\link{list}} with two elements: \code{vcov} is the
#' estimated covariance matrix and \code{mea} is the estimated non-zero mean
#' terms.}
#' \item{estimates}{If present, the estimated parameters after each iteration.}
#' \item{fun_vals}{If present, the output of \code{\link{mdgc_log_ml}} after
#' each iteration.}
#' \item{mu,lambda}{If present, the \code{mu} and \code{lambda} values at the
#' end.}
#'
#' The elements that may be present depending on the chosen \code{method}.
#'
#' @references
#' Kingma, D.P., & Ba, J. (2015). \emph{Adam: A Method for Stochastic Optimization}. abs/1412.6980.
#'
#' Johnson, R., & Zhang, T. (2013). \emph{Accelerating stochastic gradient descent using predictive variance reduction}. In Advances in neural information processing systems.
#'
#' @examples
#' \donttest{
#' # randomly mask data
#' set.seed(11)
#' masked_data <- iris
#' masked_data[matrix(runif(prod(dim(iris))) < .10, NROW(iris))] <- NA
#'
#' # use the functions in the package
#' library(mdgc)
#' obj <- get_mdgc(masked_data)
#' ptr <- get_mdgc_log_ml(obj)
#' start_vals <- mdgc_start_value(obj)
#'
#' fit <- mdgc_fit(ptr, start_vals, obj$means, rel_eps = 1e-2, maxpts = 10000L,
#'                 minvls = 1000L, use_aprx = TRUE, batch_size = 100L, lr = .001,
#'                 maxit = 100L, n_threads = 2L)
#' fit$result$vcov
#' fit$result$mea
#' }
#'
#' @importFrom stats optim
#' @importFrom utils tail
#' @export
mdgc_fit <- function(ptr, vcov, mea, lr = 1e-3, rel_eps = 1e-3,
                     maxit = 25L, batch_size = NULL,
                     method = c("svrg", "adam", "aug_Lagran"), seed = 1L,
                     epsilon = 1e-8,
                     beta_1 = .9, beta_2 = .999, n_threads = 1L,
                     do_reorder = TRUE, abs_eps = -1., maxpts = 10000L,
                     minvls = 100L, verbose = FALSE, decay = .98,
                     conv_crit = 1e-6, use_aprx = FALSE, mu = 1,
                     lambda = NULL){
  #####
  # checks
  nvars <- attr(ptr, "nvars")
  nobs <- attr(ptr, "nobs")
  idx_non_zero_mean <- attr(ptr, "idx_non_zero_mean")
  stopifnot(!is.null(nvars), !is.null(nobs), !is.null(idx_non_zero_mean))

  if(is.null(batch_size))
    batch_size <- as.integer(max(min(10L, nobs), nobs / 20))

  method <- method[1L]
  stopifnot(
    all(dim(vcov) == c(nvars, nvars)), is.numeric(vcov),
    is.numeric(lr), length(lr) == 1L, lr > 0,
    is.integer(maxit), length(maxit) == 1L, maxit > 0L,
    is.integer(batch_size), length(batch_size) == 1L, batch_size > 0L,
    batch_size <= nobs,
    is.numeric(rel_eps), length(rel_eps) == 1L, rel_eps >= 0.,
    is.character(method), method %in% c("svrg", "adam", "aug_Lagran"),
    is.integer(seed), length(seed) == 1L, is.finite(seed) || is.null(seed),
    is.integer(n_threads), length(n_threads) == 1L, n_threads > 0L,
    is.numeric(epsilon), length(epsilon) == 1L, epsilon >= 0.,
    is.numeric(beta_1), length(beta_1) == 1L, beta_1 >= 0.,
    is.numeric(beta_2), length(beta_2) == 1L, beta_2 >= 0.,
    is.numeric(abs_eps), length(abs_eps) == 1L, is.finite(abs_eps),
    is.logical(do_reorder), length(do_reorder) == 1L, !is.na(do_reorder),
    is.integer(maxpts), length(maxpts) == 1L, maxpts > 0L,
    is.integer(minvls), length(minvls) == 1L, minvls <= maxpts,
    minvls >= 0L,
    is.logical(verbose), length(verbose) == 1L, !is.na(verbose),
    is.numeric(decay), length(decay) == 1L, decay > 0,
    is.numeric(conv_crit), length(conv_crit) == 1L, is.finite(conv_crit),
    is.logical(use_aprx), length(use_aprx) == 1L,
    is.finite(use_aprx),
    is.numeric(mu), length(mu) == 1L, mu > 0,
    is.null(lambda) || is.numeric(lambda),
    length(mea) == length(idx_non_zero_mean),
    all(is.finite(mea)), is.numeric(mea))

  multinomial <- attr(ptr, "multinomial")
  any_multinomial <- length(multinomial[[1L]]) > 0

  #####
  # assign functions to use

  # we work with a lower dimensional covariance matrix when there are
  # multinomial variables. Thus, we make a few adjustments in this case
  # these are the indices for the free and fixed rows and columns (except
  # for possible scaling constraints)
  fixed_dim <- if(any_multinomial)
    multinomial[[1L]][3, ] + 1 else NULL
  free_dims <- if(any_multinomial)
    setdiff(1:NROW(vcov), fixed_dim) else 1:NROW(vcov)

  # returns the free part of the covariance matrix
  get_free_vcov <- function(x)
    x[free_dims, free_dims]

  # the reverse
  get_free_vcov_inv <- function(x){
    out <- matrix(0., nvars, nvars)
    out[free_dims, free_dims] <- x
    out
  }

  # use the function
  vcov <- get_free_vcov(vcov)

  # indices used to apply a matrix product with a get_commutation matrix
  com_vec <- get_commutation_vec(NROW(vcov), NROW(vcov), FALSE)

  # these are the rows and columns that are subject to a scaling constrain
  is_scale_constrained <- 1:NCOL(vcov)
  if(any_multinomial)
    for(i in 1:NCOL(multinomial[[1L]]))
      is_scale_constrained <- setdiff(
        is_scale_constrained,
        multinomial[[1L]][3, i] - i + 3:multinomial[[1L]][2, i])

  # this function rescales the rows and columns that are subject to a scaling
  # constrain
  rescale_vcov <- function(x){
    if(length(is_scale_constrained) == 0)
      # nothing to do
      return(x)

    sds <- sqrt(diag(x))
    if(length(is_scale_constrained) != NCOL(x))
      # there are some which are not zero
      sds[-is_scale_constrained] <- 1

    # re-scale and return
    sds <- 1/sds
    sds * x * rep(sds, each = NROW(x))
  }

  # indices of a diagonal entries with in a vector with only the lower
  # diagonal
  idx_diag <- c(1L, 1L + cumsum(NCOL(vcov):2))

  # computes the approximate log marginal likelihood.
  #
  # Args:
  #   par_vcov: log-Cholesky decomposition.
  #   par_mea: non-zero mean terms
  #   comp_derivs: logical for whether to approximate the gradient.
  #   indices: indices to use.
  #   maxpts: maximum number of samples to draw.
  par_fn <- function(par_vcov, par_mea, comp_derivs = FALSE, indices,
                     maxpts){
    if(!is.null(seed))
      set.seed(seed)
    Arg <- .get_lchol_inv(par_vcov)
    Arg <- get_free_vcov_inv(Arg)

    res <- .mdgc_log_ml(
      ptr = ptr, Arg, mu = par_mea, comp_derivs = comp_derivs,
      indices = indices, n_threads = n_threads, rel_eps = rel_eps,
      do_reorder = do_reorder, abs_eps = abs_eps, maxpts = maxpts,
      minvls = minvls, use_aprx = use_aprx)
    log_ml <- c(res)
    if(comp_derivs){
      # handle dervatives wrt vcov
      gr_vcov <- attr(res, "grad_vcov")
      gr_vcov <- get_free_vcov(gr_vcov)
      tmp <- matrix(0, NROW(vcov), NROW(vcov))
      tmp[lower.tri(tmp, TRUE)] <- par_vcov
      diag(tmp) <- exp(diag(tmp))
      gr_vcov <- gr_vcov[com_vec] + c(gr_vcov)
      gr_vcov <- x_dot_X_kron_I(x = gr_vcov, X = tmp, l = NROW(vcov))
      gr_vcov <- gr_vcov[, lower.tri(tmp, TRUE)]
      gr_vcov[idx_diag] <- gr_vcov[idx_diag] * diag(tmp)

      attr(log_ml, "grad_vcov") <- gr_vcov

      # handle derivatives wrt par_mea
      attr(log_ml, "grad_mea") <- drop(attr(res, "grad_mea"))

    }

    log_ml
  }

  maxpts_base <- exp(log(minvls / maxpts) / (maxit - 1))
  maxpts_scale <- maxpts_base^(maxit:1)
  maxpts_scale <- maxpts_scale / max(maxpts_scale)
  maxpts_use <- pmax(minvls, as.integer(floor(maxpts * maxpts_scale)))

  if(method == "adam")
    return(adam(
      par_fn = par_fn, nobs = nobs, val_vcov = .get_lchol(vcov),
      val_mea = mea,
      batch_size = batch_size, maxit = maxit, seed = seed,
      epsilon = epsilon, alpha = lr, beta_1 = beta_1, beta_2 = beta_2,
      maxpts = maxpts_use, rescale_vcov = rescale_vcov,
      get_free_vcov_inv = get_free_vcov_inv))
  else if(method == "svrg")
    return(svrg(
      par_fn = par_fn, nobs = nobs, val_vcov = .get_lchol(vcov),
      val_mea = mea,
      batch_size = batch_size, maxit = maxit, seed = seed, lr = lr,
      verbose = verbose, maxpts = maxpts_use, decay = decay,
      conv_crit = conv_crit, rel_eps = rel_eps, rescale_vcov = rescale_vcov,
      get_free_vcov_inv = get_free_vcov_inv))
  else if(method != "aug_Lagran")
    stop(sprintf("Method '%s' is not implemented", method))

  # construct function to evaluate the constraints and the derivatives
  # first assign the indices for the contraints
  constraints_idx <- matrix(rep(is_scale_constrained, 2), ncol = 2,
                            dimnames = list(NULL, c("row", "col")))
  constraints_idx <- constraints_idx - 1
  # then the wanted value
  constraints_val <- rep(1, NROW(constraints_idx))

  # then the functions
  get_cons <- function(par){
    par_exp <- par
    par_exp[idx_diag] <- exp(par[idx_diag])
    conts <- lower_tri_inner(x = par_exp, idx = constraints_idx,
                             jacob = FALSE, rhs = numeric())
    conts - constraints_val
  }
  c_func <- function(par, derivs = FALSE, lambda, mu){
    mu_scaled <- mu * nobs / 2
    conts <- get_cons(par)

    if(!derivs)
      return(-drop(lambda %*% conts) + mu_scaled * sum(conts^2))

    rhs <- -lambda + 2 * mu_scaled * conts
    par_exp <- par
    par_exp[idx_diag] <- exp(par[idx_diag])
    out <- lower_tri_inner(x = par_exp, idx = constraints_idx,
                           jacob = TRUE, rhs = rhs)
    out[idx_diag] <- out[idx_diag] * exp(par[idx_diag])
    out
  }

  # functions to return the parameters for the mean and the covariance
  # matrix
  n_mea <- length(idx_non_zero_mean)
  get_par_vcov <- function(par)
    if(n_mea > 0)
      head(par, -n_mea) else par
  get_par_mea <- function(par)
    if(n_mea > 0)
      tail(par,  n_mea) else numeric()

  # assign function for the augmented problem
  indices <- 0:(nobs - 1L)
  is_ok <- function(par){
    par <- get_par_vcov(par)
    egs <- try(eigen(.get_lchol_inv(par))$values, silent = TRUE)
    if(inherits(egs, "try-error"))
      return(FALSE)
    !(any(egs < 0) || max(egs) * .Machine$double.eps * 10 > min(egs))
  }
  fn <- function(par, mu, lambda){
    if(!is_ok(par))
      return(NA_real_)
    par_vcov <- get_par_vcov(par)
    -par_fn(par_vcov = par_vcov, par_mea = get_par_mea(par),
            comp_derivs = FALSE, indices = indices,
            maxpts = maxpts) + c_func(
              par = par_vcov, derivs = FALSE, lambda = lambda, mu = mu)
  }
  gr <- function(par, mu, lambda){
    if(!is_ok(par))
      stop("invalid 'par' in 'gr'")

    par_vcov <- get_par_vcov(par)
    t1 <- par_fn(
      par_vcov = par_vcov, par_mea = get_par_mea(par),
      comp_derivs = TRUE, indices = indices,
      maxpts = maxpts)
    t2 <- c_func(par = par_vcov, derivs = TRUE, lambda = lambda, mu = mu)

    c(-attr(t1, "grad_vcov") + t2, -attr(t1, "grad_mea"))
  }

  # perform the augmented Lagrangian method
  if(is.null(lambda))
    lambda <- numeric(length(constraints_val))
  stopifnot(length(lambda) == length(constraints_val),
            all(is.finite(lambda)))
  par <- c(.get_lchol(vcov), mea)
  const_norm <- Inf
  for(i in 1:maxit){
    if(verbose)
      cat(sprintf("\nSolving inner problem in iteration %d\n", i))
    opt_res <- optim(par, fn, gr, method = "BFGS", mu = mu, lambda = lambda,
                     control = list(trace = verbose, reltol = conv_crit,
                                    REPORT = 1, maxit = maxit))
    par <- opt_res$par

    # check constraints
    const <- get_cons(get_par_vcov(par))
    all_ok <- all(abs(const) < 1e-4)
    if(all_ok)
      break

    # update lambda, mu and par. Then repeat
    lambda <- lambda - mu * nobs * const
    new_norm <- norm(t(const), "F")
    if(verbose)
      cat(sprintf("Norm of constraint violation: %16.8f\n",
                  new_norm))
    if(new_norm > const_norm / 2)
      mu <- mu * 10
    else
      mu <- mu * 1.5
    const_norm <- new_norm
  }

  vcov_res <- rescale_vcov(.get_lchol_inv(get_par_vcov(par)))
  list(result = list(
    vcov = get_free_vcov_inv(vcov_res),
    mea = get_par_mea(par)), mu = mu, lambda = lambda)
}

#####
# performs stochastic gradient descent instead (using ADAM).
#
# Args:
#   par_fn: function to evaluate the log marginal likelihood.
#   nobs: number of observation.
#   val_vcov,val_mea: starting value for the covariance matrix and the
#                     mean.
#   batch_size: number of observations in each batch.
#   maxit: maximum number of iteration.
#   seed: seed to use.
#   epsilon, alpha, beta_1, beta_2: ADAM parameters.
#   maxpts: maximum number of samples to draw.
#   rescale_vcov: function to rescale the covariance matrix.
#   get_free_vcov_inv: function to get the full covariance matrix which includes
#   fixed parameters.
adam <- function(par_fn, nobs, val_vcov, val_mea, batch_size, maxit = 10L,
                 seed = 1L, epsilon = 1e-8, alpha = .001, beta_1 = .9,
                 beta_2 = .999, maxpts, rescale_vcov,
                 get_free_vcov_inv){
  n_blocks <- max(1L, ceiling(nobs / batch_size))
  get_blocks <- function(){
    indices <- sample.int(nobs, replace = FALSE) - 1L
    blocks <- tapply(indices, (seq_along(indices) - 1L) %% n_blocks,
                     identity, simplify = FALSE)
  }

  # assign function to get the mean and covariance matrix parameters
  n_mea <- length(val_mea)
  get_par_vcov <- function(par)
    if(n_mea > 0)
      head(par, -n_mea) else par
  get_par_mea <- function(par)
    if(n_mea > 0)
      tail(par,  n_mea) else numeric()

  # assign wrapper for par_fn
  par_fn_wrap <- function(val, comp_derivs, ...){
    out <- par_fn(par_vcov = get_par_vcov(val), par_mea = get_par_mea(val),
                  comp_derivs = comp_derivs, ...)
    if(comp_derivs)
      attr(out, "grad") <- c(attr(out, "grad_vcov"), attr(out, "grad_mea"))
    out
  }

  # assign parameter vector and other needed quantites
  val <- c(val_vcov, val_mea)
  is_vcov <- seq_along(val_vcov)

  n_par <- length(val)
  m <- v <- numeric(n_par)
  fun_vals <- numeric(maxit)
  estimates <- matrix(NA_real_, n_par, maxit)
  i <- -1L

  for(k in 1:maxit){
    blocks <- get_blocks()
    for(ii in 1:n_blocks){
      i <- i + 1L
      idx_b <- (i %% n_blocks) + 1L
      m_old <- m
      v_old <- v
      res <- par_fn_wrap(val, comp_derivs = TRUE, indices = blocks[[idx_b]],
                         maxpts[k])
      fun_vals[(i %/% n_blocks) + 1L] <-
        fun_vals[(i %/% n_blocks) + 1L] + c(res)

      gr <- attr(res, "grad")

      m <- beta_1 * m_old + (1 - beta_1) * gr
      v <- beta_2 * v_old + (1 - beta_2) * gr^2

      m_hat <- m / (1 - beta_1^(i + 1))
      v_hat <- v / (1 - beta_2^(i + 1))

      val <- val + alpha * m_hat / (sqrt(v_hat) + epsilon)

      # alter the covariance matrix part
      val[is_vcov] <- .get_lchol(rescale_vcov(.get_lchol_inv(val[is_vcov])))
    }

    estimates[, k] <- val
  }

  list(result = list(
    vcov = get_free_vcov_inv(.get_lchol_inv(get_par_vcov(val))),
    mea = get_par_mea(val)), estimates = estimates)
}

#####
# performs stochastic gradient descent instead (using SVRG).
#
# Args:
#   par_fn: function to evaluate the log marginal likelihood.
#   nobs: number of observation.
#   val_vcov,val_mea: starting value for the covariance matrix and the
#                     mean.
#   batch_size: number of observations in each batch.
#   n_threads: number of threads to use.
#   maxit: maximum number of iteration.
#   seed: seed to use.
#   lr: learning rate.
#   verbose: print output during the estimation.
#   maxpts: maximum number of samples to draw.
#   decay: numeric scalar used to decrease the learning rate.
#   conv_crit: relative convergence threshold.
#   rel_eps: relative error for each term.
#   rescale_vcov: function to rescale the covariance matrix.
#   get_free_vcov_inv: function to get the full covariance matrix which includes
#   fixed parameters.
svrg <- function(par_fn, nobs, val_vcov, val_mea, batch_size, maxit = 10L,
                 seed = 1L, lr, verbose = FALSE, maxpts, decay,
                 conv_crit, rel_eps, rescale_vcov,
                 get_free_vcov_inv){
  n_blocks <- max(1L, ceiling(nobs / batch_size))
  get_blocks <- function(){
    indices <- sample.int(nobs, replace = FALSE) - 1L
    blocks <- tapply(indices, (seq_along(indices) - 1L) %% n_blocks,
                     identity, simplify = FALSE)
  }

  # assign function to get the mean and covariance matrix parameters
  n_mea <- length(val_mea)
  get_par_vcov <- function(par)
    if(n_mea > 0)
      head(par, -n_mea) else par
  get_par_mea <- function(par)
    if(n_mea > 0)
      tail(par,  n_mea) else numeric()

  # assign wrapper for par_fn
  par_fn_wrap <- function(val, comp_derivs, ...){
    out <- par_fn(par_vcov = get_par_vcov(val), par_mea = get_par_mea(val),
                  comp_derivs = comp_derivs, ...)
    if(comp_derivs)
      attr(out, "grad") <- c(attr(out, "grad_vcov"), attr(out, "grad_mea"))
    out
  }

  # assign parameter vector and other needed quantities
  val <- c(val_vcov, val_mea)
  is_vcov <- seq_along(val_vcov)

  n_par <- length(val)
  estimates <- matrix(NA_real_, n_par, maxit + 1L)
  fun_vals <- numeric(maxit + 1L)
  estimates[, 1L] <- val

  lr_use <- lr / decay
  V_mult <- qnorm(1 - .99 / maxit)

  for(k in 1:maxit + 1L){
    blocks <- get_blocks()
    old_val <- estimates[, k - 1L]
    old_grs <- sapply(1:n_blocks - 1L, function(ii){
      idx_b <- (ii %% n_blocks) + 1L
      res_old <- par_fn_wrap(
        old_val, comp_derivs = TRUE, indices = blocks[[idx_b]],
        maxpts[k - 1L])
      c(res_old, attr(res_old, "grad"))
    })

    fun_vals[k - 1L] <- sum(old_grs[1, ])
    old_grs <- old_grs[-1L, , drop = FALSE ]
    old_gr <- rowSums(old_grs) / n_blocks

    lr_use <- lr_use * decay
    for(ii in 1:n_blocks - 1L){
      idx_b <- (ii %% n_blocks) + 1L
      res <- par_fn_wrap(val, comp_derivs = TRUE, indices = blocks[[idx_b]],
                         maxpts[k - 1L])
      fun_vals[k] <- fun_vals[k] + c(res)
      dir <- attr(res, "grad") - old_grs[, ii + 1L] + old_gr

      val <- val + lr_use * dir

      # alter the covariance matrix part
      val[is_vcov] <- .get_lchol(rescale_vcov(.get_lchol_inv(val[is_vcov])))
    }

    estimates[, k] <- val

    if(verbose)
      cat(
        sprintf("End of iteration %4d with learning rate %.8f", k - 1L,
                lr_use),
        sprintf("Log marginal likelihood approximation is %12.2f", fun_vals[k]),
        sprintf("Previous approximate gradient norm was %14.2f\n",
                n_blocks * norm(as.matrix(old_gr))),
        sep = "\n")

    old_v <- fun_vals[k - 1L]
    new_v <- fun_vals[k]
    if(new_v - old_v  < conv_crit * (abs(old_v) + conv_crit))
      break
  }

  list(result = list(
    vcov = get_free_vcov_inv(.get_lchol_inv(get_par_vcov(val))),
    mea = get_par_mea(val)), fun_vals = fun_vals[2:k],
    estimates = estimates[, 2:k, drop = FALSE])
}

# creates a matrix from a log-Cholesky decomposition.
#
# Args:
#   par: p (p + 1) / 2 elements in the log-Cholesky decomposition.
.get_lchol_inv <- function(par){
  # use log-cholesky parametrization
  p <- (sqrt(8 * length(par) + 1) - 1) / 2
  L <- matrix(0, p, p)
  L[lower.tri(L, TRUE)] <- par
  diag(L) <- exp(diag(L))
  tcrossprod(L)
}

# creates the log-Cholesky decomposition.
#
# Args:
#   par: positive definite matrix to decompose
.get_lchol <- function(par){
  lSig <- t(chol(par))
  diag(lSig) <- log(diag(lSig))
  lSig[lower.tri(lSig, TRUE)]
}

#' Impute Missing Values
#'
#' @description
#' Imputes missing values given a covariance matrix and mean vector using a
#' similar quasi-random numbers method as \code{\link{mdgc_log_ml}}.
#'
#' @param object returned object from \code{\link{get_mdgc}}.
#' @param vcov covariance matrix to condition on in the imputation.
#' @param mea vector with non-zero mean entries to condition on.
#' @param abs_eps absolute convergence threshold for each term in the approximation.
#' @param rel_eps relative convergence threshold for each term in the approximation.
#' @param maxit maximum number of samples
#' @inheritParams mdgc_fit
#' @inheritParams mdgc_log_ml
#'
#' @return
#' A list of lists
#' with imputed values for the continuous variables and a vector with
#' probabilities for each level for the ordinal, binary, and multinomial
#' variables.
#'
#' @examples
#' \donttest{
#' # randomly mask data
#' set.seed(11)
#' masked_data <- iris
#' masked_data[matrix(runif(prod(dim(iris))) < .10, NROW(iris))] <- NA
#'
#' # use the functions in the package
#' library(mdgc)
#' obj <- get_mdgc(masked_data)
#' ptr <- get_mdgc_log_ml(obj)
#' start_vals <- mdgc_start_value(obj)
#'
#' fit <- mdgc_fit(ptr, start_vals, obj$means, rel_eps = 1e-2, maxpts = 10000L,
#'                 minvls = 1000L, use_aprx = TRUE, batch_size = 100L, lr = .001,
#'                 maxit = 100L, n_threads = 2L)
#'
#' # impute using the estimated values
#' imputed <- mdgc_impute(obj, fit$result$vcov, fit$result$mea, minvls = 1000L,
#'                        maxit = 10000L, n_threads = 2L, use_aprx = TRUE)
#' imputed[1:5] # first 5 observations
#' head(masked_data, 5) # observed
#' head(iris       , 5) # truth
#' }
#'
#' @export
mdgc_impute <- function(object, vcov, mea, rel_eps = 1e-3, maxit = 10000L,
                        abs_eps = -1, n_threads = 1L, do_reorder = TRUE,
                        minvls = 1000L, use_aprx = FALSE){
  #####
  # checks
  margs <- object$margs
  multinomial <- object$multinomial
  idx_non_zero_mean <- object$idx_non_zero_mean
  nvars <- length(margs)
  if(length(multinomial[[1L]]) > 0)
    nvars <- nvars + sum(multinomial[[1L]][2, ]) - NCOL(multinomial[[1L]])
  stopifnot(
    inherits(object, "mdgc"),
    is.matrix(vcov), is.numeric(vcov), all(dim(vcov) == nvars),
    all(is.finite(vcov)),
    is.numeric(rel_eps), length(rel_eps) == 1L, is.finite(rel_eps),
    is.numeric(abs_eps), length(abs_eps) == 1L, is.finite(abs_eps),
    is.integer(maxit), length(maxit) == 1L, maxit > 0L,
    is.integer(n_threads), length(n_threads) == 1L, n_threads > 0L,
    is.logical(do_reorder), length(do_reorder) == 1L, !is.na(do_reorder),
    is.integer(minvls), length(minvls) == 1L, minvls <= maxit,
    minvls >= 0L,
    is.logical(use_aprx), length(use_aprx) == 1L,
    is.finite(use_aprx),
    length(idx_non_zero_mean) == length(mea),
    all(is.finite(mea)), is.numeric(mea))

  #####
  # perform the imputation
  passed_names <- lapply(margs, function(x){
    x <- attr(x, "levels")
    if(is.null(x))
      x <- ""
    x
  })

  margs_pass <- lapply(margs, function(f){
    if(!is.null(attr(f, "levels")))
      return(f)

    # continuos outcomes; need the inverse emp. cdf
    ev <- environment(f)
    cdf_func <- ev$cdf_func
    eps <- 1 / ev$n
    function(x)
      quantile(cdf_func, max(x, eps))
  })

  mea_arg <- numeric(nvars)
  mea_arg[idx_non_zero_mean + 1L] <- mea

  if(NCOL(multinomial[[1L]]) > 0)
    # TODO: would be nice to avoid this in the C++ code
    for(i in multinomial[[1L]][3, ] + 1L)
      vcov[i, i] <- 1e-20

  impute(lower = object$lower, upper = object$upper, code = object$code,
         Sigma = vcov, mea = mea_arg, truth = object$truth,
         margs = margs_pass, rel_eps = rel_eps, abs_eps = abs_eps,
         maxit = maxit, passed_names = passed_names,
         outer_names = rownames(object$truth), n_threads = n_threads,
         do_reorder = do_reorder, minvls = minvls,
         use_aprx = use_aprx, multinomial = multinomial)
}

#' Perform Model Estimation and Imputation
#'
#' @description
#' A convenience function to perform model estimation and imputation in one
#' call. The learning rate is likely model specific and should be altered.
#' See \code{\link{mdgc_fit}}.
#'
#' See the README at \url{https://github.com/boennecd/mdgc} for examples.
#'
#' @seealso
#' \code{\link{get_mdgc}}, \code{\link{mdgc_start_value}},
#' \code{\link{get_mdgc_log_ml}}, \code{\link{mdgc_fit}},
#' \code{\link{mdgc_impute}}
#'
#' @inheritParams mdgc_fit
#' @inheritParams get_mdgc
#' @param irel_eps relative error for each term in the imputation.
#' @param imaxit maximum number of samples to draw in the imputation.
#' @param iabs_eps absolute convergence threshold for each term in the imputation.
#' @param iminvls minimum number of samples in the imputation.
#' @param start_val starting value for the covariance matrix. Use
#' \code{NULL} if unspecified.
#'
#' @details
#' It is important that the input for \code{data} has the appropriate types and
#' classes. See \code{\link{get_mdgc}}.
#'
#' @references
#' Kingma, D.P., & Ba, J. (2015). \emph{Adam: A Method for Stochastic Optimization}. abs/1412.6980.
#'
#' Johnson, R., & Zhang, T. (2013). \emph{Accelerating stochastic gradient descent using predictive variance reduction}. In Advances in neural information processing systems.
#'
#' @return
#' A list with the following entries:
#'
#' \item{ximp}{\code{\link{data.frame}} with the observed and imputed values.}
#' \item{imputed}{output from \code{\link{mdgc_impute}}.}
#' \item{vcov}{the estimated covariance matrix.}
#' \item{mea}{the estimated non-zero mean terms.}
#'
#' Additional elements may be present depending on the chosen \code{method}.
#' See \code{\link{mdgc_fit}}.
#'
#' @examples
#' \donttest{
#' if(require(catdata)){
#'   data(retinopathy)
#'
#'   # prepare data and save true data set
#'   retinopathy$RET <- as.ordered(retinopathy$RET)
#'   retinopathy$SM <- as.logical(retinopathy$SM)
#'
#'   # randomly mask data
#'   set.seed(28325145)
#'   truth <- retinopathy
#'   for(i in seq_along(retinopathy))
#'     retinopathy[[i]][runif(NROW(retinopathy)) < .3] <- NA
#'
#'   cat("\nMasked data:\n")
#'   print(head(retinopathy, 10))
#'   cat("\n")
#'
#'   # impute data
#'   impu <- mdgc(retinopathy, lr = 1e-3, maxit = 25L, batch_size = 25L,
#'                rel_eps = 1e-3, maxpts = 5000L, verbose = TRUE,
#'                n_threads = 1L, method = "svrg")
#'
#'   # show correlation matrix
#'   cat("\nEstimated correlation matrix\n")
#'   print(impu$vcov)
#'
#'   # compare imputed and true values
#'   cat("\nObserved;\n")
#'   print(head(retinopathy, 10))
#'   cat("\nImputed values:\n")
#'   print(head(impu$ximp, 10))
#'   cat("\nTruth:\n")
#'   print(head(truth, 10))
#'
#'   # using augmented Lagrangian method
#'   cat("\n")
#'   impu_aug <- mdgc(retinopathy, maxit = 25L, rel_eps = 1e-3,
#'                    maxpts = 5000L, verbose = TRUE,
#'                    n_threads = 1L, method = "aug_Lagran")
#'
#'   # compare the log-likelihood estimate
#'   obj <- get_mdgc_log_ml(retinopathy)
#'   cat(sprintf(
#'     "Maximum log likelihood with SVRG vs. augmented Lagrangian:\n  %.2f vs. %.2f\n",
#'     mdgc_log_ml(obj, vcov = impu    $vcov, mea = impu    $mea, rel_eps = 1e-3),
#'     mdgc_log_ml(obj, vcov = impu_aug$vcov, mea = impu_aug$mea, rel_eps = 1e-3)))
#'
#'   # show correlation matrix
#'   cat("\nEstimated correlation matrix (augmented Lagrangian)\n")
#'   print(impu_aug$vcov)
#'
#'   cat("\nImputed values (augmented Lagrangian):\n")
#'   print(head(impu_aug$ximp, 10))
#' }
#' }
#'
#' @export
mdgc <- function(dat, lr = 1e-3, maxit = 25L, batch_size = NULL,
                 rel_eps = 1e-3, method = c("svrg", "adam", "aug_Lagran"),
                 seed = 1L, epsilon = 1e-8,
                 beta_1 = .9, beta_2 = .999, n_threads = 1L,
                 do_reorder = TRUE, abs_eps = -1, maxpts = 10000L,
                 minvls = 100L, verbose = FALSE, irel_eps = rel_eps,
                 imaxit = maxpts, iabs_eps = abs_eps, iminvls = 1000L,
                 start_val = NULL, decay = .98, conv_crit = 1e-5,
                 use_aprx = FALSE){
  mdgc_obj <- get_mdgc(dat)
  p <- NROW(mdgc_obj$lower)
  log_ml_ptr <- get_mdgc_log_ml(mdgc_obj)
  if(is.null(start_val))
    start_val <- mdgc_start_value(mdgc_obj, n_threads = n_threads)
  stopifnot(is.matrix(start_val), all(dim(start_val) == p),
            all(is.finite(start_val)))

  if(verbose)
    cat("Estimating the model...\n")
  fit <- mdgc_fit(
    ptr = log_ml_ptr, vcov = start_val, mea = mdgc_obj$means, lr = lr,
    rel_eps = rel_eps, maxit = maxit, batch_size = batch_size,
    method = method, seed = seed, epsilon = epsilon, beta_1 = beta_1,
    beta_2 = beta_2, n_threads = n_threads, do_reorder = do_reorder,
    abs_eps = abs_eps, maxpts = maxpts, minvls = minvls, verbose = verbose,
    decay = decay, conv_crit = conv_crit, use_aprx = use_aprx)
  vcov <- fit$result$vcov
  mea <- fit$result$mea
  colnames(vcov) <- rownames(vcov) <- rownames(mdgc_obj$lower)

  if(verbose)
    cat("Performing imputation...\n")
  impu <- mdgc_impute(mdgc_obj, vcov,
                      mea = mea,
                      rel_eps = irel_eps,
                      maxit = imaxit, abs_eps = iabs_eps,
                      n_threads = n_threads, do_reorder = do_reorder,
                      minvls = iminvls, use_aprx = use_aprx)
  out <- list(ximp = .threshold(dat, impu), imputed = impu,
              vcov = vcov, mea = mea)
  if(!is.null(fit$fun_vals))
    out$logLik <- fit$fun_vals
  out
}

.threshold <- function(org_data, imputed){
  # checks
  stopifnot(NROW(org_data) == length(imputed),
            is.list(imputed), is.data.frame(org_data))

  # threshold
  is_cont <- which(sapply(org_data, is.numeric))
  is_bin  <- which(sapply(org_data, is.logical))
  is_ord  <- which(sapply(org_data, is.ordered))
  is_mult <- which(sapply(org_data, is.factor))
  is_mult <- setdiff(is_mult, is_ord)
  stopifnot(
    length(is_cont) + length(is_bin) + length(is_ord) + length(is_mult) ==
      NCOL(org_data))
  is_cat <- c(is_bin, is_ord, is_mult)

  trans_to_df <- function(x){
    if(is.matrix(x))
      as.data.frame(t(x))
    else
      as.data.frame(  x )
  }

  out_cont <- if(length(is_cont) == 0)
    NULL else trans_to_df(sapply(imputed, function(x) unlist(x[is_cont])))
  out_cat <- if(length(is_cat) == 0)
    NULL else trans_to_df(sapply(imputed, function(x)
      sapply(x[is_cat], which.max)))
  out <-
    if(is.null(out_cont))
      out_cat else if(is.null(out_cat))
        out_cont else
          cbind(out_cont, out_cat)

  # set factor levels etc.
  out <- out[, order(c(is_cont, is_bin, is_ord, is_mult))]
  if(length(is_bin) > 0)
    out[, is_bin] <- out[, is_bin] > 1L
  if(length(is_ord) > 0)
    for(i in is_ord)
      out[[i]] <- ordered(
        unlist(out[[i]]), labels = levels(org_data[, i]))
  if(length(is_mult) > 0)
    for(i in is_mult)
      out[[i]] <- factor(
        unlist(out[[i]]), labels = levels(org_data[, i]))

  colnames(out) <- colnames(org_data)
  out
}
