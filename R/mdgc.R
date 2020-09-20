#' Get mdgc Object
#' @description
#' Creates a mdgc object which is needed for estimation of the
#' correlation matrix and to perform imputation.
#'
#' @param dat \code{\link{data.frame}} with continuous, ordinal, and binary
#' data.
#' @importFrom stats na.omit
#' @importFrom utils head
#' @export
#' @importFrom stats qnorm ecdf quantile
#'
#' @seealso
#' \code{\link{get_mdgc_log_ml}}, \code{\link{mdgc_impute}}
get_mdgc <- function(dat){
  # checks
  stopifnot(is.data.frame(dat), NROW(dat) > 0L)

  # get the types
  reals <- which(sapply(dat, is.numeric))
  bins  <- which(sapply(dat, is.logical))
  ords  <- which(sapply(dat, is.ordered))

  nc <- NCOL(dat)
  stopifnot(length(reals) + length(bins) + length(ords) == nc)

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
      stopifnot(all(phat > 0))
      lb <- c(-Inf, qnorm(cumsum(head(phat, -1))))
      ub <- c(qnorm(cumsum(head(phat, -1))), Inf)

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
      return(out)

    }

    stop("type not implemented")
  })

  # get upper and lower bounds or points for the latent variables
  Z <- mapply(function(x, func) func(x), x = dat, func = margs,
                SIMPLIFY = "array")
  lower <- t(Z[, 1L, ])
  upper <- t(Z[, 2L, ])

  # assign codes for missingness pattern
  code <- matrix(0L, NROW(lower), NCOL(lower))
  is_na <- is.na(lower) + is.na(upper)
  code[is_na == 0L] <- 2L
  code[is_na == 2L] <- 1L

  # get the true values (needed for imputation)
  truth <- t(matrix(as.numeric(unlist(dat)), NROW(dat), NCOL(dat)))

  structure(list(
    lower = lower, upper = upper, code = code, margs = margs,
    reals = reals, bins = bins, ords = ords, truth = truth), class = "mdgc")
}

#' Get Pointer to C++ Object to Approximate the Log Marginal Likelihood
#' @description
#' Creates a C++ object which is needed to approximate the log marginal
#' likelihood. The object cannot be saved.
#'
#' @param object mdgc object from \code{\link{get_mdgc}}. Ignored by the
#' default method.
#' @param lower #variables x #observation matrix with lower bounds
#' for each variable on the normal scale.
#' @param upper #variables x #observation matrix with upper bounds
#' for each variable on the normal scale.
#' @param code #variables x #observation matrix integer code for the
#' each variable on the normal scale. Zero implies an observed value (the
#' value in \code{upper}), one implies a missing value, and two implies an
#' interval.
#' @param ... used to pass arguments to S3 methods.
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
                          code = object$code)

#' @rdname get_mdgc_log_ml
#' @export
get_mdgc_log_ml.default <- function(object, lower, upper, code, ...){
  # checks
  di <- dim(lower)
  stopifnot(
    di[1] > 0L, di[2] > 0L, length(di) == 2L,
    all(di == dim(upper)), all(di == dim(code)),
    is.numeric(lower), is.numeric(upper), is.integer(code),
    all(is.finite(code)), all(range(code) %in% 0:2))

  keep <- colSums(is.na(lower) & is.na(upper)) < NROW(upper)
  out <- get_log_lm_terms_cpp(lower = lower[, keep], upper = upper[, keep],
                              code = code[, keep])
  attr(out, "nobs") <- NCOL(upper)
  attr(out, "nvars") <- NROW(lower)
  out
}

#' Evaluate the Log Marginal Likelihood and Its Derivatives
#'
#' @description
#' Approximates the log marginal likelihood and the derivatives using
#' quasi-random numbers.
#'
#' @seealso
#' \code{\link{mdgc_fit}}
#'
#' @param ptr object returned by \code{\link{get_mdgc_log_ml}}.
#' @param vcov correlation matrix.
#' @param rel_eps relative error for each term.
#' @param n_threads number of threads to use.
#' @param comp_derivs logical for whether to approximate the gradient.
#' @param indices integer vector with which terms (observations) to include.
#' Must be zero indexed. \code{NULL} yields all observations.
#' @param do_reorder logical for whether to use heuristic variable
#' reordering. \code{TRUE} is likely the best option.
#' @param maxpts maximum number of samples to draw.
#' @param abs_eps absolute convergence threshold for each term.
#' @param minvls minimum number of samples.
#' @export
mdgc_log_ml <- function(ptr, vcov, rel_eps = 1e-2, n_threads = 1L,
                        comp_derivs = FALSE, indices = NULL,
                        do_reorder = TRUE, maxpts = 100000L,
                        abs_eps = -1., minvls = 100L){
  nvars <- attr(ptr, "nvars")
  nobs <- attr(ptr, "nobs")
  stopifnot(!is.null(nvars), !is.null(nobs))
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
    minvls >= 0L)

  .mdgc_log_ml(
    ptr = ptr, vcov = vcov, indices = indices, maxpts = maxpts,
    abs_eps = abs_eps, rel_eps = rel_eps, n_threads = n_threads,
    comp_derivs = comp_derivs, do_reorder = do_reorder, minvls = minvls)
}

.mdgc_log_ml <- function(
  ptr, vcov, rel_eps, n_threads, comp_derivs, indices, do_reorder,
  maxpts, abs_eps, minvls)
  eval_log_lm_terms(
    ptr = ptr, vcov = vcov, indices = indices, maxpts = maxpts,
    abs_eps = abs_eps, rel_eps = rel_eps, n_threads = n_threads,
    comp_derivs = comp_derivs, do_reorder = do_reorder, minvls = minvls)

#' Get Starting Value for the Correlation Matrix Using a Heuristic
#'
#' @description
#' Uses a heuristic to get starting values for the correlation matrix. These
#' can be passed e.g. to \code{\link{mdgc_fit}}.
#'
#' @inheritParams get_mdgc_log_ml
#' @param n_threads number of threads to use.
#' @export
mdgc_start_value <- function(...)
  UseMethod("mdgc_start_value")

#' @rdname mdgc_start_value
#' @export
mdgc_start_value <- function(object, ...)
  mdgc_start_value.default(lower = object$lower, upper = object$upper,
                           code = object$code, ...)

#' @rdname mdgc_start_value
#' @importFrom stats cov cov2cor
#' @export
mdgc_start_value.default <- function(object, lower, upper, code,
                                     n_threads = 1L, ...){
  Z <- get_z_hat(lower, upper, code, n_threads = n_threads)
  out <- cov(t(Z), use = "pairwise.complete.obs")

  # handle non-postive definite estimates
  eg <- eigen(out, symmetric = TRUE)
  eps <- 1e-2 * max(eg$values)
  if(any(is_small <- eg$values < eps)){
    eg$values[is_small] <- eps
    out <- tcrossprod(eg$vectors * rep(eg$values, each = NROW(Z)), eg$vectors)
  }

  cov2cor(out)
}

#' Estimate the Correlation Matrix
#'
#' @description
#' Estimates the correlation matrix. The \code{lr} parameter is data
#' and \code{batch_size} dependent. Convergence should be monitored e.g.
#' using \code{verbose = TRUE} with \code{method = "svrg"}.
#'
#' See the README at \url{https://github.com/boennecd/mdgc} for examples.
#'
#' @seealso
#' \code{\link{mdgc_log_ml}}, \code{\link{mdgc_start_value}},
#' \code{\link{mdgc_impute}}.
#'
#' @param ptr returned object from \code{\link{get_mdgc_log_ml}}.
#' @param vcov starting value.
#' @param batch_size number of observations in each batch.
#' @param lr learning rate.
#' @param decay the learning rate used by SVRG is given by \code{lr * decay^iteration_number}.
#' @param method estimation method to use.
#' @param maxit maximum number of iteration.
#' @param seed fixed seed to use. Use \code{NULL} if the seed should not be
#' fixed.
#' @param epsilon,beta_1,beta_2 ADAM parameters.
#' @param conv_crit relative convergence threshold.
#' @param verbose logical for whether to print output during the estimation.
#' @inheritParams mdgc_log_ml
#'
#' @references
#' Kingma, D.P., & Ba, J. (2015). \emph{Adam: A Method for Stochastic Optimization}. abs/1412.6980.
#'
#' Johnson, R., & Zhang, T. (2013). \emph{Accelerating stochastic gradient descent using predictive variance reduction}. In Advances in neural information processing systems.
#'
#' @export
mdgc_fit <- function(ptr, vcov, lr = 1e-3, rel_eps = 1e-3,
                     maxit = 10L, batch_size = NULL,
                     method = c("svrg", "adam"), seed = 1L, epsilon = 1e-8,
                     beta_1 = .9, beta_2 = .999, n_threads = 1L,
                     do_reorder = TRUE, abs_eps = -1., maxpts = 10000L,
                     minvls = 100L, verbose = FALSE, decay = .98,
                     conv_crit = 1e-6){
  #####
  # checks
  nvars <- attr(ptr, "nvars")
  nobs <- attr(ptr, "nobs")
  stopifnot(!is.null(nvars), !is.null(nobs))

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
    is.character(method), method %in% c("svrg", "adam"),
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
    is.numeric(conv_crit), length(conv_crit) == 1L, is.finite(conv_crit))

  #####
  # assign functions to use
  # indices used to apply a matrix product with a get_commutation matrix
  com_vec <- get_commutation_vec(nvars, nvars, FALSE)

  # computes the approximate log marginal likelihood.
  #
  # Args:
  #   par: log-Cholesky decomposition.
  #   comp_derivs: logical for whether to approximate the gradient.
  #   indices: indices to use.
  #   maxpts: maximum number of samples to draw.
  par_fn <- function(par, comp_derivs = FALSE, indices,
                     maxpts){
    if(!is.null(seed))
      set.seed(seed)
    Arg <- .get_lchol_inv(par)

    res <- .mdgc_log_ml(
      ptr = ptr, Arg, comp_derivs = comp_derivs, indices = indices,
      n_threads = n_threads, rel_eps = rel_eps, do_reorder = do_reorder,
      abs_eps = abs_eps, maxpts = maxpts, minvls = minvls)
    log_ml <- c(res)
    if(comp_derivs){
      gr <- attr(res, "grad")
      tmp <- matrix(0, nvars, nvars)
      tmp[lower.tri(tmp, TRUE)] <- par
      diag(tmp) <- exp(diag(tmp))
      gr <- gr[com_vec] + c(gr)
      gr <- x_dot_X_kron_I(x = gr, X = tmp, l = nvars)
      gr <- gr[, lower.tri(tmp, TRUE)]
      idx_diag <- c(1L, 1L + cumsum(NCOL(tmp):2))
      gr[idx_diag] <- gr[idx_diag] * diag(tmp)

      attr(log_ml, "grad") <- gr

    }

    log_ml
  }

  maxpts_base <- exp(log(minvls / maxpts) / (maxit - 1))
  maxpts_scale <- maxpts_base^(maxit:1)
  maxpts_scale <- maxpts_scale / max(maxpts_scale)
  maxpts_use <- pmax(minvls, as.integer(floor(maxpts * maxpts_scale)))

  if(method == "adam")
    return(adam(
      par_fn = par_fn, nobs = nobs, val = .get_lchol(vcov),
      batch_size = batch_size, maxit = maxit, seed = seed,
      epsilon = epsilon, alpha = lr, beta_1 = beta_1, beta_2 = beta_2,
      maxpts = maxpts_use))
  else if(method == "svrg")
    return(svrg(
      par_fn = par_fn, nobs = nobs, val = .get_lchol(vcov),
      batch_size = batch_size, maxit = maxit, seed = seed, lr = lr,
      verbose = verbose, maxpts = maxpts_use, decay = decay,
      conv_crit = conv_crit, rel_eps = rel_eps))

  stop(sprintf("Method '%s' is not implemented", method))
}

#####
# performs stochastic gradient descent instead (using ADAM).
#
# Args:
#   par_fn: function to evaluate the log marginal likelihood.
#   nobs: number of observation.
#   val: starting value.
#   batch_size: number of observations in each batch.
#   maxit: maximum number of iteration.
#   seed: seed to use.
#   epsilon, alpha, beta_1, beta_2: ADAM parameters.
#   maxpts: maximum number of samples to draw.
adam <- function(par_fn, nobs, val, batch_size, maxit = 10L,
                 seed = 1L, epsilon = 1e-8, alpha = .001, beta_1 = .9,
                 beta_2 = .999, maxpts){
  indices <- sample.int(nobs, replace = FALSE) - 1L
  blocks <- tapply(indices, (seq_along(indices) - 1L) %/% batch_size,
                   identity, simplify = FALSE)

  n_blocks <- length(blocks)
  n_par <- length(val)
  m <- v <- numeric(n_par)
  fun_vals <- numeric(maxit)
  estimates <- matrix(NA_real_, n_par, maxit)
  i <- -1L

  for(k in 1:maxit){
    for(ii in 1:n_blocks){
      i <- i + 1L
      idx_b <- (i %% n_blocks) + 1L
      m_old <- m
      v_old <- v
      res <- par_fn(val, comp_derivs = TRUE, indices = blocks[[idx_b]],
                    maxpts[k])
      fun_vals[(i %/% n_blocks) + 1L] <-
        fun_vals[(i %/% n_blocks) + 1L] + c(res)

      gr <- attr(res, "grad")

      m <- beta_1 * m_old + (1 - beta_1) * gr
      v <- beta_2 * v_old + (1 - beta_2) * gr^2

      m_hat <- m / (1 - beta_1^(i + 1))
      v_hat <- v / (1 - beta_2^(i + 1))

      val <- val + alpha * m_hat / (sqrt(v_hat) + epsilon)
      val <- .get_lchol(cov2cor(.get_lchol_inv(val)))
    }

    estimates[, k] <- val
  }

  list(result = .get_lchol_inv(val), estimates = estimates)
}

#####
# performs stochastic gradient descent instead (using SVRG).
#
# Args:
#   par_fn: function to evaluate the log marginal likelihood.
#   nobs: number of observation.
#   val: starting value.
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
svrg <- function(par_fn, nobs, val, batch_size, maxit = 10L, seed = 1L, lr,
                 verbose = FALSE, maxpts, decay, conv_crit, rel_eps){
  indices <- sample.int(nobs, replace = FALSE) - 1L
  blocks <- tapply(indices, (seq_along(indices) - 1L) %/% batch_size,
                   identity, simplify = FALSE)

  n_blocks <- length(blocks)
  n_par <- length(val)
  estimates <- matrix(NA_real_, n_par, maxit + 1L)
  fun_vals <- numeric(maxit + 1L)
  estimates[, 1L] <- val

  lr_use <- lr / decay
  V_mult <- qnorm(1 - .99 / maxit)
  for(k in 1:maxit + 1L){
    old_val <- estimates[, k - 1L]
    old_grs <- sapply(1:n_blocks - 1L, function(ii){
      idx_b <- (ii %% n_blocks) + 1L
      res_old <- par_fn(
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
      res <- par_fn(val, comp_derivs = TRUE, indices = blocks[[idx_b]],
                    maxpts[k - 1L])
      fun_vals[k] <- fun_vals[k] + c(res)
      dir <- attr(res, "grad") - old_grs[, ii + 1L] + old_gr

      val <- val + lr_use * dir
      val <- .get_lchol(cov2cor(.get_lchol_inv(val)))
    }

    estimates[, k] <- val

    if(verbose)
      cat(
        sprintf("End if iteration %4d with learning rate %.8f", k - 1L,
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

  list(result = .get_lchol_inv(val), fun_vals = fun_vals[2:k],
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

#' Impute Missing Values Given Correlation Matrix
#'
#' @description
#' Imputes missing values given a correlation matrix using a similar
#' quasi-random numbers method as \code{\link{mdgc_fit}}.
#'
#' @param object returned object from \code{\link{get_mdgc}}.
#' @param vcov correlation matrix to condition on in the imputation.
#' @inheritParams mdgc_fit
#' @inheritParams mdgc_log_ml
#' @export
#'
#' @return
#' A list with imputed values for the continuous and a vector with
#' probabilities for each categorical.
mdgc_impute <- function(object, vcov, rel_eps = 1e-3, maxit = 10000L,
                        abs_eps = -1, n_threads = 1L, do_reorder = TRUE,
                        minvls = 1000L){
  #####
  # checks
  margs <- object$margs
  nvars <- length(margs)
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
    minvls >= 0L)

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

  impute(lower = object$lower, upper = object$upper, code = object$code,
         Sigma = vcov, truth = object$truth,
         margs = margs_pass, rel_eps = rel_eps, abs_eps = abs_eps,
         maxit = maxit, passed_names = passed_names,
         outer_names = rownames(object$lower), n_threads = n_threads,
         do_reorder = do_reorder, minvls = minvls)
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
#' @param start_val starting value for the correlation matrix. Use \code{NULL} if unspecified.
#' @export
#'
#' @examples
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
#'                n_threads = 1L)
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
#' }
mdgc <- function(dat, lr = 1e-3, maxit = 10L, batch_size = NULL,
                 rel_eps = 1e-3,
                 method = c("svrg", "adam"), seed = 1L, epsilon = 1e-8,
                 beta_1 = .9, beta_2 = .999, n_threads = 1L,
                 do_reorder = TRUE, abs_eps = -1, maxpts = 10000L,
                 minvls = 100L, verbose = FALSE, irel_eps = rel_eps,
                 imaxit = maxpts, iabs_eps = abs_eps, iminvls = 1000L,
                 start_val = NULL, decay = .98, conv_crit = 1e-5){
  mdgc_obj <- get_mdgc(dat)
  p <- NCOL(dat)
  log_ml_ptr <- get_mdgc_log_ml(mdgc_obj)
  if(is.null(start_val))
    start_val <- mdgc_start_value(mdgc_obj, n_threads = n_threads)
  stopifnot(is.matrix(start_val), all(dim(start_val) == p),
            all(is.finite(start_val)))

  if(verbose)
    cat("Estimating the model...\n")
  fit <- mdgc_fit(
    ptr = log_ml_ptr, vcov = start_val, lr = lr, rel_eps = rel_eps,
    maxit = maxit, batch_size = batch_size, method = method, seed = seed,
    epsilon = epsilon, beta_1 = beta_1, beta_2 = beta_2,
    n_threads = n_threads, do_reorder = do_reorder, abs_eps = abs_eps,
    maxpts = maxpts, minvls = minvls, verbose = verbose, decay = decay,
    conv_crit = conv_crit)
  vcov <- fit$result
  colnames(vcov) <- rownames(vcov) <- colnames(dat)

  if(verbose)
    cat("Performing imputation...\n")
  impu <- mdgc_impute(mdgc_obj, fit$result, rel_eps = irel_eps,
                      maxit = imaxit, abs_eps = iabs_eps,
                      n_threads = n_threads, do_reorder = do_reorder,
                      minvls = iminvls)
  out <- list(ximp = .threshold(dat, impu), imputed = impu,
              vcov = vcov)
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
  stopifnot(
    length(is_cont) + length(is_bin) + length(is_ord) == NCOL(org_data))
  is_cat <- c(is_bin, is_ord)

  trans_to_df <- function(x){
    if(is.matrix(x))
      as.data.frame(t(x))
    else
      as.data.frame(  x )
  }

  out_cont <- trans_to_df(sapply(imputed, function(x) unlist(x[is_cont])))
  out_cat <- trans_to_df(sapply(imputed, function(x)
    sapply(x[is_cat], which.max)))
  out <- cbind(out_cont, out_cat)

  # set factor levels etc.
  out <- out[, order(c(is_cont, is_bin, is_ord))]
  if(length(is_bin) > 0)
    out[, is_bin] <- out[, is_bin] > 1L
  if(length(is_ord) > 0)
    for(i in is_ord)
      out[[i]] <- ordered(
        unlist(out[[i]]), labels = levels(org_data[, i]))

  colnames(out) <- colnames(org_data)
  out
}
