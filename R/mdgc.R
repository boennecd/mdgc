#' @importFrom stats na.omit
#' @importFrom utils head
#' @export
get_mgdc <- function(dat){
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
      return(function(x) cbind(NA_real_, scal * cdf_func(pmax(x, mi))))

    } else if(is.logical(x) || is.ordered(x)){
      phat <- table(x) / length(x)
      stopifnot(all(phat > 0))
      lb <- c(0, cumsum(head(phat, -1)))
      ub <- c(cumsum(head(phat, -1)), 1)
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
      return(out)

    }

    stop("type not implemented")
  })

  # get upper and lower bounds or points for the latent variables
  Z <- mapply(function(x, func) func(x), x = dat, func = margs,
                SIMPLIFY = "array")
  lower <- qnorm(t(Z[, 1L, ]))
  upper <- qnorm(t(Z[, 2L, ]))

  # assign codes for missingness pattern
  code <- matrix(0L, NROW(lower), NCOL(lower))
  is_na <- is.na(lower) + is.na(upper)
  code[is_na == 0L] <- 2L
  code[is_na == 2L] <- 1L

  structure(list(
    lower = lower, upper = upper, code = code, margs = margs,
    reals = reals, bins = bins, ords = ords), class = "mdgc")
}

#' @export
get_mgdc_log_ml <- function(object, ...)
  UseMethod("get_mgdc_log_ml")

#' @export
get_mgdc_log_ml.mdgc <- function(object, ...)
  get_mgdc_log_ml.default(lower = object$lower, upper = object$upper,
                          code = object$code)

#' @export
get_mgdc_log_ml.default <- function(object, lower, upper, code, ...){
  # checks
  di <- dim(lower)
  stopifnot(
    di[1] > 0L, di[2] > 0L, length(di) == 2L,
    all(di == dim(upper)), all(di == dim(code)),
    is.numeric(lower), is.numeric(upper), is.integer(code),
    all(is.finite(code)), all(range(code) %in% 0:2))

  out <- get_log_lm_terms_cpp(lower = lower, upper = upper, code = code)
  attr(out, "nobs") <- NCOL(upper)
  attr(out, "nvars") <- NROW(lower)
  out
}

#' Evaluate the Log Marginal Likelihood and Its Derivatives
#' @param ptr object returned by \code{\link{get_mgdc_log_ml}}.
#' @param vcov correlation matrix.
#' @param releps relative error for each term.
#' @param n_threads number of threads.
#' @param comp_derivs logical for whether to approximate the gradient.
#' @param indices integer vector with which terms to include. Must be zero
#' indexed. \code{NULL} yields all observations.
#' @param do_reorder logical for whether to use heuristic variable reordering.
#' @param maxpts maximum number of samples to draw.
#' @param abseps absolute convergence threshold for each term.
#' @export
mgdc_log_ml <- function(ptr, vcov, releps = 1e-2, n_threads = 4L,
                        comp_derivs = FALSE, indices = NULL,
                        do_reorder = TRUE, maxpts = 100000L,
                        abseps = -1.){
  nvars <- attr(ptr, "nvars")
  nobs <- attr(ptr, "nobs")
  stopifnot(!is.null(nvars), !is.null(nobs))

  if(is.null(indices))
    indices <- 0:(nobs - 1L)

  stopifnot(
    all(indices >= 0L & indices < nobs),
    all(dim(vcov) == c(nvars, nvars)),
    is.numeric(releps), length(releps) == 1L, is.finite(releps),
    is.integer(n_threads), length(n_threads) == 1L, n_threads > 0L,
    is.logical(comp_derivs), length(comp_derivs) == 1L, !is.na(comp_derivs),
    is.logical(do_reorder), length(do_reorder) == 1L, !is.na(do_reorder),
    is.integer(maxpts), length(maxpts) == 1L, maxpts > 0L,
    is.numeric(abseps), length(abseps) == 1L, is.finite(abseps))

  eval_log_lm_terms(
    ptr = ptr, vcov = vcov, indices = indices, maxpts = maxpts,
    abseps = abseps, releps = releps, n_threads = n_threads,
    comp_derivs = comp_derivs, do_reorder = do_reorder)
}

#' @export
mgdc_start_value <- function(object, ...)
  UseMethod("mgdc_start_value")

#' @export
mgdc_start_value <- function(object, ...)
  mgdc_start_value.default(lower = object$lower, upper = object$upper,
                           code = object$code)

#' @importFrom stats cov cov2cor
#' @export
mgdc_start_value.default <- function(object, lower, upper, code,
                                     n_threads = 1L, ...){
  Z <- get_z_hat(lower, upper, code, n_threads = n_threads)
  cov2cor(cov(t(Z), use = "pairwise.complete.obs"))
}

#' Estimates the Correlation Matrix
#' @param ptr returned object from \code{\link{get_mgdc_log_ml}}.
#' @param vcov starting value.
#' @param batch_size number of observations in each batch.
#' @param n_threads number of threads to use.
#' @param lr learning rate.
#' @param releps relative error for each term in quasi-Monte Carlo method.
#' @param method estimation method to use.
#' @param maxit maximum number of iteration.
#' @param seed fixed seed to use. Use \code{NULL} if the seed should not be
#' fixed.
#' @param epsilon,beta_1,beta_2 ADAM parameters.
#' @export
mgdc_fit <- function(ptr, vcov, lr = 1e-3, releps = 1e-2,
                     maxit = 10L, batch_size = NULL,
                     method = c("svrg", "adam"), seed = 1L, epsilon = 1e-8,
                     beta_1 = .9, beta_2 = .999, n_threads = 4L){
  #####
  # checks
  nvars <- attr(ptr, "nvars")
  nobs <- attr(ptr, "nobs")
  stopifnot(!is.null(nvars), !is.null(nobs))

  if(is.null(batch_size))
    batch_size <- as.integer(max(min(10, nobs), nobs / 20))

  method <- method[1L]
  stopifnot(
    all(dim(vcov) == c(nvars, nvars)), is.numeric(vcov),
    is.numeric(lr), length(lr) == 1L, lr > 0,
    is.integer(maxit), length(maxit) == 1L, maxit > 0L,
    is.integer(batch_size), length(batch_size) == 1L, batch_size > 0L,
    batch_size <= nobs,
    is.numeric(releps), length(releps) == 1L, releps >= 0.,
    is.character(method), method %in% c("svrg", "adam"),
    is.integer(seed), length(seed) == 1L, is.finite(seed) || is.null(seed),
    is.integer(n_threads), length(n_threads) == 1L, n_threads > 0L,
    is.numeric(epsilon), length(epsilon) == 1L, epsilon >= 0.,
    is.numeric(beta_1), length(beta_1) == 1L, beta_1 >= 0.,
    is.numeric(beta_2), length(beta_2) == 1L, beta_2 >= 0.)

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
  par_fn <- function(par, comp_derivs = FALSE, indices){
    if(!is.null(seed))
      set.seed(seed)
    Arg <- .get_lchol_inv(par)

    res <- mgdc_log_ml(
      ptr = ptr, Arg, comp_derivs = comp_derivs, indices = indices,
      n_threads = n_threads, releps = releps)
    log_ml <- c(res)
    if(comp_derivs){
      gr <- attr(res, "grad")
      tmp <- matrix(0, nvars, nvars)
      tmp[lower.tri(tmp, TRUE)] <- par
      diag(tmp) <- exp(diag(tmp))
      gr <- gr[com_vec] + c(gr)
      gr <- mdgc:::x_dot_X_kron_I(x = gr, X = tmp, l = nvars)
      gr <- gr[, lower.tri(tmp, TRUE)]
      idx_diag <- c(1L, 1L + cumsum(NCOL(tmp):2))
      gr[idx_diag] <- gr[idx_diag] * diag(tmp)

      attr(log_ml, "grad") <- gr

    }

    log_ml
  }

  if(method == "adam")
    return(adam(
      par_fn = par_fn, nobs = nobs, val = .get_lchol(vcov),
      batch_size = batch_size, maxit = maxit, seed = seed,
      epsilon = epsilon, alpha = lr, beta_1 = beta_1, beta_2 = beta_2))
  else if(method == "svrg")
    return(svrg(
      par_fn = par_fn, nobs = nobs, val = .get_lchol(vcov),
      batch_size = batch_size, maxit = maxit, seed = seed, lr = lr))

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
adam <- function(par_fn, nobs, val, batch_size, maxit = 10L,
                 seed = 1L, epsilon = 1e-8, alpha = .001, beta_1 = .9,
                 beta_2 = .999){
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
      res <- par_fn(val, comp_derivs = TRUE, indices = blocks[[idx_b]])
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
svrg <- function(par_fn, nobs, val, batch_size, maxit = 10L, seed = 1L, lr){
  indices <- sample.int(nobs, replace = FALSE) - 1L
  blocks <- tapply(indices, (seq_along(indices) - 1L) %/% batch_size,
                   identity, simplify = FALSE)

  n_blocks <- length(blocks)
  n_par <- length(val)
  estimates <- matrix(NA_real_, n_par, maxit + 1L)
  fun_vals <- numeric(maxit + 1L)
  estimates[, 1L] <- val

  for(k in 1:maxit + 1L){
    old_val <- estimates[, k - 1L]
    old_grs <- sapply(1:n_blocks - 1L, function(ii){
      idx_b <- (ii %% n_blocks) + 1L
      res_old <- par_fn(
        old_val, comp_derivs = TRUE, indices = blocks[[idx_b]])
      c(res_old, attr(res_old, "grad"))
    })

    fun_vals[k - 1L] <- sum(old_grs[1, ])
    old_grs <- old_grs[-1L, , drop = FALSE ]
    old_gr <- rowSums(old_grs) / n_blocks

    for(ii in 1:n_blocks - 1L){
      idx_b <- (ii %% n_blocks) + 1L
      res <- par_fn(val, comp_derivs = TRUE, indices = blocks[[idx_b]])
      fun_vals[k] <- fun_vals[k] + c(res)
      dir <- attr(res, "grad") - old_grs[, ii + 1L] + old_gr

      val <- val + lr * dir
      val <- .get_lchol(cov2cor(.get_lchol_inv(val)))
    }

    estimates[, k] <- val
  }

  list(result = .get_lchol_inv(val), fun_vals = fun_vals[-1L],
       estimates = estimates[, -1L, drop = FALSE])
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

