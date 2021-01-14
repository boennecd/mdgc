context("Testing 'get_log_lm_terms'")

# sim_dat <- function(n, p = 4){
#   # get the covariance matrix
#   Sb <- diag(p)
#   Sb[lower.tri(Sb)] <- Sb[upper.tri(Sb)] <- .5
#   Sb <- Sb / p / 5
#   Sig <- cov2cor(drop(rWishart(1L, 5L * p, Sb)))
#
#   # draw the observations
#   Z <- truth <- crossprod(chol(Sig), matrix(rnorm(n * p), p))
#
#   # mask
#   is_mask <- matrix(runif(n * p) < .45, p)
#   is_int <- ceiling(p / 2):p
#   is_mask[is_int, ] <- is_mask[is_int, ] & Z[is_int, ] < 0
#
#   Z[ is_int, ][is_mask[ is_int, ]] <- 0
#   Z[-is_int, ][is_mask[-is_int, ]] <- NA_real_
#
#   # create matrix in the Z format to pass to c++
#   lower <- matrix(-Inf, p, n)
#   upper <- Z
#   # codes are:
#   #  0: latent Z is observed (upper is the observed point).
#   #  1: latent Z can be anything..
#   #  2: latent Z is in an interval.
#   code <- matrix(0L, p, n)
#   code[-is_int, ][is_mask[-is_int, ]] <- 1L
#   code[ is_int, ][is_mask[ is_int, ]] <- 2L
#
#   list(lower = lower, upper = upper, code = code, Sigma = Sig,
#        truth = truth)
# }
#
# set.seed(3)
# p <- 5L
# dat <- sim_dat(10L, p = p)
# saveRDS(dat, "get_log_lm_terms-test.RDS")
dat <- readRDS("get_log_lm_terms-test.RDS")

test_that("'get_log_lm_terms' gives the correct result with and without gradients", {
  multinomial <- replicate(NCOL(dat$lower),
                           matrix(0L, 0, 0), simplify = FALSE)
  idx_non_zero_mean <- 1:4
  ptr <- mdgc:::get_log_lm_terms_cpp(
    lower = dat$lower, upper = dat$upper, code = dat$code,
    multinomial = multinomial, idx_non_zero_mean = idx_non_zero_mean)

  lcov_to_mat <- function(par){
    p <- (sqrt(8 * length(par) + 1) - 1) / 2
    L <- matrix(0, p, p)
    L[lower.tri(L, TRUE)] <- par
    L[upper.tri(L)] <- t(L)[upper.tri(L)]
    L
  }

  log_ml <- function(par, rel_eps = 1e-5, n_threads = 1L,
                     comp_derivs = FALSE, seed = NULL){
    if(!is.null(seed))
      set.seed(seed)
    vcov_log_chol <- head(par, -length(idx_non_zero_mean))
    mea           <- tail(par,  length(idx_non_zero_mean))
    Arg <- lcov_to_mat(vcov_log_chol)

    mdgc:::eval_log_lm_terms(
      ptr = ptr, vcov = Arg, mu = mea,
      indices = 0:(NCOL(dat$lower) - 1L),
      maxpts = 1000000L, abs_eps = -1, rel_eps = rel_eps, n_threads = n_threads,
      comp_derivs = comp_derivs, minvls = 0L)
  }

  set.seed(1)
  pa <- c(dat$Sigma[lower.tri(dat$Sigma, TRUE)],
          numeric(length(idx_non_zero_mean)))
  val <- log_ml(pa)
  # dput(val)
  expect_equal(val, -45.8852616647973, tolerance = 1e-3)

  # with more threads
  val <- log_ml(pa, n_threads = 2)
  expect_equal(val, -45.8852616647973, tolerance = 1e-3)

  grad <- log_ml(pa, comp_derivs = TRUE)
  # truth <- numDeriv::jacobian(log_ml, pa, seed = 1L)
  # dput(truth)
  true_grad <-c(
    4.44392200186265, 8.26645723733806, 0.0770691610564516,
    1.28391853303438, -10.3954140182037, -0.48739982069906, 2.03203577758658,
    -0.105887310232136, -5.21844960723815, -4.20995247942651, -2.94589717845707,
    4.71792106473805, -1.62254861814731, 2.75032259524273, 0.70355964710631,
    0, -2.2618892073704, -2.70197542613998, -2.09928723929281)

  # dput(matrixcalc::duplication.matrix(5L))
  jac <-structure(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    1), .Dim = c(25L, 15L))
  expect_equal(
    drop(c(attr(grad, "grad_vcov")) %*% jac), head(true_grad, -4),
    tolerance = 1e-3)
  expect_equal(
    drop(attr(grad, "grad_mea")), tail(true_grad, 4),
    tolerance = 1e-3)

  # now with two threads
  grad <- log_ml(pa, comp_derivs = TRUE, n_threads = 2L)
  expect_equal(
    drop(c(attr(grad, "grad_vcov")) %*% jac), head(true_grad, -4),
    tolerance = 1e-3)
  expect_equal(
    drop(attr(grad, "grad_mea")), tail(true_grad, 4),
    tolerance = 1e-3)

  # with a different mean
  pa <- c(dat$Sigma[lower.tri(dat$Sigma, TRUE)],
          c(-1, -1, .5, 1))

  grad <- log_ml(pa, comp_derivs = TRUE)
  # truth <- numDeriv::jacobian(log_ml, pa, seed = 1L)
  # dput(truth)
  true_grad <- c(
    4.4536865067216, 8.41361831323731, 0.0677047074608158,
    1.34830327174643, -10.4850097822861, 0.715754440295721, 2.01389217846365,
    -2.26194994919321, -8.33693028534194, -5.47086432114414, -1.71520257450454,
    5.56054925709124, 0.305752972603772, 1.77182637497433, 4.35720858728649,
    0, -0.667362476405746, -4.1050781283087, -6.14566749995671)
  expect_equal(
    drop(c(attr(grad, "grad_vcov")) %*% jac), head(true_grad, -4),
    tolerance = 1e-3)
  expect_equal(
    drop(attr(grad, "grad_mea")), tail(true_grad, 4),
    tolerance = 1e-3)
})

test_that("'get_log_lm_terms' gives the correct result with and without gradients with multinomial variables", {
  skip_if_not_installed("datasets")
  library(datasets)
  dat <- iris[c(1:10, 51:60, 101:110), 3:5]
  dat[seq(1, NROW(dat), by = 4), 1] <- NA
  dat[seq(2, NROW(dat), by = 4), 2] <- NA
  dat[seq(3, NROW(dat), by = 4), 3] <- NA

  obj <- get_mdgc(dat)
  ptr <- get_mdgc_log_ml(obj)
  start <- mdgc_start_value(obj)

  lcov_to_mat <- function(par){
    p <- (sqrt(8 * length(par) + 1) - 1) / 2
    L <- matrix(0, p, p)
    L[lower.tri(L, TRUE)] <- par
    L[upper.tri(L)] <- t(L)[upper.tri(L)]
    L
  }

  log_ml <- function(par, rel_eps = 1e-5, n_threads = 1L,
                     comp_derivs = FALSE, seed = NULL){
    if(!is.null(seed))
      set.seed(seed)
    vcov_log_chol <- head(par, -2)
    mea           <- tail(par,  2)
    Arg <- lcov_to_mat(vcov_log_chol)

    mdgc_log_ml(ptr = ptr, vcov = Arg, mea = mea,
                maxpts = 1000000L, abs_eps = -1, rel_eps = rel_eps,
                n_threads = n_threads, comp_derivs = comp_derivs,
                minvls = 0L)
  }

  set.seed(1)
  pa <- c(start[lower.tri(start, TRUE)], obj$means)
  val <- log_ml(pa)
  # dput(val)
  expect_equal(val, -62.496433767902, tolerance = 1e-3)

  grad <- log_ml(pa, comp_derivs = TRUE)
  # truth <- numDeriv::jacobian(log_ml, pa, seed = 1L)
  # dput(truth)
  true_grad <-c(
    -10.8490252758386, 11.3412044896911, -8.09698535456282,
    0.313894770695329, 7.78309064063369, -14.9804883174154, -10.4020331443358,
    1.3207808003981, 9.08125234387075, -7.69933557281718, 1.2228223849216,
    14.1758487591125, -0.538639488002858, -0.145543412030162, -7.01515267411361,
    -0.262549150275879, -1.05471744166441)

  # dput(matrixcalc::duplication.matrix(5L))
  jac <-structure(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    1), .Dim = c(25L, 15L))
  expect_equal(
    drop(c(attr(grad, "grad_vcov")) %*% jac), head(true_grad, -2),
    tolerance = 1e-3)
  expect_equal(
    drop(attr(grad, "grad_mea")), tail(true_grad, 2),
    tolerance = 1e-3)

  # now with two threads
  grad <- log_ml(pa, comp_derivs = TRUE, n_threads = 2L)
  expect_equal(
    drop(c(attr(grad, "grad_vcov")) %*% jac), head(true_grad, -2),
    tolerance = 1e-3)
  expect_equal(
    drop(attr(grad, "grad_mea")), tail(true_grad, 2),
    tolerance = 1e-3)
})

# sim_dat <- function(n, p = 4, n_lvls = 5L){
#   # get the covariance matrix
#   Sb <- diag(p)
#   Sb[lower.tri(Sb)] <- Sb[upper.tri(Sb)] <- .5
#   Sb <- Sb / p / 5
#   Sig <- cov2cor(drop(rWishart(1L, 5L * p, Sb)))
#
#   # draw the observations
#   Z <- truth <- crossprod(chol(Sig), matrix(rnorm(n * p), p))
#
#   # determine the type
#   type <- rep(1:3, each = floor((p + 3 - 1) / 3))[1:p]
#   is_con <- type == 1L
#   is_bin <- type == 2L
#   is_ord <- type == 3L
#
#   # mask data
#   is_mask <- matrix(runif(n * p) < .33, p)
#   Z[is_mask] <- NA_real_
#
#   # set code arugment. Codes are:
#   #  0: latent Z is observed (upper is the observed point).
#   #  1: latent Z can be anything..
#   #  2: latent Z is in an interval.
#   code <- matrix(0L, p, n)
#   code[is_ord | is_bin, ] <- 2L
#   code[is_mask] <- 1L
#
#   # create upper and lower bounds
#   upper <- lower <- matrix(NA_real_, p, n)
#   lower[is_con, ] <- upper[is_con, ] <- truth[is_con, ]
#
#   bs <- c(-Inf, 0, Inf)
#   lower[is_bin, ] <- head(bs, -1)[cut(truth[is_bin, ], breaks = bs)]
#   upper[is_bin, ] <- bs[-1]      [cut(truth[is_bin, ], breaks = bs)]
#
#   bs <- qnorm(seq(0, 1, length.out = n_lvls + 1L))
#   lower[is_ord, ] <- head(bs, -1)[cut(truth[is_ord, ], breaks = bs)]
#   upper[is_ord, ] <- bs[-1]      [cut(truth[is_ord, ], breaks = bs)]
#
#   lower[is_mask] <- upper[is_mask] <- NA_real_
#
#   list(lower = lower, upper = upper, code = code, Sigma = Sig,
#        truth = truth)
# }
# set.seed(3)
# p <- 6L
# dat <- sim_dat(10L, p = p)
# saveRDS(dat, "get_log_lm_terms-test-ord-bin.RDS")
dat <- readRDS("get_log_lm_terms-test-ord-bin.RDS")

test_that("'get_log_lm_terms' gives the correct result with and without gradients (ordinal and binary data)", {
  multinomial <- replicate(NCOL(dat$lower),
                           matrix(0L, 0, 0), simplify = FALSE)
  idx_non_zero_mean <- 1:5
  ptr <- mdgc:::get_log_lm_terms_cpp(
    lower = dat$lower, upper = dat$upper, code = dat$code,
    multinomial = multinomial, idx_non_zero_mean = idx_non_zero_mean)

  lcov_to_mat <- function(par){
    p <- (sqrt(8 * length(par) + 1) - 1) / 2
    L <- matrix(0, p, p)
    L[lower.tri(L, TRUE)] <- par
    L[upper.tri(L)] <- t(L)[upper.tri(L)]
    L
  }

  n_mea <- length(idx_non_zero_mean)
  log_ml <- function(par, rel_eps = 1e-5, n_threads = 1L,
                     comp_derivs = FALSE, seed = NULL){
    if(!is.null(seed))
      set.seed(seed)
    Arg <- lcov_to_mat(head(par, -n_mea))
    mea <- tail(par, n_mea)

    mdgc:::eval_log_lm_terms(
      ptr = ptr, vcov = Arg, mu = mea, indices = 0:(NCOL(dat$lower) - 1L),
      maxpts = 1000000L, abs_eps = -1, rel_eps = rel_eps, n_threads = n_threads,
      comp_derivs = comp_derivs, minvls = 0L)
  }

  set.seed(1)
  par <- c(dat$Sigma[lower.tri(dat$Sigma, TRUE)], numeric(n_mea))
  val <- log_ml(par)
  # dput(val)
  expect_equal(val, -50.633131628025, tolerance = 1e-3)

  val <- log_ml(par, n_threads = 2)
  expect_equal(val, -50.633131628025, tolerance = 1e-3)

  grad <- log_ml(par, comp_derivs = TRUE)
  # truth <- numDeriv::jacobian(log_ml, par, seed = 1L)
  # dput(truth)
  true_grad <- c(
    2.18828846268184, 1.14819853573079, -0.833095815058999,
    -1.64775972436418, -3.02462052022083, 2.38038644284807, 4.02535849827106,
    -0.6823072284458, -1.52397119919305, -5.83385872168417, -0.526004633937689,
    -0.772721040056757, -2.15498992108154, 4.76760990982455, -0.156298143778767,
    0.415973653298547, -0.190202356593625, 2.98216021926465, 1.03082094344473,
    0.816776256035582, -1.55995256546344,
    0, -0.151148246224951, 3.97666677411463,
    -3.30058060270238, 2.99870854999467)

  # dput(matrixcalc::duplication.matrix(5L))
  jac <-structure(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 1), .Dim = c(36L, 21L))
  expect_equal(drop(c(attr(grad, "grad_vcov")) %*% jac),
               head(true_grad, -n_mea), tolerance = 1e-3)
  expect_equal(drop(attr(grad, "grad_mea")), tail(true_grad, n_mea),
               tolerance = 1e-3)

  grad <- log_ml(par, comp_derivs = TRUE, n_threads = 2L)
  expect_equal(drop(c(attr(grad, "grad_vcov")) %*% jac),
               head(true_grad, -n_mea), tolerance = 1e-3)
  expect_equal(drop(attr(grad, "grad_mea")), tail(true_grad, n_mea),
               tolerance = 1e-3)
})

test_that("mdgc_log_ml gives the correct result with a single log marginal likelihood term", {
  n_lvl <- c(3L, 4L, 5L)
  n_latent <- n_lvl[1] + n_lvl[3] + 4L
  # dput(cov2cor(drop(rWishart(1L, 2 * n_latent, diag(n_latent)))))
  Sig <- structure(c(1, -0.139976486533139, -0.252394041490724, 0.0328479391884292, 0.33209161389653, -0.0802581313618109, 0.142279177667647, 0.114400486140885, 0.0084305016074294, -0.208392048350881, 0.0241710440085548, 0.0692588360372489, -0.139976486533139, 1, 0.15006370015725, 0.0387966113834449, 0.040653657443624, -0.163867176796534, -0.163776074146943, -0.00677335530339495, -0.0962742888220896, -0.363172104047447, -0.0345447459494092, 0.017633351661753, -0.252394041490724, 0.15006370015725, 1, -0.2576125675092, -0.120204683054427, 0.325392375081554, 0.0617473677009837, -0.323965717346053, 0.172023461202382, -0.0483353723590527, 0.00743318236676117, 0.188726831046407, 0.0328479391884292, 0.0387966113834449, -0.2576125675092, 1, 0.369589014205649, -0.226449781307351, -0.232296759237575, 0.119732751860623, -0.341216668163613, -0.117472102794596, -0.208383502686246, -0.177819043633142, 0.33209161389653, 0.040653657443624, -0.120204683054427, 0.369589014205649, 1, -0.0319803022599951, -0.329435492022031, 0.0999548164047468, 0.0113250918738276, 0.0845154592053214, -0.12843295462476, -0.0529800713738894, -0.0802581313618109, -0.163867176796534, 0.325392375081554, -0.226449781307351, -0.0319803022599951, 1, 0.0540971106518889, -0.263844638490977, 0.10682435039763, 0.0102738072433955, -0.115033109688909, 0.0885163680527274, 0.142279177667647, -0.163776074146943, 0.0617473677009837, -0.232296759237575, -0.329435492022031, 0.0540971106518889, 1, 0.282127993763662, 0.196370386133481, -0.262622833662247, 0.0673683870070675, -0.0443091951152617, 0.114400486140885, -0.00677335530339495, -0.323965717346053, 0.119732751860623, 0.0999548164047468, -0.263844638490977, 0.282127993763662, 1, -0.156072363656726, 0.0454778524087398, 0.205042507925608, -0.0150339689622075, 0.0084305016074294, -0.0962742888220896, 0.172023461202382, -0.341216668163613, 0.0113250918738276, 0.10682435039763, 0.196370386133481, -0.156072363656726, 1, 0.156828343541278, -0.0063302709111426, 0.0523227195166851, -0.208392048350881, -0.363172104047447, -0.0483353723590527, -0.117472102794596, 0.0845154592053214, 0.0102738072433955, -0.262622833662247, 0.0454778524087398, 0.156828343541278, 1, 0.191638870428741, -0.0762526882470861, 0.0241710440085548, -0.0345447459494092, 0.00743318236676117, -0.208383502686246, -0.12843295462476, -0.115033109688909, 0.0673683870070675, 0.205042507925608, -0.0063302709111426, 0.191638870428741, 1, 0.184869417557784, 0.0692588360372489, 0.017633351661753, 0.188726831046407, -0.177819043633142, -0.0529800713738894, 0.0885163680527274, -0.0443091951152617, -0.0150339689622075, 0.0523227195166851, -0.0762526882470861, 0.184869417557784, 1), .Dim = c(12L, 12L ))
  # dput(c(rnorm(n_lvl[1]), 0, 0, rnorm(n_latent - n_lvl[1] - 2L)))
  mu <- c(1.37414772401687, 0.287678504785475, -0.262694211676933, 0, 0, 1.89638937266271, 0.0806366010191161, 0.853503953969782, -1.30949921135887, -0.639517569806596, 0.724772050061993, -0.46190976626267)
  # dput(mvtnorm::rmvnorm(1L, mu, Sig))
  obs_x <- c(0.698464938402946, -0.91241587967978, -0.588384537744868, -0.566074730636501, 1.22396370835824, 1.04393663847435, 0.134511857961685, 1.7761587334823, -0.731218046435423, 0.28578514459854, -0.083723324156426, -1.25132416246612)

  O_lvls <- qnorm(seq(0, 1, length.out = n_lvl[2L] + 1L))
  obs <- data.frame(
    M1 = factor(which.min(obs_x[1:n_lvl[1L]]), levels = 1:n_lvl[1L]),
    C1 = obs_x[n_lvl[1L] + 1L],
    C2 = obs_x[n_lvl[1L] + 2L],
    B1 = obs_x[n_lvl[1L] + 3L] > 0,
    O1 = cut(obs_x[n_lvl[1L] + 4L], breaks = O_lvls),
    M2 = factor(which.min(tail(obs_x, n_lvl[3L])), levels = 1:n_lvl[3L]))
  obs$O1 <- ordered(obs$O1, labels = levels(obs$O1),
                    levels = levels(obs$O1))

  # setup
  lower <- rep(-Inf, n_latent) # default to -Inf
  upper <- numeric(n_latent)   # default to zero

  # handle continuous
  upper[n_lvl[1L] + 1:2] <- obs_x[n_lvl[1L] + 1:2]

  # handle the binary
  if(obs$B1){
    lower[n_lvl[1L] + 3L] <- -mu[n_lvl[1L] + 3L]
    upper[n_lvl[1L] + 3L] <- Inf
  } else
    upper[n_lvl[1L] + 3L] <- -mu[n_lvl[1L] + 3L]

  # handle the ordinal
  lower[n_lvl[1L] + 4L] <-
    O_lvls[as.integer(obs$O1)     ] - mu[n_lvl[1L] + 4L]
  upper[n_lvl[1L] + 4L] <-
    O_lvls[as.integer(obs$O1) + 1L] - mu[n_lvl[1L] + 4L]

  # the multinomial
  upper[1:n_lvl[1L]] <- -mu[1:n_lvl[1L]]
  upper[n_latent - n_lvl[3L]:1 + 1L] <- -tail(mu, n_lvl[3L])
  code <- c(rep(2L, n_lvl[1L]), 0L, 0L, rep(2L, n_latent - n_lvl[1L] - 2L))

  # setup the multinomial argument
  multinomial <- list(rbind(
    c(as.integer(obs$M1), as.integer(obs$M2)), n_lvl[c(1, 3)],
    c(0L, n_latent - n_lvl[3])))

  # compute value to compare with
  is_con <- n_lvl[1L] + 1:2
  is_int <- setdiff(1:n_latent, is_con)
  # dput(mvtnorm::dmvnorm(obs_x[is_con], mu[is_con], Sig[is_con, is_con], log = TRUE))
  f1 <- -3.11413974714096

  Sig_int <- Sig[is_int, is_int] - Sig[is_int, is_con] %*% solve(
    Sig[is_con, is_con], Sig[is_con, is_int])
  mu_int <- drop(mu[is_int] + Sig[is_int, is_con] %*%
    solve(Sig[is_con, is_con], obs_x[is_con] - mu[is_con]))
  n_int <- length(mu_int)

  D <- matrix(0., n_int - 2L, 2L)
  D[1:(n_lvl[1L] - 1), 1L] <- 1
  D[n_int - 2L - (n_lvl[3L] - 1):1 + 1, 2L] <- 1
  idx_cat_obs <- c(as.integer(obs$M1),
                   n_lvl[1] + 2L + as.integer(obs$M2))
  idx_not_cat_obs <- setdiff(1:n_int, idx_cat_obs)
  Sig_use <- Sig_int[idx_not_cat_obs, idx_not_cat_obs] +
    D %*% Sig_int[idx_cat_obs, idx_cat_obs] %*% t(D) -
    D %*% Sig_int[idx_cat_obs, idx_not_cat_obs] -
    Sig_int[idx_not_cat_obs, idx_cat_obs] %*% t(D)
  mu_use <- drop(mu_int[idx_not_cat_obs] - D %*% mu_int[idx_cat_obs])

  lw <- rep(-Inf, n_int - 2L)
  up <- numeric(n_int - 2L)

  # handle the binary
  if(obs$B1){
    lw[n_lvl[1L] - 1L + 1L] <- 0
    up[n_lvl[1L] - 1L + 1L] <- Inf
  }

  # handle the ordinal
  lw[n_lvl[1L] - 1L + 2L] <- O_lvls[as.integer(obs$O1)     ]
  up[n_lvl[1L] - 1L + 2L] <- O_lvls[as.integer(obs$O1) + 1L]

  # dput(log(mvtnorm::pmvnorm(lw, up, mu_use, sigma = Sig_use, algorithm = mvtnorm::GenzBretz(1e6, 0, 1e-6))))
  f2 <- -6.8975911813762

  #####
  # permutations should not matter so we check different permutations
  var_idx <- c(list(1:n_lvl[1]), as.list(n_lvl[1] + 1:4),
               list(n_lvl[1] + 4L + 1:n_lvl[3L]))

  # get_perms <- function(x){
  #   stopifnot(is.atomic(x)) # for the matrix call to make sense
  #   out <- as.matrix(expand.grid(
  #     replicate(length(x), x, simplify = FALSE), stringsAsFactors = FALSE))
  #   out[apply(out,1, anyDuplicated) == 0, ]
  # }
  # perms <- apply(get_perms(1:6), 1, list)
  # perms <- sample(perms, 15)
  # dput(perms, control = c("keepInteger"))
  perms <- list(list(c(4L, 6L, 3L, 2L, 5L, 1L)), list(c(3L, 5L, 2L, 4L, 1L, 6L)), list(c(3L, 4L, 5L, 6L, 2L, 1L)), list(c(6L, 5L, 4L, 2L, 1L, 3L)), list(c(5L, 3L, 6L, 1L, 2L, 4L)), list(c(4L, 2L, 6L, 3L, 5L, 1L)), list(c(4L, 3L, 5L, 1L, 6L, 2L)), list(c(2L, 5L, 4L, 1L, 3L, 6L)), list(c(1L, 3L, 4L, 6L, 5L, 2L)), list(c(1L, 6L, 2L, 5L, 4L, 3L)), list(c(1L, 4L, 2L, 3L, 5L, 6L)), list(c(3L, 4L, 2L, 5L, 1L, 6L)), list(c(3L, 2L, 1L, 6L, 5L, 4L)), list(c(6L, 1L, 4L, 2L, 5L, 3L)), list(c(3L, 1L, 6L, 5L, 2L, 4L)))

  for(perm in perms){
    ord <- unlist(perm)
    indicies <- unlist(var_idx[ord])

    # compute the result
    cate_arg <- multinomial[[1L]]
    if(which(ord == 1L) > which(ord == 6L)){
      cate_arg <- cbind(cate_arg[, 2], cate_arg[, 1])
      cate_arg[3, 1] <- which(indicies == var_idx[[6L]][1L]) - 1L
      cate_arg[3, 2] <- which(indicies == var_idx[[1L]][1L]) - 1L
    } else {
      cate_arg[3, 1] <- which(indicies == var_idx[[1L]][1L]) - 1L
      cate_arg[3, 2] <- which(indicies == var_idx[[6L]][1L]) - 1L
    }

    cate_arg <- list(cate_arg)

    cpp_obj <- get_mdgc_log_ml(
      lower = matrix(lower[indicies]),
      upper = matrix(upper[indicies]),
      code  = matrix(code [indicies]),
      multinomial = cate_arg,
      idx_non_zero_mean = integer())

    res <- mdgc_log_ml(cpp_obj, vcov = Sig[indicies, indicies],
                       rel_eps = 1e-5, maxpts = 10000L, mea = numeric())
    expect_equal(res, f1 + f2, tolerance = 1e-4)
  }

  #####
  # missing values
  do_drop <- c(1:n_lvl[1], n_lvl[1] + 2L)
  code[do_drop] <- 1L

  # compute the true value
  is_con <- setdiff(is_con, do_drop)
  is_int <- setdiff(is_int, do_drop)
  # dput(mvtnorm::dmvnorm(obs_x[is_con], mu[is_con], Sig[is_con, is_con, drop = FALSE], log = TRUE))
  f1 <- -1.07915883353727

  Sig_int <- Sig[is_int, is_int] - Sig[is_int, is_con] %*% solve(
    Sig[is_con, is_con, drop = FALSE], Sig[is_con, is_int, drop = FALSE])
  mu_int <- drop(mu[is_int] + Sig[is_int, is_con, drop = FALSE] %*%
                   solve(Sig[is_con, is_con, drop = FALSE],
                         obs_x[is_con, drop = FALSE] - mu[is_con]))
  n_int <- length(mu_int)

  D <- matrix(0., n_int - 1L, 1L)
  D[n_int - 1L - (n_lvl[3L] - 1):1 + 1, 1L] <- 1
  idx_cat_obs <- 2L + as.integer(obs$M2)
  idx_not_cat_obs <- setdiff(1:n_int, idx_cat_obs)
  Sig_use <- Sig_int[idx_not_cat_obs, idx_not_cat_obs] +
    D %*% Sig_int[idx_cat_obs, idx_cat_obs] %*% t(D) -
    D %*% Sig_int[idx_cat_obs, idx_not_cat_obs] -
    Sig_int[idx_not_cat_obs, idx_cat_obs] %*% t(D)
  mu_use <- drop(mu_int[idx_not_cat_obs] - D %*% mu_int[idx_cat_obs])

  lw <- rep(-Inf, n_int - 1L)
  up <- numeric(n_int - 1L)

  # handle the binary
  if(obs$B1){
    lw[1L] <- 0
    up[1L] <- Inf
  }

  # handle the ordinal
  lw[2L] <- O_lvls[as.integer(obs$O1)     ]
  up[2L] <- O_lvls[as.integer(obs$O1) + 1L]

  # dput(log(mvtnorm::pmvnorm(lw, up, mu_use, sigma = Sig_use, algorithm = mvtnorm::GenzBretz(1e6, 0, 1e-6))))
  f2 <- -3.90255011045022

  # check the result
  for(perm in perms){
    ord <- unlist(perm)
    indicies <- unlist(var_idx[ord])

    # compute the result
    cate_arg <- multinomial[[1L]]
    if(which(ord == 1L) > which(ord == 6L)){
      cate_arg <- cbind(cate_arg[, 2], cate_arg[, 1])
      cate_arg[3, 1] <- which(indicies == var_idx[[6L]][1L]) - 1L
      cate_arg[3, 2] <- which(indicies == var_idx[[1L]][1L]) - 1L
    } else {
      cate_arg[3, 1] <- which(indicies == var_idx[[1L]][1L]) - 1L
      cate_arg[3, 2] <- which(indicies == var_idx[[6L]][1L]) - 1L
    }

    cate_arg <- list(cate_arg)

    cpp_obj <- get_mdgc_log_ml(
      lower = matrix(lower[indicies]),
      upper = matrix(upper[indicies]),
      code  = matrix(code [indicies]),
      multinomial = cate_arg, idx_non_zero_mean = integer())

    res <- mdgc_log_ml(cpp_obj, vcov = Sig[indicies, indicies],
                       rel_eps = 1e-5, maxpts = 10000L,
                       mea = numeric())
    expect_equal(res, f1 + f2, tolerance = 1e-4)
  }
})
