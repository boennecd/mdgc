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
  ptr <- mdgc:::get_log_lm_terms_cpp(lower = dat$lower, upper = dat$upper,
                                     code = dat$code)

  lcov_to_mat <- function(par){
    p <- (sqrt(8 * length(par) + 1) - 1) / 2
    L <- matrix(0, p, p)
    L[lower.tri(L, TRUE)] <- par
    L[upper.tri(L)] <- t(L)[upper.tri(L)]
    L
  }

  log_ml <- function(vcov_log_chol, releps = 1e-5, n_threads = 1L,
                     comp_derivs = FALSE, seed = NULL){
    if(!is.null(seed))
      set.seed(seed)
    Arg <- lcov_to_mat(vcov_log_chol)

    mdgc:::eval_log_lm_terms(
      ptr = ptr, vcov = Arg, indices = 0:(NCOL(dat$lower) - 1L),
      maxpts = 1000000L, abseps = -1, releps = releps, n_threads = n_threads,
      comp_derivs = comp_derivs)
  }

  set.seed(1)
  val <- log_ml(dat$Sigma[lower.tri(dat$Sigma, TRUE)])
  # dput(val)
  expect_equal(val, -45.8852616647973, tolerance = 1e-3)

  val <- log_ml(dat$Sigma[lower.tri(dat$Sigma, TRUE)], n_threads = 2)
  expect_equal(val, -45.8852616647973, tolerance = 1e-3)

  grad <- log_ml(dat$Sigma[lower.tri(dat$Sigma, TRUE)], comp_derivs = TRUE)
  # truth <- numDeriv::jacobian(log_ml, dat$Sigma[lower.tri(dat$Sigma, TRUE)], seed = 1L)
  # dput(truth)
  true_grad <- structure(c(
    4.44392200186265, 8.26645723733806, 0.0770691610564516,
    1.28391853303438, -10.3954140181992, -0.487395978873399,
    2.03205557640943, -0.105889973274042, -5.218467741155,
    -4.20994886759362, -2.94592203652891, 4.71790766215746,
    -1.62255015484385, 2.75034658786021, 0.703564051111291),
    .Dim = c(1L, 15L))

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
  expect_equal(c(attr(grad, "grad")) %*% jac, true_grad,
               tolerance = 1e-3)

  grad <- log_ml(dat$Sigma[lower.tri(dat$Sigma, TRUE)], comp_derivs = TRUE,
                 n_threads = 2L)
  expect_equal(c(attr(grad, "grad")) %*% jac, true_grad,
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
  ptr <- mdgc:::get_log_lm_terms_cpp(lower = dat$lower, upper = dat$upper,
                                     code = dat$code)

  lcov_to_mat <- function(par){
    p <- (sqrt(8 * length(par) + 1) - 1) / 2
    L <- matrix(0, p, p)
    L[lower.tri(L, TRUE)] <- par
    L[upper.tri(L)] <- t(L)[upper.tri(L)]
    L
  }

  log_ml <- function(vcov_log_chol, releps = 1e-5, n_threads = 1L,
                     comp_derivs = FALSE, seed = NULL){
    if(!is.null(seed))
      set.seed(seed)
    Arg <- lcov_to_mat(vcov_log_chol)

    mdgc:::eval_log_lm_terms(
      ptr = ptr, vcov = Arg, indices = 0:(NCOL(dat$lower) - 1L),
      maxpts = 1000000L, abseps = -1, releps = releps, n_threads = n_threads,
      comp_derivs = comp_derivs)
  }

  set.seed(1)
  val <- log_ml(dat$Sigma[lower.tri(dat$Sigma, TRUE)])
  # dput(val)
  expect_equal(val, -50.633131628025, tolerance = 1e-3)

  val <- log_ml(dat$Sigma[lower.tri(dat$Sigma, TRUE)], n_threads = 2)
  expect_equal(val, -50.633131628025, tolerance = 1e-3)

  grad <- log_ml(dat$Sigma[lower.tri(dat$Sigma, TRUE)], comp_derivs = TRUE)
  # truth <- numDeriv::jacobian(log_ml, dat$Sigma[lower.tri(dat$Sigma, TRUE)], seed = 1L)
  # dput(truth)
  true_grad <- structure(c(
    2.18828846268184, 1.14819853573079, -0.833095815058999,
    -1.64775972436418, -3.02462052022083, 2.38038644284807, 4.02535849827106,
    -0.6823072284458, -1.52397119919305, -5.83385872168417, -0.526004633937689,
    -0.772721040056757, -2.15498992108154, 4.76760990982455, -0.156298143778767,
    0.415973653298547, -0.190202356593625, 2.98216021926465, 1.03082094344473,
    0.816776256035582, -1.55995256546344),
    .Dim = c(1L, 21L))

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
  expect_equal(c(attr(grad, "grad")) %*% jac, true_grad,
               tolerance = 1e-3)

  grad <- log_ml(dat$Sigma[lower.tri(dat$Sigma, TRUE)], comp_derivs = TRUE,
                 n_threads = 2L)
  expect_equal(c(attr(grad, "grad")) %*% jac, true_grad,
               tolerance = 1e-3)
})
