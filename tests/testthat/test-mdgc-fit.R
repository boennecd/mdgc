context("testing mdgc_fit")

# sim_dat <- function(n, p = 4, n_lvls = 5L){
#   # get the covariance matrix
#   Sb <- diag(p)
#   Sb[lower.tri(Sb)] <- Sb[upper.tri(Sb)] <- .5
#   Sb <- Sb / p / 5
#   Sig <- cov2cor(drop(rWishart(1L, 5L * p, Sb)))
#
#   # draw the observations
#   truth <- matrix(rnorm(n * p), n) %*% chol(Sig)
#
#   # determine the type
#   type <- rep(1:3, each = floor((p + 3 - 1) / 3))[1:p]
#   is_con <- type == 1L
#   is_bin <- type == 2L
#   is_ord <- type == 3L
#
#   # sample which are masked data
#   is_mask <- matrix(runif(n * p) < .3, n)
#
#   # create observed data
#   truth_obs <- data.frame(truth)
#   truth_obs[, is_con] <- qexp(pnorm(as.matrix(truth_obs[, is_con])))
#
#   bs_bin <- c(-Inf, 0., Inf)
#   truth_obs[, is_bin] <- truth_obs[, is_bin] > bs_bin[2]
#
#   bs_ord <- qnorm(seq(0, 1, length.out = n_lvls + 1L))
#   truth_obs[, is_ord] <- as.integer(cut(truth[, is_ord], breaks = bs_ord))
#   for(i in which(is_ord)){
#     truth_obs[, i] <- ordered(truth_obs[, i])
#     levels(truth_obs[, i]) <-
#       LETTERS[seq_len(length(unique(truth_obs[, i])))]
#   }
#
#   # mask the data
#   seen_obs <- truth_obs
#   seen_obs[is_mask] <- NA
#
#   list(truth = truth, truth_obs = truth_obs, seen_obs = seen_obs,
#        Sigma = Sig)
# }
#
#
# set.seed(1)
# dat <- sim_dat(n = 200, p = 5)
# dat <- dat$seen_obs
# saveRDS(dat, "mdgc-fit-test.RDS")
dat <- readRDS("mdgc-fit-test.RDS")

test_that("ADAM gives the same", {
  mdgc_obj <- get_mdgc(dat)
  log_ml_ptr <- get_mdgc_log_ml(mdgc_obj)
  start_val <- mdgc_start_value(mdgc_obj)

  set.seed(1L)
  fit_adam <- mdgc_fit(
    ptr = log_ml_ptr, vcov = start_val, mea = mdgc_obj$means,
    n_threads = 1L, lr = 1e-2, maxit = 5L, batch_size = 10L,
    method = "adam", minvls = 1000L)
  expect_known_value(fit_adam, file = "mdgc_fit-ADAM-res.RDS",
                     tolerance = 1e-6, update = FALSE)

  set.seed(1L)
  fit_adam_2 <- mdgc_fit(
    ptr = log_ml_ptr, vcov = start_val, mea = mdgc_obj$means,
    n_threads = 2L, lr = 1e-2, maxit = 5L, batch_size = 10L,
    method = "adam", minvls = 1000L)
  expect_equal(fit_adam, fit_adam_2, tolerance = 1e-2)

  set.seed(1L)
  fit_adam_3 <- mdgc_fit(
    ptr = log_ml_ptr, vcov = start_val, mea = mdgc_obj$means,
    n_threads = 1L, lr = 1e-2, maxit = 5L, batch_size = 10L,
    method = "adam", minvls = 1000L, use_aprx = TRUE)
  expect_equal(fit_adam, fit_adam_3, tolerance = 1e-6)
})

test_that("svrg gives the same", {
  mdgc_obj <- get_mdgc(dat)
  log_ml_ptr <- get_mdgc_log_ml(mdgc_obj)
  start_val <- mdgc_start_value(mdgc_obj)

  set.seed(1L)
  fit_svrg <- mdgc_fit(
    ptr = log_ml_ptr, vcov = start_val, mea = mdgc_obj$means,
    n_threads = 1L, lr = 1e-2, maxit = 5L, batch_size = 10L,
    method = "svrg", minvls = 1000L)
  expect_known_value(fit_svrg, file = "mdgc_fit-SVRG-res.RDS",
                     tolerance = 1e-6, udpate = FALSE)

  set.seed(1L)
  fit_svrg_2 <- mdgc_fit(
    ptr = log_ml_ptr, vcov = start_val, mea = mdgc_obj$means,
    n_threads = 2L, lr = 1e-2, maxit = 5L, batch_size = 10L,
    method = "svrg", minvls = 1000L)
  expect_equal(fit_svrg, fit_svrg_2, tolerance = 1e-2)

  set.seed(1L)
  fit_svrg_3 <- mdgc_fit(
    ptr = log_ml_ptr, vcov = start_val, mea = mdgc_obj$means,
    n_threads = 1L, lr = 1e-2, maxit = 5L, batch_size = 10L,
    method = "svrg", minvls = 1000L, use_aprx = TRUE)
  expect_equal(fit_svrg, fit_svrg_3, tolerance = 1e-6)
})
