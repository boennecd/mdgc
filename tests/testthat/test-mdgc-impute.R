context("testing mdgc_impute")

dat <- readRDS("mdgc-fit-test.RDS")

threshold <- function(org_data, imputed){
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

  out_cont <- as.data.frame(
    t(sapply(imputed, function(x) unlist(x[is_cont]))))
  out_cat <- as.data.frame(t(sapply(imputed, function(x)
    sapply(x[is_cat], which.max))))
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

test_that("mdgc_impute gives the same as before and with and without reordering", {
  p <- 5L
  Sig <- diag(p)
  Sig[lower.tri(Sig)] <- Sig[upper.tri(Sig)] <- .5

  dat <- head(dat, 25)
  mdgc_obj <- get_mdgc(dat)
  set.seed(1)
  imp_no_reorder <- mdgc_impute(
    object = mdgc_obj, vcov = Sig, mea = mdgc_obj$means, rel_eps = 1e-5,
    maxit = 1000000L, do_reorder = FALSE, n_threads = 1L, minvls = 100000L)

  imp_reorder <- mdgc_impute(
    object = mdgc_obj, vcov = Sig, mea = mdgc_obj$means, rel_eps = 1e-5,
    maxit = 1000000L, do_reorder = TRUE, n_threads = 1L, minvls = 100000L)

  imp_reorder_two <- mdgc_impute(
    object = mdgc_obj, vcov = Sig, mea = mdgc_obj$means, rel_eps = 1e-5,
    maxit = 1000000L, do_reorder = TRUE, n_threads = 2L, minvls = 100000L)

  expect_equal(imp_no_reorder, imp_reorder, tolerance = 1e-2)
  expect_equal(imp_reorder_two, imp_reorder, tolerance = 1e-2)

  set.seed(1)
  imp_no_reorder_arpx <- mdgc_impute(
    object = mdgc_obj, vcov = Sig, mea = mdgc_obj$means, rel_eps = 1e-5,
    maxit = 1000000L, do_reorder = FALSE, n_threads = 1L, minvls = 100000L,
    use_aprx = TRUE)
  expect_equal(imp_no_reorder, imp_no_reorder_arpx, tolerance = 1e-6)

  # observed values match
  expect_true({
    thresh <- threshold(dat, imp_reorder)
    all((thresh == dat)[!is.na(dat)])
  })
  expect_true({
    thresh <- threshold(dat, imp_no_reorder)
    all((thresh == dat)[!is.na(dat)])
  })

  # compare with saved results
  expect_known_value(imp_reorder, "imp_reorder.RDS",
                     tolerance = 1e-4, update = FALSE)
})

test_that("imputation with ordinal variables yields the correct probabilities", {
  skip_on_cran()

  #####
  # test with two binary variables and a ordinal variable
  probs <- c(.2, .3, .3, .2)
  bs <- qnorm(cumsum(c(0, probs)))
  ps <- c(.6, .4)
  Sig <- matrix(c(1, .5, 0, .5, 1, .5, 0, .5, 1), 3)

  # set.seed(1)
  # y <- c(TRUE, TRUE)
  # Z <- matrix(rnorm(3 * 10000000), ncol = 3) %*% chol(Sig)
  # Z <- Z[if(y[1]) Z[, 1] > -qnorm(ps[1]) else Z[, 1] <=  -qnorm(ps[1]), ]
  # Z <- Z[if(y[2]) Z[, 2] > -qnorm(ps[2]) else Z[, 2] <=  -qnorm(ps[2]), ]
  #
  # ords <- cut(Z[, 3], breaks = bs)
  # dput(table(ords) / length(ords), control = c("keepNA", "keepInteger"))
  truth <- rbind(c(0.244451111011956, 0.337242768857552, 0.287465661112002, 0.130840459018489),
                 c(0.0315335077048528, 0.159682268705399, 0.348076230714025, 0.460707992875723),
                 c(0.330976519542751, 0.361282334883063, 0.235098249872969, 0.0726428957012181),
                 c(0.0825080347201134, 0.244722074602693, 0.35860632402065, 0.314163566656544))

  lower <- rbind(c(-Inf, -Inf,    0, 0),
                 c(-Inf,    0, -Inf, 0),
                 c(NA, NA, NA, NA))
  upper <- rbind(c(0,   0, Inf, Inf),
                 c(0, Inf,   0, Inf),
                 c(NA, NA, NA, NA))
  code  <- rbind(c(2L, 2L, 2L, 2L), c(2L, 2L, 2L, 2L),
                 c(1L, 1L, 1L, 1L))

  margs <- list(structure(list(), borders = c(-Inf, 0, Inf), mean = 10),
                structure(list(), borders = c(-Inf, 0, Inf), mean = -90),
                structure(list(), levels = letters[1:4], borders = bs))

  passed_names <- list(c("FALSE", "TRUE"), c("FALSE", "TRUE"), letters[1:4])
  outer_names <- c("B1", "B2", "O1")
  multinomial <- replicate(4L, matrix(0L, 0, 0), simplify = FALSE)

  set.seed(old_seed <- 1)
  res <- mdgc:::impute(
    lower = lower, upper = upper, code = code, Sigma = Sig, mea = c(qnorm(ps), 0),
    truth = rbind(c(0., 0., 1., 1.),
                  c(0., 1., 0., 1.), c(NA, NA, NA, NA)),
    margs = margs, multinomial = multinomial, rel_eps = 1e-10, abs_eps = -1,
    maxit = 100000L, minvls = 10000L, passed_names = passed_names,
    outer_names = outer_names, n_threads = 1L, do_reorder = TRUE,
    use_aprx = TRUE)

  expect_equal(t(sapply(res, `[[`, "O1")), truth, tolerance = 1e-2,
               check.attributes = FALSE)

  # should give the same
  set.seed(old_seed)
  res_again <- mdgc:::impute(
    lower = lower, upper = upper, code = code, Sigma = Sig, mea = c(qnorm(ps), 0),
    truth = rbind(c(0., 0., 1., 1.),
                  c(0., 1., 0., 1.), c(NA, NA, NA, NA)),
    margs = margs, multinomial = multinomial, rel_eps = 1e-10, abs_eps = -1,
    maxit = 100000L, minvls = 10000L, passed_names = passed_names,
    outer_names = outer_names, n_threads = 1L, do_reorder = TRUE,
    use_aprx = TRUE)
  expect_equal(res_again, res)

  #####
  # test with a continuous variable and two ordinal variables
  p1 <- c(.2, .3, .3, .2)
  b1 <- qnorm(cumsum(c(0, p1)))
  p2 <- c(.2, .3, .5)
  b2 <- qnorm(cumsum(c(0, p2)))
  Sig <- matrix(c(1, .5, .33, .5, 1, .5, .33, .5, 1), 3)
  zs <- c(-1.5, .5, 2.5)

  # truth <- mapply(function(i, j){
  #   Sig_use <- Sig[2:3, 2:3] - Sig[2:3, 1] %o% Sig[1, 2:3] / Sig[1, 1]
  #   mu_use <- Sig[2:3, 1] * zs[i] / Sig[1, 1]
  #   set.seed(1)
  #   Z <- matrix(rnorm(2 * 10000000), ncol = 2) %*% chol(Sig_use)
  #   Z <- sweep(Z, 2L, mu_use, `+`)
  #   Z <- Z[b2[j] < Z[, 2] & Z[, 2] <= b2[j + 1], ]
  #
  #   ords <- cut(Z[, 1], breaks = b1)
  #   table(ords) / length(ords)
  # }, i = rep(1:3, each = 3), j = rep(1:3, 3))
  # dput(truth, control = c("keepNA", "keepInteger"))

  truth <- c(0.633478389442718, 0.282795504521967, 0.0760021134745213, 0.00772399256079393,
             0.440965798379131, 0.377536850586787, 0.157091383862209, 0.0244059671718739,
             0.267778738657526, 0.394576415388093, 0.264320571709374, 0.0733242742450074,
             0.254419357249157, 0.395106914470519, 0.27421070310956, 0.076263025170765,
             0.132552929350638, 0.344658631827303, 0.36453743075579, 0.15825100806627,
             0.0514078518523202, 0.222816532648261, 0.390369571868182, 0.335406043631237,
             0.0489844344428253, 0.224577331287157, 0.399491894819648, 0.326946339450369,
             0.0187626934728291, 0.135225388784816, 0.361196135666787, 0.484815782075568,
             0.00381527516768236, 0.0461127271967623, 0.214518466290286, 0.73555353134527)

  lower <- rbind(rep(NA_real_, 9),
                 rep(NA_real_, 9),
                 rep(head(b2, -1), 3))
  upper <- rbind(rep(zs, each = 3),
                 rep(NA_real_, 9),
                 rep(b2[-1], 3))
  code  <- rbind(rep(0L, 9), rep(1L, 9), rep(2L, 9))

  margs <- list(list(),
                structure(list(), levels = letters[1:4], borders = b1),
                structure(list(), levels = letters[1:3], borders = b2))

  passed_names <- list("", letters[1:4], letters[1:3])
  outer_names <- c("C1", "O1", "O2")
  multinomial <- replicate(9L, matrix(0L, 0, 0), simplify = FALSE)

  set.seed(1)
  res <- mdgc:::impute(
    lower = lower, upper = upper, code = code, Sigma = Sig, mea = numeric(2),
    truth = rbind(upper[1, ], rep(NA, 9), rep(1:3, 3)),
    margs = margs, multinomial = multinomial, rel_eps = 1e-10, abs_eps = -1,
    maxit = 100000L, minvls = 10000L, passed_names = passed_names,
    outer_names = outer_names, n_threads = 1L, do_reorder = TRUE,
    use_aprx = TRUE)

  expect_equal(sapply(res, `[[`, "O1"), truth, tolerance = 1e-3,
               check.attributes = FALSE)
})
