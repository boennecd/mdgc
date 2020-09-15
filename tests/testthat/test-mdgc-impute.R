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

test_that("mdgc_impute gives the same as before and with and withour reordering", {
  p <- 5L
  Sig <- diag(p)
  Sig[lower.tri(Sig)] <- Sig[upper.tri(Sig)] <- .5

  dat <- head(dat, 25)
  mdgc_obj <- get_mdgc(dat)
  set.seed(1)
  imp_no_reorder <- mdgc_impute(
    object = mdgc_obj, vcov = Sig, rel_eps = 1e-5, maxit = 1000000L,
    do_reorder = FALSE, n_threads = 1L, minvls = 100000L)

  imp_reorder <- mdgc_impute(
    object = mdgc_obj, vcov = Sig, rel_eps = 1e-5, maxit = 1000000L,
    do_reorder = TRUE, n_threads = 1L, minvls = 100000L)

  imp_reorder_two <- mdgc_impute(
    object = mdgc_obj, vcov = Sig, rel_eps = 1e-5, maxit = 1000000L,
    do_reorder = TRUE, n_threads = 2L, minvls = 100000L)

  expect_equal(imp_no_reorder, imp_reorder, tolerance = 1e-2)
  expect_equal(imp_reorder_two, imp_reorder, tolerance = 1e-2)

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
