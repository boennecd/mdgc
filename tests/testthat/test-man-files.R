context("testing examples in man files")

test_that("mdgc examples gives the same", {
  skip_if_not_installed("catdata")

  data(retinopathy, envir = environment(), package = "catdata")
  retinopathy$RET <- as.ordered(retinopathy$RET)
  retinopathy$SM <- as.logical(retinopathy$SM)

  set.seed(28325145)
  truth <- retinopathy
  for(i in seq_along(retinopathy))
    retinopathy[[i]][runif(NROW(retinopathy)) < .3] <- NA

  impu <- mdgc(retinopathy, lr = 1e-3, maxit = 25L, batch_size = 25L,
               rel_eps = 1e-3, maxpts = 5000L, n_threads = 1L,
               method = "svrg")

  expect_known_value(impu, "mdgc-man-test.RDS", update = FALSE,
                     tolerance = 1e-5)

  impu <- mdgc(retinopathy, maxit = 25L, rel_eps = 1e-3, maxpts = 5000L,
               n_threads = 1L, method = "aug_Lagran")

  expect_known_value(impu, "mdgc-man-Lagran-test.RDS", update = FALSE,
                     tolerance = 1e-5)
})
