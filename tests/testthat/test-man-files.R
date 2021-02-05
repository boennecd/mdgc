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

test_that("iris examples give the same", {
  # # randomly mask data
  # set.seed(11)
  # masked_data <- iris
  # masked_data[matrix(runif(prod(dim(iris))) < .10, NROW(iris))] <- NA
  # saveRDS(masked_data, "iris-masked.RDS")
  masked_data <- readRDS("iris-masked.RDS")

  obj <- get_mdgc(masked_data)
  ptr <- get_mdgc_log_ml(obj)
  start_vals <- mdgc_start_value(obj)

  # dput(start_vals)
  start_vals_berfore <- structure(c(1, -0.0243963788661465, 0.875816118921817, 0.775978643248133, 5.31132077880212e-19, 0.202567943996632, 0.676494518590346, -0.0243963788661465, 1, -0.18158929404252, -0.108750459396995, 7.22002430143096e-17, -0.565744498651241, -0.195257263139898, 0.875816118921817, -0.181589294042521, 1, 0.832652098028913, -1.37508084725307e-16, 0.129128733940724, 0.845735541031366, 0.775978643248133, -0.108750459396995, 0.832652098028913, 1, -2.55524531818353e-17, 0.0784289851037759, 0.832742423549245, 5.31132077880212e-19, 7.22002430143096e-17, -1.37508084725307e-16, -2.55524531818353e-17, 0, 4.07359095805304e-16, 5.11297740121521e-16, 0.202567943996632, -0.565744498651241, 0.129128733940724, 0.0784289851037759, 4.07359095805304e-16, 1, -0.138440497661705, 0.676494518590346, -0.195257263139898, 0.845735541031366, 0.832742423549245, 5.1129774012152e-16, -0.138440497661705, 1), .Dim = c(7L, 7L))

  expect_equal(start_vals, start_vals_berfore)

  # dput(mdgc_log_ml(ptr, start_vals, obj$means, rel_eps = 1e-5, maxpts = 1000000L))
  truth <- -622.142234979514
  expect_equal(mdgc_log_ml(ptr, start_vals, obj$means),
               truth, tolerance = 1e-4)
  expect_equal(mdgc_log_ml(ptr, start_vals, obj$means, use_aprx = TRUE),
               truth, tolerance = 1e-4)

  # dput(mdgc_log_ml(ptr, start_vals, obj$means, rel_eps = 1e-5, maxpts = 1000000L, comp_derivs = TRUE))
  truth <- structure(-622.142645597963, grad_vcov = structure(c(-5.50752025304173, 15.05467527095, 32.3454449002899, -26.7646832253561, -11.1776068406575, 14.171203475721, -2.99359663506347, 15.05467527095, -21.2939666619648, -15.6146261227432, 7.02228963690592, 21.4732013893823, -14.7665379582272, -6.70666343115514, 32.3454449002899, -15.6146261227432, -121.728737390317, 63.0756920691633, -24.2768074014849, -3.6592268442394, 27.9360342457243, -26.764683225356, 7.02228963690594, 63.0756920691633, -68.2207423432056, -17.2944580410574, -0.413288187640495, 17.7077462286979, -11.1776068406575, 21.4732013893823, -24.2768074014849, -17.2944580410574, -59.77809709442, 17.5365091626148, 42.2415879318052, 14.171203475721, -14.7665379582272, -3.6592268442394, -0.413288187640495, 17.5365091626148, -12.1892006533897, -5.34730850922514, -2.99359663506347, -6.70666343115514, 27.9360342457243, 17.7077462286979, 42.2415879318052, -5.34730850922514, -36.89427942258), .Dim = c(7L, 7L)), grad_mea = structure(c(9.0354807072327, -3.16800905422099), .Dim = 2:1))
  res <- mdgc_log_ml(ptr, start_vals, obj$means, use_aprx = TRUE, comp_derivs = TRUE)
  expect_equal(res, truth, tolerance = 1e-2)

  fit <- mdgc_fit(ptr, start_vals, obj$means, rel_eps = 1e-2, maxpts = 10000L,
                  minvls = 1000L, use_aprx = TRUE, batch_size = 100L, lr = .001,
                  maxit = 100L, n_threads = 2L)
  expect_known_value(fit, "iris-fit.RDS", update = FALSE, tolerance = 1e-2)

  imputed <- mdgc_impute(obj, fit$result$vcov, fit$result$mea, minvls = 1000L,
                         maxit = 10000L, n_threads = 2L, use_aprx = TRUE)
  expect_known_value(imputed[1:15], "iris-imputed.RDS", update = FALSE,
                     tolerance = 1e-2)
})
