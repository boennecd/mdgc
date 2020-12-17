context("testing util functions for multinomial outcomes")

test_that("eval_multinomial_prob gives the correct result", {
  mu <- c(0, c(0.55, -0.97, -2.6, 0.78, 0.27))

  # library(mvtnorm)
  # p <- length(mu)
  # true_vals <- sapply(seq_along(mu), function(idx)
  #   pmvnorm(lower = rep(-Inf, p - 1L), upper = mu[idx] - mu[-idx],
  #           sigma = diag(p - 1) + 1, algorithm = GenzBretz(
  #             maxpts = 1e7, abseps = 0, releps = 1e-8)))
  # dput(true_vals)

  true_vals <- c(0.126705424292455, 0.281696569205857, 0.0215825509358431, 0.000306603937292752, 0.379021425440399, 0.190687426031255)
  output <- sapply(seq_along(mu) - 1, mdgc:::eval_multinomial_prob,
                   means = mu[-1])

  expect_equal(output, true_vals)
})

test_that("multinomial_find_means gives the correct result", {
  p_vals <- c(0.0476190476190476, 0.0952380952380952, 0.142857142857143, 0.19047619047619, 0.238095238095238, 0.285714285714286)
  mu_vals <- mdgc:::multinomial_find_means(p_vals)
  expect_known_value(mu_vals, "multinomial_find_means-means.RDS")

  output <- sapply(seq_along(p_vals) - 1, mdgc:::eval_multinomial_prob,
                   means = mu_vals)

  expect_equal(output, p_vals)
})
