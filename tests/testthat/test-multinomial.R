context("testing util functions for multinomial outcomes")

test_that("eval_multinomial_prob and eval_multinomial_prob_gr gives the correct result", {
  skip_on_os("solaris")
  mu <- c(0, c(0.55, -0.97, -2.6, 0.78, 0.27))

  # library(mvtnorm)
  # p <- length(mu)
  # true_vals <- sapply(seq_along(mu), function(idx){
  #   Sig <- diag(p - 1)
  #   if(idx > 1){
  #     Sig[1, 1] <- 0
  #     Sig <- Sig + 1
  #   }
  #
  #   pmvnorm(lower = rep(-Inf, p - 1L), upper = mu[idx] - mu[-idx],
  #           sigma = Sig, algorithm = GenzBretz(
  #             maxpts = 1e7, abseps = -1, releps = 1e-8))
  # })
  # dput(true_vals)

  # GL_func <- function(mu, icase){
  #   gw <- pracma::gaussLaguerre(100L)
  #   if(icase == 1L)
  #     return(prod(pnorm(-mu[-1L])))
  #
  #   fval <- sapply(gw$x, function(x)
  #     x + dnorm(x, mu[icase], log = TRUE) +
  #       sum(pnorm(x, mu[-c(1L, icase)], log.p = TRUE)))
  #   sum(exp(fval) * gw$w)
  # }
  # all.equal(true_vals, sapply(seq_along(mu), GL_func, mu = mu))

  true_vals <- c(0.0207080219507315, 0.315745869382916, 0.0258216323301918,
                 0.000365873445112673, 0.421195801033678, 0.216162802754332)
  output <- sapply(seq_along(mu) - 1, mdgc:::eval_multinomial_prob,
                   means = mu[-1])

  expect_equal(output, true_vals)

  # num_grads <- sapply(seq_along(mu), function(icase){
  #   out <- numDeriv::grad(GL_func, mu, icase = icase)
  #   c(GL_func(mu, icase), out[-1])
  # })
  # dput(num_grads)

  num_grads <- structure(c(0.0207080219507315, -0.0243910445019943, -0.00618843591731022,  -0.000282593647183361, -0.0279954191264598, -0.0202388339533849,  0.315745870711258, 0.28579011711952, -0.0149322047858618, -0.00028346459410506,  -0.154793407366863, -0.0913899958753983, 0.0258216330490677,  -0.0149322047887554, 0.0507897558101159, -5.12626363552435e-05,  -0.0181724482578314, -0.011445404209538, 0.000365873459721941,  -0.000283464594023833, -5.12626362501077e-05, 0.00117844428747922,  -0.000337014964323351, -0.000224108445706233, 0.421195796521844,  -0.154793407364547, -0.0181724482607582, -0.000337014964066006,  0.315117364964993, -0.113819075241132, 0.216162804307017, -0.0913899958697994,  -0.0114454042105515, -0.000224108445544844, -0.113819075245085,  0.23711741772536), .Dim = c(6L, 6L))
  output <- sapply(seq_along(mu) - 1, function(icase) {
    out <- mdgc:::eval_multinomial_prob_gr(icase, mu[-1])
    c(attr(out, "prob"), out)
  })

  expect_equal(output, num_grads)
})

test_that("multinomial_find_means gives the correct result", {
  skip_on_os("solaris")
  p_vals <- c(0.0476190476190476, 0.0952380952380952, 0.142857142857143, 0.19047619047619, 0.238095238095238, 0.285714285714286)
  mu_vals <- mdgc:::multinomial_find_means(p_vals)
  expect_known_value(mu_vals, "multinomial_find_means-means.RDS")

  output <- sapply(seq_along(p_vals) - 1, mdgc:::eval_multinomial_prob,
                   means = mu_vals)

  expect_equal(output, p_vals)
})
