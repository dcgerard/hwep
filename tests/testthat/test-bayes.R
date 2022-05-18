test_that("dirichlet is correct", {

  set.seed(1)
  itermax <- 10000
  dirmat <- matrix(NA_real_, nrow = itermax, ncol = 2)
  for (i in seq_len(itermax)) {
    dirmat[i, ] <- rdirichlet1(alpha = c(1, 2))
  }

  bvec <- rbeta(n = itermax, shape1 = 1, shape2 = 2)

  # qqplot(dirmat[, 1], bvec)
  # abline(0, 1, lty = 2, col = 2)

  expect_true(stats::ks.test(dirmat[, 1], bvec)$p.value > 0.01)
})


test_that("samp_gametes is correct", {
  set.seed(1)
  ploidy <- 8
  x <- 1:(ploidy + 1)
  p <- stats::runif(ploidy / 2 + 1)
  p <- p / sum(p)
  y <- samp_gametes(x = x, p = p)
  expect_equal(sum(y), sum(x) * 2)

  ## 10 microseconds on average
  # bench::mark(
  #   samp_gametes(x = x, p = p)
  # )
})

test_that("dmultinom_cpp is correct", {
  expect_equal(
    dmultinom_cpp(x = c(0, 1, 30, 11), p = c(0.1, 0.5, 0.1, 0.3), lg = TRUE),
    dmultinom(x = c(0, 1, 30, 11), prob = c(0.1, 0.5, 0.1, 0.3), log = TRUE)
  )

  expect_equal(
    dmultinom_cpp(x = c(0, 0, 0, 0), p = c(0, 0.6, 0.1, 0.3), lg = TRUE),
    dmultinom(x = c(0, 0, 0, 0), prob = c(0, 0.6, 0.1, 0.3), log = TRUE)
  )

  expect_equal(
    dmultinom_cpp(x = c(0, 1, 0, 0), p = c(0, 0.6, 0.1, 0.3), lg = TRUE),
    dmultinom(x = c(0, 1, 0, 0), prob = c(0, 0.6, 0.1, 0.3), log = TRUE)
  )

  expect_equal(
    dmultinom_cpp(x = c(1, 0, 0, 0), p = c(0, 0.6, 0.1, 0.3), lg = TRUE),
    dmultinom(x = c(1, 0, 0, 0), prob = c(0, 0.6, 0.1, 0.3), log = TRUE)
  )

  expect_equal(
    dmultinom_cpp(x = c(1, 1, 0, 0), p = c(0, 0.6, 0.1, 0.3), lg = TRUE),
    dmultinom(x = c(1, 1, 0, 0), prob = c(0, 0.6, 0.1, 0.3), log = TRUE)
  )

})

test_that("dirichlet pdf is correct", {
  # set.seed(1)
  # p <- runif(10)
  # p <- p / sum(p)
  # alpha <- abs(rnorm(10)) / 10
  # expect_equal(
  #   ddirichlet(x = p, alpha = alpha, lg = FALSE),
  #   gtools::ddirichlet(x = p, alpha = alpha)
  # )

  p <- c(0.6, 0.4)
  alpha <- c(2, 4)
  expect_equal(
    ddirichlet(x = p, alpha = alpha, lg = TRUE),
    stats::dbeta(x = p[[1]], shape1 = alpha[[1]], shape2 = alpha[[2]], log = TRUE)
  )

  expect_error(ddirichlet(p, alpha = c(-1, 10)))
  expect_error(ddirichlet(c(1, 2), alpha = c(1, 10)))
  expect_error(ddirichlet(c(0.4, 0.6), alpha = c(1)))
})


test_that("marginal likelihood is 1 with no data", {
  expect_equal(
    gibbs_known(x = rep(0, 3), alpha = rep(1, 2), lg = TRUE)$mx,
    0
  )

  expect_equal(
    gibbs_known(x = rep(0, 5), alpha = rep(1, 3), lg = TRUE)$mx,
    0
  )

  expect_equal(
    gibbs_known(x = rep(0, 7), alpha = rep(1, 4), lg = TRUE)$mx,
    0
  )

  expect_equal(
    gibbs_known(x = rep(0, 9), alpha = rep(1, 5), lg = TRUE)$mx,
    0
  )

  ## median of 50 ms.
  # bdf <- bench::mark(
  #   gibbs_known(x = c(1, 2, 5, 1, 2), alpha = c(1, 1, 1), lg = TRUE)
  # )
  # bdf$median
})


test_that("plq and gl_alt_marg are the same", {
  gl <- matrix(-abs(rnorm(110)), nrow = 10)
  beta <- 1:11

  expect_equal(
    plq(gl = gl, beta = beta, lg = TRUE),
    gl_alt_marg(gl = gl, beta = beta, lg = TRUE)
  )

  ## plq() is about 10 times faster
  # temp <- bench::mark(
  #   plq(gl = gl, beta = beta, lg = TRUE),
  #   gl_alt_marg(gl = gl, beta = beta, lg = TRUE)
  # )

})
