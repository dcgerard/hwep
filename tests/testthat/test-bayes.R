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


test_that("gibbs sampler for known genotypes give good estimates", {
  set.seed(1)
  ploidy <- 8

  ## Simulate under the null
  p <- stats::runif(ploidy / 2 + 1)
  p <- p / sum(p)
  q <- stats::convolve(p, rev(p), type = "open")

  ## See BF increase
  nvec <- c(stats::rmultinom(n = 1, size = 100000, prob = q))

  gout <- gibbs_known(x = nvec, alpha = rep(1, ploidy / 2 + 1), more = TRUE, lg = TRUE)

  ## plot(cumsum(gout$post) / seq_along(gout$post), type = "l")
  # plot(gout$p[, 1], type = "l")
  # plot(gout$p[, 2], type = "l")
  # plot(gout$p[, 3], type = "l")

  expect_equal(colMeans(gout$p), p, tolerance = 0.01)
})

test_that("gibbs sampler using genotype log-likelihoods give good estimates", {
  set.seed(1)
  ploidy <- 8

  ## Simulate under the null
  p <- stats::runif(ploidy / 2 + 1)
  p <- p / sum(p)
  q <- stats::convolve(p, rev(p), type = "open")

  ## See BF increase
  nvec <- c(stats::rmultinom(n = 1, size = 10, prob = q))
  gl <- simgl(nvec = nvec, sig = 0.01)

  gout <- gibbs_gl(gl = gl, alpha = rep(1, ploidy / 2 + 1), more = TRUE, lg = TRUE)

  # gout$mx
  # plot(cumsum(gout$post) / seq_along(gout$post), type = "l")
  # plot(gout$p[, 1], type = "l")
  # plot(gout$p[, 2], type = "l")
  # plot(gout$p[, 3], type = "l")
  # plot(gout$p[, 4], type = "l")
  # plot(gout$p[, 5], type = "l")
  #
  # table(round(colMeans(gout$z)))
  # table(apply(gl, 1, which.max) - 1)
  # nvec


  # plot(apply(gl, 1, which.max) - 1, colMeans(gout$z))
  # abline(0, 1)

  expect_equal(apply(gl, 1, which.max) - 1, colMeans(gout$z))

  # colMeans(gout$p)
  # p
})


test_that("alt gibbs sampler using genotype log-likelihoods give good estimates", {
  skip("alt not ready for primetime")
  set.seed(1)
  ploidy <- 8

  ## Simulate under the null
  q <- stats::runif(ploidy + 1)
  q <- q / sum(q)

  ## See BF increase
  nvec <- c(stats::rmultinom(n = 1, size = 10000, prob = q))
  gl <- simgl(nvec = nvec, sig = 0.2)

  gout <- gibbs_gl_alt(gl = gl, beta = rep(1, ploidy + 1), more = TRUE, lg = TRUE)

  i <- 1000
  gout$x[i, ] / sum(gout$x[i, ])
  round(gout$q[i, ], digits = 3)

  gout$mx
  plot(cumsum(gout$post) / seq_along(gout$post), type = "l")
  plot(gout$q[, 1], type = "l")
  plot(gout$q[, 2], type = "l")
  plot(gout$q[, 3], type = "l")
  plot(gout$q[, 4], type = "l")
  plot(gout$q[, 5], type = "l")

  table(round(colMeans(gout$z)))
  table(apply(gl, 1, which.max) - 1)
  nvec

  plot(apply(gl, 1, which.max) - 1, colMeans(gout$z))
  abline(0, 1)

  expect_equal(apply(gl, 1, which.max) - 1, colMeans(gout$z))

  plot(colMeans(gout$q), q)
  abline(0, 1)
})


test_that("sample_z() works", {
  set.seed(1)
  ploidy <- 8

  ## Simulate under the null
  q <- stats::runif(ploidy + 1)
  q <- q / sum(q)

  ## See BF increase
  nvec <- c(stats::rmultinom(n = 1, size = 10, prob = q))
  gl <- simgl(nvec = nvec, sig = 0.2)

  sample_z(gl = gl, q = q)

  zmat <- t(replicate(n = 100000, expr = {
    sample_z(gl = gl, q = q)
  }))

  probmat <- gl + rep(log(q), each = nrow(gl))
  probmat <- exp(probmat - apply(probmat, 1, log_sum_exp))
  rowSums(probmat)

  expect_equal(colMeans(zmat == 0), probmat[, 1], tolerance = 0.05)
  expect_equal(colMeans(zmat == 1), probmat[, 2], tolerance = 0.05)
  expect_equal(colMeans(zmat == 2), probmat[, 3], tolerance = 0.05)
  expect_equal(colMeans(zmat == 3), probmat[, 4], tolerance = 0.05)
  expect_equal(colMeans(zmat == 4), probmat[, 5], tolerance = 0.05)
  expect_equal(colMeans(zmat == 5), probmat[, 6], tolerance = 0.05)
  expect_equal(colMeans(zmat == 6), probmat[, 7], tolerance = 0.05)
  expect_equal(colMeans(zmat == 7), probmat[, 8], tolerance = 0.05)
  expect_equal(colMeans(zmat == 8), probmat[, 9], tolerance = 0.05)
})


test_that("posterior matrix is filled properly", {
  set.seed(1)
  n <- 10
  ploidy <- 4
  gl <- matrix(-abs(rnorm(n * (ploidy + 1))), nrow = n)
  q <- runif(ploidy + 1)
  q <- q / sum(q)

  postmat <- matrix(NA_real_, nrow = n, ncol = ploidy + 1)

  mod_postmat(postmat = postmat, gl = gl, q = q)

  p2 <- exp(gl) * rep(q, each = n)
  p2 <- p2 / rowSums(p2)

  expect_equal(postmat, p2)

  # tdf <- bench::mark(
  #   mod_postmat(postmat = postmat, gl = gl, q = q),
  #   {
  #   p2 <- exp(gl) * rep(q, each = n)
  #   p2 <- p2 / rowSums(p2)
  #   },
  #   check = FALSE
  # )
  # View(tdf)

})


test_that("sample_int returns proper proportions", {
  set.seed(1)
  q <- runif(4)
  q <- q / sum(q)

  qhat <- table(replicate(n = 100000, expr = sample_int(probs = q)))
  expect_true(chisq.test(x = qhat, p = q)$p.value > 0.01)
  qhat <- prop.table(qhat)
  expect_equal(q, qhat, tolerance = 0.01, ignore_attr = TRUE)
})


test_that("dpairs from Levene (1949) is correct", {
  A <- matrix(c(1, 0, 0,
                0, 0, 0,
                0, 0, 0),
              ncol = 3,
              byrow = TRUE)
  y <- c(2, 0, 0)
  expect_equal(dpairs(A, y), 1)

  A <- matrix(c(0, 1, 0,
                0, 0, 1,
                0, 0, 0),
              ncol = 3,
              byrow = TRUE)
  y <- c(1, 2, 1)
  expect_equal(dpairs(A, y), 2/3)

  A <- matrix(c(0, 0, 1,
                0, 1, 0,
                0, 0, 0),
              ncol = 3,
              byrow = TRUE)
  y <- c(1, 2, 1)
  expect_equal(dpairs(A, y), 1/3)
})

test_that("gam_from_pairs() works", {
  A <- matrix(c(0, 0, 1,
                0, 1, 0,
                0, 0, 0),
              ncol = 3,
              byrow = TRUE)
  y <- c(1, 2, 1)
  expect_equal(gam_from_pairs(A), y)
})

test_that("exact marginal calculation is correct", {
  x <- c(4, 3, 0, 1, 2)
  alpha <- 1:3

  expect_equal(
    gibbs_known(x = x, alpha = alpha, lg = TRUE)$mx,
    tetra_rm_marg(x = x, alpha = alpha, lg = TRUE)
  )

  x <- c(4, 3, 10, 1, 2)
  alpha <- 1:3

  expect_equal(
    gibbs_known(x = x, alpha = alpha, lg = TRUE)$mx,
    tetra_rm_marg(x = x, alpha = alpha, lg = TRUE),
    tolerance = 0.01
  )

  # temp <- bench::mark(
  #   gibbs_known(x = x, alpha = alpha, lg = TRUE)$mx,
  #   tetra_rm_marg(x = x, alpha = alpha, lg = TRUE),
  #   check = FALSE
  # )


   x <- c(4, 3, 0, 0, 0, 1, 2)
  alpha <- 1:4

  expect_equal(
    gibbs_known(x = x, alpha = alpha, lg = TRUE)$mx,
    hexa_rm_marg(x = x, alpha = alpha, lg = TRUE)
  )

  x <- c(4, 3, 13, 14, 13, 1, 2)
  alpha <- 1:4

  expect_equal(
    gibbs_known(x = x, alpha = alpha, lg = TRUE)$mx,
    hexa_rm_marg(x = x, alpha = alpha, lg = TRUE),
    tolerance = 0.01
  )

  # temp <- bench::mark(
  #   gibbs_known(x = x, alpha = alpha, lg = TRUE)$mx,
  #   hexa_rm_marg(x = x, alpha = alpha, lg = TRUE),
  #   check = FALSE
  # )

})


test_that("dirichlet-multinomial sums to 1", {
  tupmat <- all_multinom(n = 4, k = 3)
  probvec <- rep(NA_real_, length.out = nrow(tupmat))
  alpha <- 1:3

  for (i in seq_along(probvec)) {
    probvec[[i]] <- ddirmult(x = tupmat[i, ], alpha = alpha, lg = TRUE)
  }

  expect_equal(log_sum_exp(probvec), 0)
})


test_that("beta_from_alpha same as conc_default", {
  expect_equal(
    beta_from_alpha(alpha = rep(1, 3)),
    conc_default(ploidy = 4)$beta
  )

  expect_equal(
    beta_from_alpha(alpha = rep(1, 4)),
    conc_default(ploidy = 6)$beta
  )

  expect_equal(
    beta_from_alpha(alpha = rep(1, 5)),
    conc_default(ploidy = 8)$beta
  )

  expect_equal(
    beta_from_alpha(alpha = rep(1, 6)),
    conc_default(ploidy = 10)$beta
  )
})

test_that("beta_from_alpha more generally gives same BF", {
  alpha <- c(1, 2, 4)
  beta <- beta_from_alpha(alpha = alpha)

  expect_equal(
    rmbayes(nvec = c(1, 0, 0, 0, 0), lg = TRUE, alpha = alpha, beta = beta),
    0
  )

  expect_equal(
    rmbayes(nvec = c(0, 1, 0, 0, 0), lg = TRUE, alpha = alpha, beta = beta),
    0
  )

  expect_equal(
    rmbayes(nvec = c(0, 0, 1, 0, 0), lg = TRUE, alpha = alpha, beta = beta),
    0
  )

  expect_equal(
    rmbayes(nvec = c(0, 0, 0, 1, 0), lg = TRUE, alpha = alpha, beta = beta),
    0
  )

  expect_equal(
    rmbayes(nvec = c(0, 0, 0, 0, 1), lg = TRUE, alpha = alpha, beta = beta),
    0
  )
})

test_that("rmbayesgl() works without errors", {
  load("./glshir.RData")
  beta <- conc_default(ploidy = 6)$beta
  expect_true(!is.na(rmbayesgl(gl = glshir, nburn = 2, niter = 5)))
})
