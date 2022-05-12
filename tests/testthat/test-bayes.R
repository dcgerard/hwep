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
