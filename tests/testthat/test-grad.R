test_that("dg_dp works", {

  set.seed(18)
  p <- stats::runif(5)
  p <- p / sum(p)

  fn <- function(p) {
    stats::convolve(p, rev(p), type = "open")
  }

  expect_equal(
    dg_dp(p = p),
    pracma::jacobian(f = fn, x0 = p),
    tolerance = 10^-6
  )

  # microbenchmark::microbenchmark(
  #   dg_dp(p = p),
  #   pracma::jacobian(f = fn, x0 = p)
  # )

})

test_that("dp_dB works", {
  set.seed(19)
  q <- stats::runif(11)

  B <- matrix(stats::runif(66), nrow = 11)
  B <- B / rowSums(B)

  fn <- function(B) {
    ## q is loaded because this is R
    t(B) %*% q
  }

  expect_equal(
    dp_dB(q = q),
    pracma::jacobian(f = fn, x0 = B),
    tolerance = 10^-6
  )

  # microbenchmark::microbenchmark(
  #   dp_dB(q = q),
  #   pracma::jacobian(f = fn, x0 = B)
  # )
})

test_that("dB_dalpha works", {
  alpha <- c(0.3, 0.1)
  ploidy <- 10

  expect_equal(
    dB_dalpha(ploidy = ploidy),
    pracma::jacobian(f = gsegmat2, x0 = alpha, ploidy = ploidy),
    tolerance = 10^-6
  )

  # microbenchmark::microbenchmark(
  #   dB_dalpha(ploidy = ploidy),
  #   pracma::jacobian(f = gsegmat2, x0 = alpha, ploidy = ploidy)
  # )

})

test_that("dfreqnext_dalpha works", {
  set.seed(17)
  alpha <- c(0.3, 0.1)
  ploidy <- 10
  freq <- runif(ploidy + 1)
  freq <- freq / sum(freq)

  expect_equal(
    dfreqnext_dalpha(freq = freq, alpha = alpha),
    pracma::jacobian(f = freqnext, x0 = alpha, freq = freq),
    tolerance = 10^-6
  )

  # microbenchmark::microbenchmark(
  #   dfreqnext_dalpha(freq = freq, alpha = alpha),
  #   pracma::jacobian(f = freqnext, x0 = alpha, freq = freq)
  # )
})

test_that("dpearsondiv_dalpha works", {
  set.seed(21)
  ploidy <- 10
  alpha <- c(0.2, 0.14)
  q <- runif(ploidy + 1)
  q <- q / sum(q)
  nvec <- c(rmultinom(n = 1, size = 100, prob = q))

  expect_equal(
    dpearsondiv_dalpha(nvec = nvec, alpha = alpha),
    c(pracma::jacobian(f = pearsondiv, x0 = alpha, nvec = nvec)),
    tolerance = 10^-6
  )

  expect_equal(
    dpearsondiv_dalpha(nvec = nvec, alpha = alpha, ngen = 2),
    c(pracma::jacobian(f = pearsondiv, x0 = alpha, nvec = nvec, ngen = 2)),
    tolerance = 10^-6
  )

  # microbenchmark::microbenchmark(
  #   dpearsondiv_dalpha(nvec = nvec, alpha = alpha),
  #   c(pracma::jacobian(f = pearsondiv, x0 = alpha, nvec = nvec))
  # )
})


test_that("dneymandiv_dalpha works", {
  set.seed(22)
  ploidy <- 10
  alpha <- c(0.2, 0.14)
  q <- runif(ploidy + 1)
  q <- q / sum(q)
  nvec <- c(rmultinom(n = 1, size = 100, prob = q))

  expect_equal(
    dneymandiv_dalpha(nvec = nvec, alpha = alpha),
    c(pracma::jacobian(f = neymandiv, x0 = alpha, nvec = nvec)),
    tolerance = 10^-6
  )

  expect_equal(
    dneymandiv_dalpha(nvec = nvec, alpha = alpha, ngen = 2),
    c(pracma::jacobian(f = neymandiv, x0 = alpha, nvec = nvec, ngen = 2)),
    tolerance = 10^-6
  )

  # microbenchmark::microbenchmark(
  #   dneymandiv_dalpha(nvec = nvec, alpha = alpha),
  #   c(pracma::jacobian(f = neymandiv, x0 = alpha, nvec = nvec))
  # )
})

test_that("dgdiv_dalpha works", {
  set.seed(23)
  ploidy <- 10
  alpha <- c(0.2, 0.14)
  q <- runif(ploidy + 1)
  q <- q / sum(q)
  nvec <- c(rmultinom(n = 1, size = 100, prob = q))

  expect_equal(
    dgdiv_dalpha(nvec = nvec, alpha = alpha),
    c(pracma::jacobian(f = gdiv, x0 = alpha, nvec = nvec)),
    tolerance = 10^-6
  )

  expect_equal(
    dgdiv_dalpha(nvec = nvec, alpha = alpha, ngen = 2),
    c(pracma::jacobian(f = gdiv, x0 = alpha, nvec = nvec, ngen = 2)),
    tolerance = 10^-6
  )

  # microbenchmark::microbenchmark(
  #   dgdiv_dalpha(nvec = nvec, alpha = alpha),
  #   c(pracma::jacobian(f = gdiv, x0 = alpha, nvec = nvec))
  # )
})


test_that("dfreqnext_dalpha_ngen works", {
  set.seed(1)
  ploidy <- 4
  freq <- runif(ploidy + 1)
  freq <- freq / sum(freq)
  alpha <- 0.1
  ngen <- 3

  expect_equal(
    pracma::jacobian(f = freqnext_ngen, x0 = alpha, freq = freq, ngen = ngen),
    dfreqnext_dalpha_ngen(freq = freq, alpha = alpha, ngen = ngen)
  )

  # microbenchmark::microbenchmark(
  #   pracma::jacobian(f = freqnext_ngen, x0 = alpha, freq = freq, ngen = ngen),
  #   dfreqnext_dalpha_ngen(freq = freq, alpha = alpha, ngen = ngen)
  # )

})


test_that("duobj_dalpha works", {
  set.seed(75)
  ploidy <- 10
  alpha <- c(0.2, 0.14)
  q <- runif(ploidy + 1)
  q <- q / sum(q)
  nvec <- c(rmultinom(n = 1, size = 100, prob = q))
  which_keep <- sample(x = c(TRUE, FALSE), size = ploidy + 1, replace = TRUE)


  expect_equal(
    duobj_dalpha(nvec = nvec, alpha = alpha, omega = NULL, which_keep = NULL),
    c(pracma::jacobian(f = uobj, x0 = alpha, nvec = nvec, omega = NULL, which_keep = NULL))
  )

  # microbenchmark::microbenchmark(
  #   duobj_dalpha(nvec = nvec, alpha = alpha, omega = NULL, which_keep = NULL),
  #   c(pracma::jacobian(f = uobj, x0 = alpha, nvec = nvec, omega = NULL, which_keep = NULL))
  # )

  expect_equal(
    duobj_dalpha(nvec = nvec, alpha = alpha, omega = NULL, which_keep = which_keep),
    c(pracma::jacobian(f = uobj, x0 = alpha, nvec = nvec, omega = NULL, which_keep = which_keep))
  )

  A <- matrix(rnorm((ploidy + 1) * (ploidy - 1)), ncol = ploidy - 1)
  omega <- tcrossprod(A)
  expect_equal(
    duobj_dalpha(nvec = nvec, alpha = alpha, omega = omega, which_keep = NULL),
    c(pracma::jacobian(f = uobj, x0 = alpha, nvec = nvec, omega = omega, which_keep = NULL))
  )

  numkeep <- sum(which_keep)
  A <- matrix(rnorm((numkeep + 1) * (numkeep - 1)), ncol = numkeep - 1)
  omega <- tcrossprod(A)
  expect_equal(
    duobj_dalpha(nvec = nvec, alpha = alpha, omega = omega, which_keep = which_keep),
    c(pracma::jacobian(f = uobj, x0 = alpha, nvec = nvec, omega = omega, which_keep = which_keep))
  )

  # microbenchmark::microbenchmark(
  #   duobj_dalpha(nvec = nvec, alpha = alpha, omega = omega, which_keep = which_keep),
  #   c(pracma::jacobian(f = uobj, x0 = alpha, nvec = nvec, omega = omega, which_keep = which_keep))
  # )
})

