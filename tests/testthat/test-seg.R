test_that("Segregation probabilities are same as known proportions in tetraploids", {
  set.seed(1)
  alpha <- runif(n = 1, min = 0, max = 1)

  expect_equal(
    dgamete(x = 0:2, alpha = alpha, G = 0, ploidy = 4),
    c(1, 0, 0)
  )

  expect_equal(
    dgamete(x = 0:2, alpha = alpha, G = 1, ploidy = 4),
    c(2 + alpha, 2 * (1 - alpha), alpha) / 4
  )

  expect_equal(
    dgamete(x = 0:2, alpha = alpha, G = 2, ploidy = 4),
    c(1 + 2 * alpha, 4 * (1 - alpha), 1 + 2 * alpha) / 6
  )

  expect_equal(
    dgamete(x = 0:2, alpha = alpha, G = 3, ploidy = 4),
    c(alpha, 2 * (1 - alpha), 2 + alpha) / 4
  )

  expect_equal(
    dgamete(x = 0:2, alpha = alpha, G = 4, ploidy = 4),
    c(0, 0, 1)
  )

})


test_that("Segregation probabilities are same as known proportions in hexaploids", {
  set.seed(2)
  alpha <- runif(n = 1, min = 0, max = 1)

  expect_equal(
    dgamete(x = 0:3, alpha = alpha, G = 0, ploidy = 6),
    c(1, 0, 0, 0)
  )

  expect_equal(
    dgamete(x = 0:3, alpha = alpha, G = 1, ploidy = 6),
    c(3 + alpha, 3 - 2 * alpha, alpha, 0) / 6
  )

  expect_equal(
    dgamete(x = 0:3, alpha = alpha, G = 2, ploidy = 6),
    c(3 + 3 * alpha, 9 - 5 * alpha, 3 + alpha, alpha) / 15
  )

  expect_equal(
    dgamete(x = 0:3, alpha = alpha, G = 3, ploidy = 6),
    c(1 + 3 * alpha, 9 - 3 * alpha, 9 - 3 * alpha, 1 + 3 * alpha) / 20
  )

  expect_equal(
    dgamete(x = 0:3, alpha = alpha, G = 4, ploidy = 6),
    c(alpha, 3 + alpha, 9 - 5 * alpha, 3 + 3 * alpha) / 15
  )

  expect_equal(
    dgamete(x = 0:3, alpha = alpha, G = 5, ploidy = 6),
    c(0, alpha, 3 -  2 * alpha, 3 + alpha) / 6
  )

  expect_equal(
    dgamete(x = 0:3, alpha = alpha, G = 6, ploidy = 6),
    c(0, 0, 0, 1)
  )

})

test_that("gsegmat() and gsegmat2() are same as special cases", {
  alpha <- 0.112

  expect_equal(
    gsegmat(alpha = NULL, ploidy = 2),
    gsegmat_diploid()
  )

  expect_equal(
    gsegmat(alpha = alpha, ploidy = 4),
    gsegmat_tetraploid(alpha = alpha)
  )

  expect_equal(
    gsegmat2(alpha = alpha, ploidy = 4),
    gsegmat_tetraploid(alpha = alpha),
    ignore_attr = TRUE
  )

  expect_equal(
    gsegmat(alpha = alpha, ploidy = 6),
    gsegmat_hexaploid(alpha = alpha)
  )

  expect_equal(
    gsegmat2(alpha = alpha, ploidy = 6),
    gsegmat_hexaploid(alpha = alpha),
    ignore_attr = TRUE
  )

  # microbenchmark::microbenchmark(
  #   gsegmat2(alpha = alpha, ploidy = 6),
  #   gsegmat_hexaploid(alpha = alpha)
  # )
})


test_that("zsegarray() works", {

  segarray <- zsegarray(alpha = c(0.14, 0.1), ploidy = 10)

  expect_equal(
    segarray[1, 4, 1:6],
    dgamete(x = 0:5, alpha = c(0.14, 0.1), G = 3, ploidy = 10),
    ignore_attr = TRUE
  )

})

test_that("hwe probs are asymptotically correct", {
  freq1 <- hwefreq(r = 0.3, alpha = c(0, 0), ploidy = 10, tol = sqrt(.Machine$double.eps))
  freq2 <- stats::dbinom(x = 0:10, size = 10, prob = 0.3)
  expect_equal(freq1, freq2, tolerance = 10^-4)

  freq3 <- hwefreq(r = 0.3, alpha = c(0.5, 0.1), ploidy = 10, tol = sqrt(.Machine$double.eps))

})

test_that("freqnext() and freqnext2() give same results", {
  set.seed(8)
  ploidy <- 8
  freq <- runif(ploidy + 1)
  freq <- freq / sum(freq)
  alpha <- c(0.3, 0.1)

  expect_equal(
    freqnext(freq = freq, alpha = alpha),
    freqnext2(freq = freq, alpha = alpha),
    ignore_attr = TRUE
  )

  # microbenchmark::microbenchmark(
  #   freqnext(freq = freq, alpha = alpha),
  #   freqnext2(freq = freq, alpha = alpha)
  # )
})

test_that("gsegmat() and gsegmat2() give same results", {
  set.seed(1)
  ploidy <- 2
  alpha <- NULL
  expect_equal(
    gsegmat(alpha = alpha, ploidy = ploidy),
    gsegmat2(alpha = alpha, ploidy = ploidy),
    tolerance = 10^-6,
    ignore_attr = TRUE
  )

  ploidy <- 4
  alpha <- stats::runif(n = floor(ploidy / 4))
  alpha <- alpha / sum(alpha) * 0.5
  expect_equal(
    gsegmat(alpha = alpha, ploidy = ploidy),
    gsegmat2(alpha = alpha, ploidy = ploidy),
    tolerance = 10^-6,
    ignore_attr = TRUE
  )

  ploidy <- 6
  alpha <- stats::runif(n = floor(ploidy / 4))
  alpha <- alpha / sum(alpha) * 0.5
  expect_equal(
    gsegmat(alpha = alpha, ploidy = ploidy),
    gsegmat2(alpha = alpha, ploidy = ploidy),
    tolerance = 10^-6,
    ignore_attr = TRUE
  )

  ploidy <- 8
  alpha <- stats::runif(n = floor(ploidy / 4))
  alpha <- alpha / sum(alpha) * 0.5
  expect_equal(
    gsegmat(alpha = alpha, ploidy = ploidy),
    gsegmat2(alpha = alpha, ploidy = ploidy),
    tolerance = 10^-6,
    ignore_attr = TRUE
  )

  ploidy <- 10
  alpha <- stats::runif(n = floor(ploidy / 4))
  alpha <- alpha / sum(alpha) * 0.5
  expect_equal(
    gsegmat(alpha = alpha, ploidy = ploidy),
    gsegmat2(alpha = alpha, ploidy = ploidy),
    tolerance = 10^-6,
    ignore_attr = TRUE
  )

  ploidy <- 12
  alpha <- stats::runif(n = floor(ploidy / 4))
  alpha <- alpha / sum(alpha) * 0.5
  expect_equal(
    gsegmat(alpha = alpha, ploidy = ploidy),
    gsegmat2(alpha = alpha, ploidy = ploidy),
    tolerance = 10^-6,
    ignore_attr = TRUE
  )

  # microbenchmark::microbenchmark(
  #   gsegmat(alpha = alpha, ploidy = ploidy),
  #   gsegmat2(alpha = alpha, ploidy = ploidy)
  # )
})

test_that("freqnext gives parental gametes", {
  set.seed(8)
  ploidy <- 10
  alpha <- c(0.3, 0.1)
  q <- runif(ploidy + 1)
  q <- q / sum(q)

  fout <- freqnext(freq = q, alpha = alpha, more = TRUE)

  expect_equal(
    fout$q,
    stats::convolve(fout$p, rev(fout$p), type = "open")
  )
})

test_that("hwefreq gives parental gametes", {
  set.seed(8)
  ploidy <- 10
  alpha <- c(0.3, 0.1)
  r <- 0.1

  hout <- hwefreq(r = r, alpha = alpha, ploidy = ploidy, more = TRUE)
  expect_equal(
    hout$q,
    stats::convolve(hout$p, rev(hout$p), type = "open")
  )

  alpha <- c(0, 0)
  hout <- hwefreq(r = r, alpha = alpha, ploidy = ploidy, more = TRUE)
  expect_equal(
    hout$q,
    stats::convolve(hout$p, rev(hout$p), type = "open")
  )

  ploidy <- 4
  alpha <- 1/4
  hout <- hwefreq(r = r, alpha = alpha, ploidy = ploidy, more = TRUE)
  expect_equal(
    hout$q,
    stats::convolve(hout$p, rev(hout$p), type = "open")
  )

  ploidy <- 4
  alpha <- 1/4
  hout <- hwefreq(r = r, alpha = alpha, ploidy = ploidy, more = TRUE, niter = 1)
  expect_equal(
    hout$q,
    stats::convolve(hout$p, rev(hout$p), type = "open")
  )
})


test_that("freqnext and freqnext_ngen give same results at ngen = 1", {
  set.seed(7)
  ploidy <- 4
  freq <- runif(ploidy + 1)
  freq <- freq / sum(freq)
  alpha <- 0.1

  expect_equal(
    freqnext(freq = freq, alpha = alpha),
    freqnext_ngen(freq = freq, alpha = alpha, ngen = 1)
  )
})


test_that("parental indices are transposable", {
  A <- zsegarray(alpha = 0.1, ploidy = 6)
  expect_equal(A, aperm(A, c(2, 1, 3)),
               tolerance = 10^-6,
               ignore_attr = TRUE)
})


test_that("gsegmat_symb works", {

  ploidy <- 4
  a <- 0.1
  b <- 0.9
  symbout <- gsegmat_symb(ploidy = ploidy, out = "exp")
  segmat <- gsegmat(alpha = b, ploidy = ploidy)

  for (i in 0:ploidy) {
    for (j in 0:(ploidy / 2)) {
      expect_equal(
        eval(symbout[i+1, j+1]),
        segmat[i+1, j+1]
      )
    }
  }

})
