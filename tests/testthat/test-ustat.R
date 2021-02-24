test_that("ginv works", {
  set.seed(98)
  n <- 7
  A <- matrix(rnorm(n * (n-1)), nrow = n)
  omega <- tcrossprod(A)

  gout <- ginv(omega)
  gomega <- gout$mat

  expect_equal(gomega, t(gomega))

  eout <- eigen(gomega %*% omega, symmetric = TRUE)

  expect_equal(eout$values[1:(n-1)], rep(1, n-1))
  expect_equal(eout$values[n], 0)
})


test_that("hweustat() is ok", {

  nvec <- c(7363L, 2297L, 316L, 24L, 0L, 0L, 0L)
  expect_error(hweustat(nvec = nvec, thresh = 0), NA)
  expect_error(hweustat(nvec = nvec), NA)

  nvec <- c(7319L, 2366L, 293L, 21L, 1L, 0L, 0L)
  expect_error(hweustat(nvec = nvec), NA)

  nvec <- c(7379L, 2287L, 310L, 23L, 1L, 0L, 0L)
  expect_error(hweustat(nvec = nvec), NA)

  nvec <- c(7301L, 2313L, 344L, 42L, 0L, 0L, 0L)
  expect_error(hweustat(nvec = nvec), NA)

  nvec <- c(524L, 359L, 96L, 18L, 3L, 0L, 0L)
  expect_error(hweustat(nvec = nvec), NA)

  nvec <- c(18L, 89L, 225L, 340L, 227L, 78L, 23L)
  expect_error(hweustat(nvec = nvec), NA)

  nmat <- structure(c(84L, 79L, 3L, 5L, 11L, 6L, 2L, 8L, 0L, 2L, 0L, 0L,
0L, 0L, 0L, 0L, 0L, 0L), .Dim = c(2L, 9L))
  expect_error(hwefit(nmat, type = "ustat"), NA)
})
