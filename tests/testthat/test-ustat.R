test_that("ginv works", {
  set.seed(98)
  n <- 7
  A <- matrix(rnorm(n * (n-1)), nrow = n)
  omega <- tcrossprod(A)

  gomega <- ginv(omega)

  expect_equal(gomega, t(gomega))

  eout <- eigen(gomega %*% omega, symmetric = TRUE)

  expect_equal(eout$values[1:(n-1)], rep(1, n-1))
  expect_equal(eout$values[n], 0)
})
