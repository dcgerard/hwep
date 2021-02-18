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


test_that("hweustat() is ok", {

  ## Simulation setting
  # ploidy <- 6
  # size <- 10000
  # r <- 0.05
  # alpha <- 0

  thresh <- 100
  nvec <- c(7363L, 2297L, 316L, 24L, 0L, 0L, 0L)
  u1 <- hweustat(nvec = nvec, thresh_mult = Inf, thresh_tot = 0)
  u1 <- hweustat(nvec = nvec)
  u1$chisq_hwe
  u1$df_hwe
  u1$p_hwe

  nvec <- c(7319L, 2366L, 293L, 21L, 1L, 0L, 0L)
  u2 <- hweustat(nvec = nvec)
  u2$chisq_hwe
  u2$df_hwe
  u2$p_hwe

  nvec <- c(7379L, 2287L, 310L, 23L, 1L, 0L, 0L)
  u3 <- hweustat(nvec = nvec)
  u3$chisq_hwe
  u3$p_hwe

  nvec <- c(7301L, 2313L, 344L, 42L, 0L, 0L, 0L)
  u4 <- hweustat(nvec = nvec)
  u4$chisq_hwe
  u4$p_hwe

})
