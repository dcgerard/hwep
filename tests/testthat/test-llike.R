test_that("rmem works", {
  set.seed(1)
  pvec <- stats::runif(6)
  pvec <- pvec / sum(pvec)

  ## Genotype frequencies
  qvec <- stats::convolve(pvec, rev(pvec), type = "open")

  ## Generate data
  nvec <- round(100000000 * qvec)

  ## Estimate pvec
  pest <- rmem(nvec = nvec, tol = 10^-6)

  expect_equal(pest, pvec, tolerance = 10^-4)
})
