test_that("hwemom() works", {
  set.seed(11)
  alpha <- stats::runif(1)
  nvec <- round(hwefreq(r = 0.5,
                        alpha = alpha,
                        ploidy = 6,
                        niter = 100,
                        tol = -Inf) * 10000000)
  hout <- hwemom(nvec)
  expect_equal(hout$alpha, alpha, tolerance = 10^-4)
})
