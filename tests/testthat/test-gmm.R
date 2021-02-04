test_that("simplex transform in chisqdiv works", {
  y <- -1
  x <- real_to_simplex(y)
  nvec <- c(1, 2, 3, 1, 2)
  expect_equal(
    obj_reals(y = y, nvec = nvec),
    chisqdiv(x[-1], nvec = nvec)
  )

  ## make sure diploid case works
  y <- NULL
  nvec <- c(1, 2, 3)
  dipout <- obj_reals(y = y, nvec = nvec)
})


test_that("drbounds works", {
  expect_equal(drbounds(4), 1/6)
  expect_equal(drbounds(6), 3/10)
  expect_equal(drbounds(8), c(27/70, 3/140))
  expect_equal(drbounds(10), c(55/126, 5/84))
  expect_equal(drbounds(12), c(285/616, 65/616, 5/1848))
})
