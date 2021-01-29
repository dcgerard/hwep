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
