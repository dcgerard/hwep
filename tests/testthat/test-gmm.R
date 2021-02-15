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

test_that("pearson and neyman work", {
  nvec <- 1:7
  alpha <- 0.1
  expect_equal(
    chisqdiv(nvec = nvec, alpha = alpha, denom = "expected"),
    pearsondiv(nvec = nvec, alpha = alpha)
  )
  expect_equal(
    chisqdiv(nvec = nvec, alpha = alpha, denom = "observed"),
    neymandiv(nvec = nvec, alpha = alpha)
  )
})

test_that("gdiv works on diploids", {
  expect_error(gdiv(nvec = 0:2, alpha = NULL), NA)
})


test_that("hwemom gives no errors on hard data", {
  nvec <- c(0,0,0,0,0,1,1,1,99999)
  expect_error(hwemom(nvec = nvec, ngen = 5, obj = "g"), NA)
  expect_error(hwemom(nvec = nvec, ngen = 1, obj = "g"), NA)

  expect_error(hwemom(nvec = nvec, ngen = 5, obj = "pearson"), NA)
  expect_error(hwemom(nvec = nvec, ngen = 1, obj = "pearson"), NA)

  expect_error(hwemom(nvec = nvec, ngen = 5, obj = "neyman"), NA)
  expect_error(hwemom(nvec = nvec, ngen = 1, obj = "neyman"), NA)
})

test_that("hwetetra gives no errors on hard data", {
  nvec <- c(0,0,0,99,1)
  expect_error(hwetetra(nvec = nvec), NA)
})
