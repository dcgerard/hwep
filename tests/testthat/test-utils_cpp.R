test_that("log-sum-exp works", {
  set.seed(12)
  expect_equal(
    log_sum_exp_cpp(c(-Inf, -Inf, -Inf)),
    -Inf
  )

  expect_equal(
    log_sum_exp_2_cpp(-Inf, -Inf),
    -Inf
  )

  x <- log(runif(20))
  expect_equal(
    log_sum_exp(x),
    log_sum_exp_cpp(x)
  )

  x[[1]] <- -Inf
  expect_equal(
    log_sum_exp(x),
    log_sum_exp_cpp(x)
  )

  expect_equal(
    log_sum_exp_2_cpp(-10, -Inf),
    -10
  )

  expect_equal(
    log_sum_exp_2_cpp(-10, -11),
    log(exp(-10) + exp(-11))
  )
})

test_that("conv_cpp is same as stats::convolve()", {
  set.seed(99)
  x <- runif(5)
  y <- runif(11)

  expect_equal(
    stats::convolve(x, rev(y), type = "open"),
    conv_cpp(x, y)
  )

  expect_equal(
    stats::convolve(y, rev(x), type = "open"),
    conv_cpp(y, x)
  )

  expect_equal(
    stats::convolve(1, rev(x), type = "open"),
    conv_cpp(1, x)
  )

  expect_equal(
    stats::convolve(1, 2, type = "open"),
    conv_cpp(1, 2)
  )

  # bench::mark(
  #   stats::convolve(x, rev(y), type = "open"),
  #   conv_cpp(x, y)
  # )

})



