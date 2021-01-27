test_that("log_sum_exp() works", {
  set.seed(1)
  x <- rnorm(10)
  x[[2]] <- -Inf

  expect_equal(log(sum(exp(x))), log_sum_exp(x))

  x[[1]] <- NA_real_
  expect_equal(NA_real_, log_sum_exp(x))

  expect_equal(log(sum(exp(x), na.rm = TRUE)), log_sum_exp(x, na.rm = TRUE))
})
