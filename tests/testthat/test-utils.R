test_that("log_sum_exp() works", {
  set.seed(1)
  x <- rnorm(10)
  x[[2]] <- -Inf

  expect_equal(log(sum(exp(x))), log_sum_exp(x))

  x[[1]] <- NA_real_
  expect_equal(NA_real_, log_sum_exp(x))

  expect_equal(log(sum(exp(x), na.rm = TRUE)), log_sum_exp(x, na.rm = TRUE))
})

test_that("logit/expit are inverses", {
  set.seed(1)
  pvec <- stats::runif(100)
  expect_equal(expit(logit(pvec)), pvec)
  expect_equal(logit(expit(pvec)), pvec)
})

test_that("simplex transform works", {
  set.seed(1)

  ## Check they are inverses ----
  x <- stats::runif(5)
  x <- x / sum(x)
  expect_equal(real_to_simplex(simplex_to_real(x)), x)

  y <- stats::rnorm(5)
  x <- real_to_simplex(y)
  expect_equal(sum(x), 1)
  expect_equal(simplex_to_real(x), y)

  ## Check two and one case works
  x <- c(0.3, 0.7)
  expect_equal(real_to_simplex(simplex_to_real(x)), x)

})
