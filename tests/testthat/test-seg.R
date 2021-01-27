test_that("Segregation probabilities are same as known proportions in tetraploids", {
  set.seed(1)
  alpha <- runif(n = 1, min = 0, max = 1)

  expect_equal(
    dgamete(x = 0:2, alpha = alpha, G = 0, ploidy = 4),
    c(1, 0, 0)
  )

  expect_equal(
    dgamete(x = 0:2, alpha = alpha, G = 1, ploidy = 4),
    c(2 + alpha, 2 * (1 - alpha), alpha) / 4
  )

  expect_equal(
    dgamete(x = 0:2, alpha = alpha, G = 2, ploidy = 4),
    c(1 + 2 * alpha, 4 * (1 - alpha), 1 + 2 * alpha) / 6
  )

  expect_equal(
    dgamete(x = 0:2, alpha = alpha, G = 3, ploidy = 4),
    c(alpha, 2 * (1 - alpha), 2 + alpha) / 4
  )

  expect_equal(
    dgamete(x = 0:2, alpha = alpha, G = 4, ploidy = 4),
    c(0, 0, 1)
  )

})


test_that("Segregation probabilities are same as known proportions in hexaploids", {
  set.seed(2)
  alpha <- runif(n = 1, min = 0, max = 1)

  expect_equal(
    dgamete(x = 0:3, alpha = alpha, G = 0, ploidy = 6),
    c(1, 0, 0, 0)
  )

  expect_equal(
    dgamete(x = 0:3, alpha = alpha, G = 1, ploidy = 6),
    c(3 + alpha, 3 - 2 * alpha, alpha, 0) / 6
  )

  expect_equal(
    dgamete(x = 0:3, alpha = alpha, G = 2, ploidy = 6),
    c(3 + 3 * alpha, 9 - 5 * alpha, 3 + alpha, alpha) / 15
  )

  expect_equal(
    dgamete(x = 0:3, alpha = alpha, G = 3, ploidy = 6),
    c(1 + 3 * alpha, 9 - 3 * alpha, 9 - 3 * alpha, 1 + 3 * alpha) / 20
  )

  expect_equal(
    dgamete(x = 0:3, alpha = alpha, G = 4, ploidy = 6),
    c(alpha, 3 + alpha, 9 - 5 * alpha, 3 + 3 * alpha) / 15
  )

  expect_equal(
    dgamete(x = 0:3, alpha = alpha, G = 5, ploidy = 6),
    c(0, alpha, 3 -  2 * alpha, 3 + alpha) / 6
  )

  expect_equal(
    dgamete(x = 0:3, alpha = alpha, G = 6, ploidy = 6),
    c(0, 0, 0, 1)
  )

})

