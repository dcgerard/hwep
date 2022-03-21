test_that("aggfun works", {
  expect_equal(choose_agg(x = c(10, 11), thresh = 5),
               c(TRUE, TRUE))

  expect_equal(choose_agg(x = c(10, 11, 0), thresh = 5),
               c(FALSE, TRUE, FALSE))

  expect_equal(choose_agg(x = c(10, 11, 3, 3), thresh = 5),
               c(TRUE, TRUE, FALSE, FALSE))
})
