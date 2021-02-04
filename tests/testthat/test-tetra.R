test_that("gamete and alpha correspondence holds", {
  a <- 0.11
  r <- 0.78
  p <- gam_from_a(a = a, r = r)

  expect_equal(sum(p), 1)

  ar <- a_from_gam(p = p)

  expect_equal(ar[["a"]], a)
  expect_equal(ar[["r"]], r)
})


test_that("sim and closed form are same", {
  sim <- p_from_alpha(alpha = 0.1, p = 0.1, ploidy = 4)
  cf <- gam_from_a(a = 0.1, r = 0.1)
  expect_equal(sim, cf, tolerance = 10^-4, ignore_attr = TRUE)
})
