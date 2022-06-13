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

test_that("random mating same for tetraploids", {
  nvec <- c(1, 4, 2, 3, 1)

  tout <- hwetetra(nvec = nvec)
  rout <- rmlike(nvec = nvec, nstarts = 0)

  expect_equal(rout$chisq_rm, tout$chisq_rm)
})


test_that("early return is same as late return", {
  nvec <- c(1, 0, 0, 0, 0)
  nvec_full <- c(1, 1, 1, 1, 1)

  hout1 <- hwetetra(nvec = nvec, more = TRUE)
  hout2 <- hwetetra(nvec = nvec_full, more = TRUE)
  expect_equal(
    names(hout1),
    names(hout2)
  )

  hout1 <- hwetetra(nvec = nvec, more = FALSE)
  hout2 <- hwetetra(nvec = nvec_full, more = FALSE)
  expect_equal(
    names(hout1),
    names(hout2)
  )
})
