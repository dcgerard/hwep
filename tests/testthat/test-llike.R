test_that("rmem works", {
  set.seed(1)
  pvec <- stats::runif(6)
  pvec <- pvec / sum(pvec)

  ## Genotype frequencies
  qvec <- stats::convolve(pvec, rev(pvec), type = "open")

  ## Generate data
  nvec <- round(100000000 * qvec)

  ## Estimate pvec
  pest <- rmem(nvec = nvec, tol = 10^-6)

  expect_equal(pest, pvec, tolerance = 10^-4)
})

test_that("hwelike works", {
  # hard nvec
  nvec <- c(421L, 390L, 159L, 27L, 3L, 0L, 0L, 0L, 0L)
  expect_error(hwelike(nvec = nvec), NA)

  nvec <- c(522L, 356L, 110L, 10L, 2L, 0L, 0L)
  expect_error(hwelike(nvec = nvec), NA)

  nvec <- c(48194L, 31292L, 14514L, 4696L, 1067L, 205L, 30L, 2L, 0L)
  expect_error(hwelike(nvec = nvec), NA)

  nvec <- c(25L, 0L, 0L, 0L, 0L)
  expect_error(rmlike(nvec = nvec), NA)
})

test_that("hwenodr works", {
  nvec <- c(421L, 390L, 159L, 27L, 3L, 0L, 0L, 0L, 0L)
  expect_error(hwenodr(nvec = nvec), NA)

  nvec <- c(522L, 356L, 110L, 10L, 2L, 0L, 0L)
  expect_error(hwenodr(nvec = nvec), NA)

  nvec <- c(48194L, 31292L, 14514L, 4696L, 1067L, 205L, 30L, 2L, 0L)
  expect_error(hwenodr(nvec = nvec), NA)
})
