test_that("Huang theoretical genotypes are correct", {
  set.seed(7)
  alpha <- runif(1)
  r <- runif(1)

  tout <- theofreq(alpha = alpha, r = r, ploidy = 4)
  hout <- hwefreq(r = r, alpha = alpha, ploidy = 4, more = TRUE)

  expect_equal(
    tout$q,
    hout$q,
    tolerance = 10^-4
  )

  tout <- theofreq(alpha = alpha, r = r, ploidy = 6)
  hout <- hwefreq(r = r, alpha = alpha, ploidy = 6, more = TRUE)

  expect_equal(
    tout$q,
    hout$q,
    tolerance = 10^-4
  )


  alpha <- runif(2, max = 0.5)

  tout <- theofreq(alpha = alpha, r = r, ploidy = 8)
  hout <- hwefreq(r = r, alpha = alpha, ploidy = 8, more = TRUE)

  expect_equal(
    tout$q,
    hout$q,
    tolerance = 10^-4
  )

  tout <- theofreq(alpha = alpha, r = r, ploidy = 10)
  hout <- hwefreq(r = r, alpha = alpha, ploidy = 10, more = TRUE)

  expect_equal(
    tout$q,
    hout$q,
    tolerance = 10^-4
  )
})
