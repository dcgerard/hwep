##################
## Gradients for objective function
##################

#' Derivative of convolution
#'
#' Derivative of `p*p` with respect to `p`, where `*` is discrete linear
#' convolution. I.e. stats::convolve(p, rev(p), type = "open")
#'
#' @param p The gamete dosage probabilities
#'
#' @author David Gerard
#'
#' @noRd
dg_dp <- function(p) {
  ploidy <- 2 * (length(p) - 1)
  pracma::Toeplitz(a = c(2 * p, rep(0, ploidy / 2)),
                   b = c(2 * p[[1]], rep(0, ploidy / 2)))
}

#' Derivative of total probability law.
#'
#' Derivative of p = B'q with respect to B.
#' Does not depend on B, so only need q (the dosage probabilities).
#'
#' @param q The dosage probabilities.
#'
#' @author David Gerard
#'
#' @noRd
dp_dB <- function(q) {
  ploidy <- length(q) - 1
  kronecker(diag(ploidy / 2 + 1), matrix(q, nrow = 1))
}

#' Derivative of segregation matrix with respect to DR parameters
#'
#' Gradient of \code{\link{gsegmat}()} with respect to alpha.
#' Does not depend on alpha because linear in alpha.
#'
#' @param ploidy The ploidy of the species.
#'
#' @author David Gerard
#'
#' @noRd
dB_dalpha <- function(ploidy) {
  stopifnot(ploidy %% 2 == 0, length(ploidy) == 1, ploidy >= 4)
  ibdr <- floor(ploidy / 4)

  rmat <- expand.grid(ell = 0:ploidy, k = 0:(ploidy / 2))
  ellvec <- rmat[, "ell"]
  kvec <- rmat[, "k"]

  subval <- -exp(
    lchoose(ellvec, kvec) +
      lchoose(ploidy - ellvec, ploidy / 2 - kvec) -
      lchoose(ploidy, ploidy / 2)
  )

  gradmat <- matrix(rep(subval, times = ibdr), nrow = nrow(rmat), ncol = ibdr)

  for (i in seq_len(ibdr)) {
    for (j in 0:i) {
      gradmat[, i] <- gradmat[, i] +
        exp(
          lchoose(ellvec, j) +
          lchoose(ellvec - j, kvec - 2 * j) +
          lchoose(ploidy - ellvec, i - j) +
          lchoose(ploidy - ellvec - (i - j), ploidy / 2 - kvec - 2 * (i - j)) -
          lchoose(ploidy, i) -
          lchoose(ploidy - i, ploidy / 2 - 2 * i)
        )
    }
  }

  return(gradmat)
}

#' Gradient of \code{\link{freqnext}()} with respect to alpha
#'
#' @inheritParams freqnext
#'
#' @author David Gerard
#'
#' @noRd
dfreqnext_dalpha <- function(freq, alpha) {
  ploidy <- length(freq) - 1
  stopifnot(ploidy %% 2 == 0, ploidy >=4, length(ploidy) == 1)
  stopifnot(length(alpha) == floor(ploidy / 4))
  stopifnot(alpha >= 0, sum(alpha) <= 1)
  ibdr <- floor(ploidy / 4)
  B <- gsegmat2(alpha = alpha, ploidy = ploidy)
  p <- c(t(freq) %*% B)

  dg_dp(p = p) %*% dp_dB(q = freq) %*% dB_dalpha(ploidy = ploidy)
}

#' Derivative of pearsondiv() with respect to alpha
#'
#' @param nvec The counts
#' @param alpha The DR parameters
#'
#' @author David Gerard
#'
#' @noRd
dpearsondiv_dalpha <- function(nvec, alpha) {
  n <- sum(nvec)
  q <- nvec / n
  fq <- freqnext(freq = q, alpha = alpha)
  grad_f <- dfreqnext_dalpha(freq = q, alpha = alpha)

  n * colSums((-2 * (q - fq) / fq - (q - fq)^2 / fq^2) * grad_f)
}

#' Derivative of neymandiv() with respect to alpha
#'
#' @param nvec The counts
#' @param alpha The DR parameters
#'
#' @author David Gerard
#'
#' @noRd
dneymandiv_dalpha <- function(nvec, alpha) {
  n <- sum(nvec)
  q <- nvec / n
  fq <- freqnext(freq = q, alpha = alpha)
  grad_f <- dfreqnext_dalpha(freq = q, alpha = alpha)

  n * colSums(-2 * (q - fq) / q  * grad_f)
}

#' Derivative of gdiv() with respect to alpha
#'
#' @param nvec The counts
#' @param alpha The DR parameters
#'
#' @author David Gerard
#'
#' @noRd
dgdiv_dalpha <- function(nvec, alpha) {
  n <- sum(nvec)
  q <- nvec / n
  fq <- freqnext(freq = q, alpha = alpha)
  grad_f <- dfreqnext_dalpha(freq = q, alpha = alpha)

  n * colSums(-2 * q / fq * grad_f)
}
