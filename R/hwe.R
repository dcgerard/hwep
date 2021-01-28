##########################
## Functions to simulate HWE genotype proportions
##########################

#' Update genotype frequencies after one generation
#'
#' After one generation of random mating, update the genotype
#' frequencies.
#'
#' @inheritParams dgamete
#' @param freq The current genotype frequencies. This should be a
#'     vector of length K+1, where K is the ploidy of the species.
#'     \code{freq[i]} could contain the proportion of individuals
#'     that have genotype \code{i-1}.
#' @param segarray The output of \code{\link{zsegarray}()}. We will calculate
#'     it if \code{segarray = NULL}. It is just an option so we don't
#'     have to recalculate it very often if \code{alpha} is known.
#'
#' @return A vector of length \code{lenght(freq)} that contains the
#'     updated genotype frequencies after one generation of random mating.
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' freq <- c(0.5, 0, 0, 0, 0.5)
#' freqnext(freq = freq, alpha = 0)
#'
freqnext <- function(freq, alpha, segarray = NULL) {
  ploidy <- length(freq) - 1
  stopifnot(ploidy %% 2 == 0)
  stopifnot(length(alpha) == floor(ploidy / 4))
  stopifnot(alpha >= 0, sum(alpha) <= 1)

  if (is.null(segarray)) {
    segarray <- zsegarray(alpha = alpha, ploidy = ploidy)
  } else {
    stopifnot(dim(segarray) == rep(ploidy + 1, 3))
  }

  mararray <- sweep(x = segarray,
                    MARGIN = c(1, 2),
                    STATS = tcrossprod(freq),
                    FUN = `*`)

  freqnew <- apply(mararray, 3, sum)
  freqnew <- freqnew / sum(freqnew) ## to resolve numerical issues

  return(freqnew)
}

#' Generate HWE genotype frequencies
#'
#' Generate genotype frequencies under Hardy-Weinberg equilibrium
#' given the allele frequency of the reference allele (\code{p}),
#' the double reduction parameter (\code{alpha}), and the ploidy
#' of the species (\code{ploidy}).
#'
#' If \code{alpha} is not all 0, then this function repeatedly
#' applies \code{\link{freqnext}()} to simulate genotype frequencies
#' under HWE. Otherwise, it uses \code{\link[stats]{dbinom}()}.
#'
#' @inheritParams dgamete
#' @param p The allele frequency of the reference allele.
#' @param niter The maximum number of iterations to simulate.
#' @param tol The stopping criterion on the Chi-square divergence between
#'     old and new genotype frequencies.
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' freq1 <- hwefreq(p = 0.5, alpha = 0, ploidy = 4)
#' freq2 <- hwefreq(p = 0.5, alpha = 1/4, ploidy = 4)
#'
#' plot(x = 0:4,
#'      y = freq1,
#'      type = "h",
#'      ylim = c(0, 0.4),
#'      xlab = "dosage",
#'      ylab = "Pr(dosage)")
#' plot(x = 0:4,
#'      y = freq2,
#'      type = "h",
#'      ylim = c(0, 0.4),
#'      xlab = "dosage",
#'      ylab = "Pr(dosage)")
#'
hwefreq <- function(p, alpha, ploidy, niter = 100, tol = 10^-4) {
  stopifnot(length(p) == 1L, length(ploidy) == 1L, length(niter) == 1L)
  stopifnot(ploidy %% 2 == 0)
  stopifnot(ploidy > 1)
  stopifnot(length(alpha) == floor(ploidy / 4))
  stopifnot(alpha >= 0, sum(alpha) <= 1)
  stopifnot(p >= 0, p <= 1)
  stopifnot(niter >= 1)

  ## Return theoretical result when no double reduction ----
  if (all(alpha < sqrt(.Machine$double.eps))) {
    freq <- stats::dbinom(x = 0:ploidy, size = ploidy, prob = p)
    return(freq)
  }

  ## Iterate freqnext() if double reduction ----
  segarray <- zsegarray(alpha = alpha, ploidy = ploidy)

  freq <- c(1 - p, rep(0, length.out = ploidy - 1), p)
  i <- 1
  err <- Inf
  while (i < niter && err > tol) {
    oldfreq <- freq
    freq <- freqnext(freq = freq, alpha = alpha, segarray = segarray)
    i <- i + 1
    err <- sum((oldfreq - freq) ^ 2 / freq) * (ploidy + 1)
  }

  return(freq)
}

