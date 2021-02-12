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
#' @examples
#' freq <- c(0.5, 0, 0, 0, 0.5)
#' freqnext2(freq = freq, alpha = 0)
#'
#' @noRd
freqnext2 <- function(freq, alpha, segarray = NULL) {
  ploidy <- length(freq) - 1
  stopifnot(ploidy %% 2 == 0, ploidy > 0)
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
#' @param segmat You can provide your own segregation matrix.
#'     \code{segmat[i, j]} is the probability that a parent with
#'     dosage \code{i-1} produces a gamete with dosage \code{j-1}.
#' @param more A logical. Should we return more output (\code{TRUE}) or
#'     less (\code{FALSE}). See the Value section for details.
#'
#' @return If \code{more = FALSE}, then returns a vector of length
#'     \code{lenght(freq)} that contains the updated genotype frequencies
#'     after one generation of random mating. If \code{more = TRUE}, then
#'     returns a list with these genotype frequencies (\code{q}) as well as the
#'     parental gamete frequencies (\code{p}).
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' freq <- c(0.5, 0, 0, 0, 0.5)
#' freqnext(freq = freq, alpha = 0)
#'
freqnext <- function(freq, alpha, segmat = NULL, more = FALSE) {
  ploidy <- length(freq) - 1
  stopifnot(ploidy %% 2 == 0, ploidy > 0)
  stopifnot(length(alpha) == floor(ploidy / 4))
  stopifnot(alpha >= 0, sum(alpha) <= 1)
  stopifnot(is.logical(more), length(more) == 1)

  if (is.null(segmat)) {
    segmat <- gsegmat(alpha = alpha, ploidy = ploidy)
  } else {
    stopifnot(dim(segmat) == c(ploidy + 1, ploidy / 2 + 1))
    stopifnot(abs(rowSums(segmat) - 1) < 10^-6)
  }

  p <- c(t(freq) %*% segmat)

  freqnew <- stats::convolve(p, rev(p), type = "open")

  ## resolve numerical issues
  freqnew[freqnew < 0] <- 0
  freqnew <- freqnew / sum(freqnew)

  if (more) {
    return(list(q = freqnew, p = p))
  } else {
    return(freqnew)
  }
}

#' Generate HWE genotype frequencies
#'
#' Generate genotype frequencies under Hardy-Weinberg equilibrium
#' given the allele frequency of the reference allele (\code{r}),
#' the double reduction parameter (\code{alpha}), and the ploidy
#' of the species (\code{ploidy}).
#'
#' If \code{alpha} is not all 0, then this function repeatedly
#' applies \code{\link{freqnext}()} to simulate genotype frequencies
#' under HWE. Otherwise, it uses \code{\link[stats]{dbinom}()}.
#'
#' @inheritParams dgamete
#' @param r The allele frequency of the reference allele.
#' @param niter The maximum number of iterations to simulate.
#' @param tol The stopping criterion on the Chi-square divergence between
#'     old and new genotype frequencies.
#' @param more A logical. Should we return more output (\code{TRUE}) or
#'     less (\code{FALSE}). See the Value section for details.
#'
#' @return If \code{more = FALSE}, then returns just the genotype frequencies
#'     after \code{niter} generations of random mating. If \code{more = TRUE},
#'     then returns a list with these genotype frequencies, as well as
#'     the parental gamete frequencies.
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' freq1 <- hwefreq(r = 0.5, alpha = 0, ploidy = 4)
#' freq2 <- hwefreq(r = 0.5, alpha = 1/6, ploidy = 4)
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
hwefreq <- function(r,
                    alpha,
                    ploidy,
                    niter = 100,
                    tol = sqrt(.Machine$double.eps),
                    more = FALSE) {
  stopifnot(length(r) == 1L, length(ploidy) == 1L, length(niter) == 1L)
  stopifnot(ploidy %% 2 == 0)
  stopifnot(ploidy > 1)
  stopifnot(length(alpha) == floor(ploidy / 4))
  stopifnot(alpha >= 0, sum(alpha) <= 1)
  stopifnot(r >= 0, r <= 1)
  stopifnot(niter >= 1)
  stopifnot(is.logical(more), length(more) == 1)

  ## Return theoretical result when no double reduction and large niter ----
  if (all(alpha < sqrt(.Machine$double.eps)) & niter >= 6) {
    freq <- stats::dbinom(x = 0:ploidy, size = ploidy, prob = r)
    if (more) {
      pgam <- stats::dbinom(x = 0:(ploidy / 2), size = ploidy / 2, prob = r)
      return(list(q = freq, p = pgam))
    } else {
      return(freq)
    }
  }

  ## Special code for tetraplolids at equilibrium and large niter ----
  if (ploidy == 4 & niter >= 6) {
    pgam <- gam_from_a(a = alpha, r = r)
    names(pgam) <- NULL
    freq <- stats::convolve(pgam, rev(pgam), type = "open")
    if (more) {
      return(list(q = freq, p = pgam))
    } else {
      return(freq)
    }
  }

  ## Create segregation matrix so don't need to remake it each iteration
  segmat <- gsegmat(alpha = alpha, ploidy = ploidy)

  ## Iterate freqnext() if double reduction ----
  freq <- c(1 - r, rep(0, length.out = ploidy - 1), r)
  i <- 1
  err <- Inf
  while (i <= niter && err > tol) {
    oldfreq <- freq
    freqlist <- freqnext(freq = freq, alpha = alpha, segmat = segmat, more = TRUE)
    freq <- freqlist$q
    pos <- freq > sqrt(.Machine$double.eps)
    i <- i + 1
    err <- sum((oldfreq[pos] - freq[pos]) ^ 2 / freq[pos]) * (ploidy + 1)
  }

  if (more) {
    return(freqlist)
  } else {
    return(freq)
  }
}


#' Obtain gamete frequencies at equilibrium given rates of double reduction.
#'
#' Given the rate of doulbe reduction and the major allele frequency, this
#' function will calculate the gametic frequencies.
#'
#' @inheritParams dgamete
#' @param p The allele frequency of the major allele.
#' @param ploidy The ploidy of the species.
#'
#' @return A numeric vector of length \code{ploidy / 2 + 1}, where element
#' \code{i} is the probability that a gamete carries \code{i-1} copies of
#' the major allele.
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' p_from_alpha(0.2, 0.5, 4)
#'
p_from_alpha <- function(alpha, p, ploidy) {
  stopifnot(length(ploidy) == 1, length(p) == 1)
  stopifnot(ploidy %% 2 == 0)
  stopifnot(length(alpha) == floor(ploidy / 4))
  stopifnot(alpha > 0)
  stopifnot(sum(alpha) <= 1)
  stopifnot(p >= 0, p <= 1)

  q <- hwefreq(r = p, alpha = alpha, ploidy = ploidy)

  pgamete <- rep(0, ploidy / 2 + 1)
  for (i in 0:ploidy) {
    pgamete <- pgamete + dgamete(x = 0:(ploidy/2),
                                 alpha = alpha,
                                 G = i,
                                 ploidy = ploidy,
                                 log_p = FALSE) * q[[i + 1]]
  }
  return(pgamete)
}

