###################################
## Functions to calculate segregation probabilities
###################################

#' Gamete dosage probability
#'
#' Estimates the probability of a gamete dosage given the parent dosage
#' (\code{G}), the parent ploidy (\code{ploidy}), and the double reduction
#' parameter (\code{alpha}). This is for biallelic loci.
#'
#' @param x A vector of numerics in \code{seq(0, ploidy/2)}. The dosage of the
#'     gametes.
#' @param alpha A numeric vector containing the double reduction parameter(s).
#'     This should be a
#'     vector of length \code{floor(ploidy/4)} where \code{alpha[i]}
#'     is the probability of exactly \code{i} pairs of IBDR alleles
#'     being in the gamete. Note that \code{sum(alpha)} should be less than
#'     1, as \code{1 - sum(alpha)} is the probability of no double reduction.
#' @param G The dosage of the parent. Should be an integer between \code{0}
#'     and \code{ploidy}.
#' @param ploidy The ploidy of the species. This should be an even positive
#'     integer.
#' @param log_p A logical. Should we return the log-probability (\code{TRUE})
#'     or not (\code{FALSE})? Defaults to \code{FALSE}.
#'
#' @return A vector of length \code{length(x)}, containing the (log)
#'     probabilities of a gamete carrying a dosage of \code{x} from a
#'     parent of dosage \code{G} who has ploidy \code{ploidy} and a
#'     double reduction rate \code{alpha}.
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' dgamete(x = 0:2, alpha = 0, G = 2, ploidy = 4)
#'
dgamete <- function(x, alpha, G, ploidy, log_p = FALSE) {
  ## Check input ----
  stopifnot(length(ploidy) == 1L, length(G) == 1L, length(log_p) == 1L)
  stopifnot(is.logical(log_p))
  stopifnot(ploidy %% 2 == 0)
  stopifnot(0 <= G, G <= ploidy)
  stopifnot(0 <= x, x <= ploidy / 2)
  ibdr <- floor(ploidy / 4)
  stopifnot(length(alpha) == ibdr)
  stopifnot(alpha >= 0, sum(alpha) <= 1)


  ## Faster if ploidy is 2 or if no DR to use hypergeometric dist ----
  if (ploidy == 2) {
    retvec <- stats::dhyper(x = x,
                            m = G,
                            n = ploidy - G,
                            k = ploidy / 2,
                            log = log_p)
    return(retvec)
  }
  if (all(alpha < sqrt(.Machine$double.eps))) {
    retvec <- stats::dhyper(x = x,
                            m = G,
                            n = ploidy - G,
                            k = ploidy / 2,
                            log = log_p)
    return(retvec)
  }

  ## Get sequence of indices to sum over ----
  ijmat <- cbind(utils::combn(x = 0:ibdr, m = 2),
                 matrix(rep(0:ibdr, each = 2), nrow = 2))
  jvec <- ijmat[1, ]
  ivec <- ijmat[2, ]
  alpha <- c(1 - sum(alpha), alpha)
  alphavec <- alpha[ivec + 1]

  ## Calculate probs ----
  retvec <- rep(NA_real_, length = length(x))
  for (k in seq_along(x)) {
    retvec[[k]] <- log_sum_exp(
      lchoose(G, jvec) +
        lchoose(G - jvec, x[[k]] - 2 * jvec) +
        lchoose(ploidy - G, ivec - jvec) +
        lchoose(ploidy - G - (ivec - jvec), ploidy / 2 - x[[k]] - 2 * (ivec - jvec)) -
        lchoose(ploidy, ivec) -
        lchoose(ploidy - ivec, ploidy / 2 - 2 * ivec) +
        log(alphavec)
    )
  }
  retvec[is.nan(retvec)] <- -Inf ## just means multiple -Inf in the above

  if (!log_p) {
    retvec <- exp(retvec)
  }

  return(retvec)
}

#' Segregation probabilities of gametes
#'
#' Produces the segregation probabilities for gamete dosages given
#' parental dosages and the double reduction rate.
#'
#' @inheritParams dgamete
#'
#' @return A matrix of dimension \code{ploidy + 1} by \code{ploidy / 2 + 1}.
#'     Element (i, j) is the probability that a parent carying dosage
#'     j - 1 produces a gamete with dosage i - 1.
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' gsegmat(alpha = 1/6, ploidy = 4)
#'
gsegmat <- function(alpha, ploidy) {
  stopifnot(ploidy %% 2 == 0, ploidy > 0)
  stopifnot(length(alpha) == floor(ploidy / 4))
  stopifnot(alpha >= 0, sum(alpha) <= 1)
  segmat <- matrix(NA_real_, nrow = ploidy + 1, ncol = ploidy / 2 + 1)
  for (i in 0:ploidy) {
    segmat[i + 1, ] <- dgamete(x = 0:(ploidy / 2),
                               alpha = alpha,
                               G = i,
                               ploidy = ploidy)
  }
  return(segmat)
}


#' Alternative way to do \code{\link{gsegmat}()}
#'
#' @inheritParams gsegmat
#'
#' @author David Gerard
#'
#' @noRd
gsegmat2 <- function(alpha, ploidy) {
  stopifnot(ploidy %% 2 == 0, ploidy >=4, length(ploidy) == 1)
  stopifnot(length(alpha) == floor(ploidy / 4))
  stopifnot(alpha >= 0, sum(alpha) <= 1)
  ibdr <- floor(ploidy / 4)
  alphavec <- c(1 - sum(alpha), alpha)
  ellvec <- 0:ploidy
  kvec <- 0:(ploidy / 2)

  segmat <- matrix(0, nrow = ploidy + 1, ncol = ploidy / 2 + 1)
  for (i in 0:ibdr) {
    for (j in 0:i) {
      segmat <- segmat +
        exp(
          outer(X = lchoose(ellvec, j), Y = rep(1, ploidy / 2 + 1), FUN = `*`) +
            outer(X = ellvec, Y = kvec, FUN = function(x, y) lchoose(x - j, y - 2 * j)) +
            outer(X = lchoose(n = ploidy - ellvec, i - j), Y = rep(1, ploidy / 2 + 1), FUN = `*`) +
            outer(X = ellvec, Y = kvec, FUN = function(x, y) lchoose(ploidy - x - (i - j), ploidy / 2 - y - 2 * (i - j))) -
            matrix(lchoose(ploidy, i), nrow = ploidy + 1, ncol = ploidy / 2 + 1) -
            matrix(lchoose(ploidy - i, ploidy / 2 - 2 * i), nrow = ploidy + 1, ncol = ploidy / 2 + 1) +
            log(alphavec[[i + 1]])
        )
    }
  }
  return(segmat)
}

#' Zygote dosage probabiltites.
#'
#' Calculates the distribution of an offspring dosages given
#' parental dosages (\code{G1} and \code{G2}), the ploidy of the
#' species (\code{ploidy}), and the double reduction parameter
#' (\code{alpha}).
#'
#' @inheritParams dgamete
#' @param G1 The dosage of parent 1. Should be an integer between \code{0}
#'     and \code{ploidy}.
#' @param G2 The dosage of parent 2. Should be an integer between \code{0}
#'     and \code{ploidy}.
#'
#' @return A vector of probabilities. The \code{i}th element is the
#'     probability that the offspring will have dosage \code{i-1}.
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' zygdist(alpha = c(0.5, 0.1), G1 = 4, G2 = 5, ploidy = 8)
#'
zygdist <- function(alpha, G1, G2, ploidy) {
  ## Check parameters ----
  stopifnot(length(G1) == 1L,
            length(G2) == 1L,
            length(ploidy) == 1L)
  stopifnot(ploidy %% 2 == 0)
  stopifnot(0 <= G1, G1 <= ploidy)
  stopifnot(0 <= G2, G2 <= ploidy)
  ibdr <- floor(ploidy / 4)
  stopifnot(length(alpha) == ibdr)
  stopifnot(alpha >= 0, sum(alpha) <= 1)

  ## Get gamete probs ----
  p1gamprob <- dgamete(x = 0:(ploidy / 2),
                       alpha = alpha,
                       G = G1,
                       ploidy = ploidy,
                       log_p = FALSE)
  p2gamprob <- dgamete(x = 0:(ploidy / 2),
                       alpha = alpha,
                       G = G2,
                       ploidy = ploidy,
                       log_p = FALSE)

  ## Convolve to get zygote probabilities ----
  zygdist <- stats::convolve(p1gamprob, rev(p2gamprob), type = "open")

  return(zygdist)
}

#' Zygote segregation distributions.
#'
#' Obtains offspring genotype probabilities given parental probabilities,
#' the ploidy of the species, and the overdispersion parameter,
#' for all possible parental genotypes.
#'
#' @inheritParams dgamete
#'
#' @return An array of probabilities. Element (i, j, k) contains the
#'     probability of offspring dosage k-1 given parental dosages
#'     i-1 and j-1.
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' ploidy <- 10
#' alpha <- c(0.5, 0.1)
#' p1 <- 4
#' p2 <- 3
#' segarray <- zsegarray(alpha = alpha, ploidy = ploidy)
#' graphics::plot(x = 0:10,
#'                y = segarray[p1 + 1, p2 + 1, ],
#'                type = "h",
#'                ylab = "Pr(dosage)",
#'                xlab = "dosage")
#' graphics::mtext(paste0("P1 dosage = ",
#'                        p1,
#'                        ", ",
#'                        "P2 dosage = ",
#'                        p2))
#'
zsegarray <- function(alpha, ploidy) {
  stopifnot(length(ploidy) == 1L)
  stopifnot(ploidy %% 2 == 0)
  stopifnot(length(alpha) == floor(ploidy / 4))
  stopifnot(alpha >= 0, sum(alpha) <= 1)

  segarray <- array(data = NA_real_, dim = rep(ploidy + 1, length.out = 3))
  for (G1 in 0:ploidy) {
    for (G2 in 0:ploidy) {
      segarray[G1 + 1, G2 + 1, ] <- zygdist(alpha = alpha,
                                            G1 = G1,
                                            G2 = G2,
                                            ploidy = ploidy)
    }
  }

  return(segarray)
}
