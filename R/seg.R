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
  ## stopifnot(alpha >= 0, sum(alpha) <= 1)

  ## Faster if ploidy is 2 or if no DR to use hypergeometric dist ----
  if (ploidy == 2) {
    retvec <- stats::dhyper(x = x,
                            m = G,
                            n = ploidy - G,
                            k = ploidy / 2,
                            log = log_p)
    return(retvec)
  }
  if (all(abs(alpha) < sqrt(.Machine$double.eps))) {
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
    retvec[[k]] <-
      sum(
        exp(
          lchoose(G, jvec) +
            lchoose(G - jvec, x[[k]] - 2 * jvec) +
            lchoose(ploidy - G, ivec - jvec) +
            lchoose(ploidy - G - (ivec - jvec), ploidy / 2 - x[[k]] - 2 * (ivec - jvec)) -
            lchoose(ploidy, ivec) -
            lchoose(ploidy - ivec, ploidy / 2 - 2 * ivec)
          ) * alphavec
        )

    ## old way
    # exp(
    #   log_sum_exp(
    #     lchoose(G, jvec) +
    #       lchoose(G - jvec, x[[k]] - 2 * jvec) +
    #       lchoose(ploidy - G, ivec - jvec) +
    #       lchoose(ploidy - G - (ivec - jvec), ploidy / 2 - x[[k]] - 2 * (ivec - jvec)) -
    #       lchoose(ploidy, ivec) -
    #       lchoose(ploidy - ivec, ploidy / 2 - 2 * ivec) +
    #       log(alphavec)
    #     )
    #   )
  }
  retvec[is.nan(retvec)] <- -Inf ## just means multiple -Inf in the above

  if (log_p) {
    retvec <- log(retvec)
  }

  return(retvec)
}

#' Less computationally efficient version of \code{\link{gsegmat}()}
#'
#' Produces the segregation probabilities for gamete dosages given
#' parental dosages and the double reduction rate.
#'
#' @inheritParams dgamete
#'
#' @return A matrix of dimension \code{ploidy + 1} by \code{ploidy / 2 + 1}.
#'     Element (i, j) is the probability that a parent carrying dosage
#'     j - 1 produces a gamete with dosage i - 1.
#'
#' @author David Gerard
#'
#' @noRd
#'
#' @examples
#' gsegmat(alpha = 1/6, ploidy = 4)
#'
gsegmat2 <- function(alpha, ploidy) {
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


#' Segregation probabilities of gametes
#'
#' Produces the segregation probabilities for gamete dosages given
#' parental dosages and the double reduction rate.
#'
#' @inheritParams dgamete
#'
#' @return A matrix of dimension \code{ploidy + 1} by \code{ploidy / 2 + 1}.
#'     Element (i, j) is the probability that a parent carrying dosage
#'     j - 1 produces a gamete with dosage i - 1.
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' gsegmat(alpha = NULL, ploidy = 2)
#'
#' gsegmat(alpha = 1/6, ploidy = 4)
#'
#' gsegmat(alpha = 0.3, ploidy = 6)
#'
#' gsegmat(alpha = c(0.35, 0.02), ploidy = 8)
#'
#' gsegmat(alpha = c(0.4, 0.05), ploidy = 10)
#'
gsegmat <- function(alpha, ploidy) {

  ## Check lengths ----
  stopifnot(ploidy %% 2 == 0, ploidy >= 2, length(ploidy) == 1)
  stopifnot(length(alpha) == floor(ploidy / 4))

  ## Early return for special cases ----
  if (ploidy == 2) {
    return(gsegmat_diploid())
  } else if (ploidy == 4) {
    return(gsegmat_tetraploid(alpha = alpha))
  } else if (ploidy == 6) {
    return(gsegmat_hexaploid(alpha = alpha))
  }

  ## Keep checking ----
  ## stopifnot(alpha >= 0, sum(alpha) <= 1)
  ibdr <- floor(ploidy / 4)
  alphavec <- c(1 - sum(alpha), alpha)
  ellvec <- 0:ploidy
  kvec <- 0:(ploidy / 2)

  ## Calculate segregation matrix for higher ploidies ----
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
            matrix(lchoose(ploidy - i, ploidy / 2 - 2 * i), nrow = ploidy + 1, ncol = ploidy / 2 + 1)
        ) * alphavec[[i + 1]]
    }
  }

  colnames(segmat) <- 0:(ploidy / 2)
  rownames(segmat) <- 0:ploidy
  return(segmat)
}


#' Symbolic representation of the segregation probability matrix
#'
#' Two alleles are identical-by-double-reduction (IBDR) if they originate from
#' the same (by origin) allele in the parent. We let "a" be the probability of
#' zero IBDR alleles, "b" be the probability of one IBDR pair,
#' "c" be the probability of two IBDR pairs, etc...
#'
#' @param ploidy The ploidy of the species
#' @param out Should we return a character matrix
#'     (\code{"str"}) or an expression matrix (\code{"exp"})?
#'
#' @return A character or expression matrix containing the mathematical
#'     form for the segregation matrix. Element (i, j) is the probability
#'     a parent with dosage i-1 produces a gamete with dosage j-1.
#'
#' @seealso \code{\link{gsegmat}()} for numerical expressions.
#'
#' @export
#'
#' @author David Gerard
#'
#' @examples
#' gsegmat_symb(4)
#' gsegmat_symb(6)
#' gsegmat_symb(8)
gsegmat_symb <- function(ploidy, out = c("str", "exp")) {
  stopifnot(ploidy %% 2 == 0, ploidy >= 2, length(ploidy) == 1)
  out <- match.arg(out)
  ibdr <- floor(ploidy / 4)

  ## Early stopping for diploids
  if (ploidy == 2) {
    segmat <- matrix(c("1", "0",
                       "1/2", "1/2",
                       "0", "1"),
                     ncol = 2,
                     byrow = TRUE)
    if (out == "exp") {
      segmat <- matrix(str2expression(segmat), nrow = ploidy + 1)
    }
    return(segmat)
  }

  avec <- letters[0:ibdr + 1]
  ellvec <- 0:ploidy
  kvec <- 0:(ploidy / 2)

  segmat <- matrix("", nrow = ploidy + 1, ncol = ploidy / 2 + 1)
  for (i in 0:ibdr) {
    for (j in 0:i) {
      nummat_num <- outer(X = choose(ellvec, j), Y = rep(1, ploidy / 2 + 1), FUN = `*`) *
        outer(X = ellvec, Y = kvec, FUN = function(x, y) choose(x - j, y - 2 * j)) *
        outer(X = choose(n = ploidy - ellvec, i - j), Y = rep(1, ploidy / 2 + 1), FUN = `*`) *
        outer(X = ellvec, Y = kvec, FUN = function(x, y) choose(ploidy - x - (i - j), ploidy / 2 - y - 2 * (i - j)))

      denmat_num <- matrix(choose(ploidy, i), nrow = ploidy + 1, ncol = ploidy / 2 + 1) *
        matrix(choose(ploidy - i, ploidy / 2 - 2 * i), nrow = ploidy + 1, ncol = ploidy / 2 + 1)

      gcdmat <- matrix(mapply(nummat_num, denmat_num, FUN = gcd), nrow = ploidy + 1)

      nummat_num <- nummat_num / gcdmat
      denmat_num <- denmat_num / gcdmat

      nummat <- nummat_num
      class(nummat) <- "character"
      denmat <- denmat_num
      class(denmat) <- "character"

      addmat <- matrix(paste0("(", nummat, "/", denmat, ")"), nrow = ploidy + 1)

      addmat[nummat_num == denmat_num] <- ""
      addmat <- matrix(paste0("+", addmat, "*", avec[[i + 1]]), nrow = ploidy + 1)
      addmat[nummat_num == 0] <- ""

      segmat <- matrix(paste0(segmat, addmat), nrow = ploidy + 1)
    }
  }

  segmat[segmat == ""] <- "0"
  segmat <- gsub(pattern = "^\\+", replacement = "", x = segmat)
  segmat <- gsub(pattern = "^\\*", replacement = "", x = segmat)
  segmat <- gsub(pattern = "\\+\\*", replacement = "\\+", x = segmat)
  if (out == "str") {
    segmat <- gsub(pattern = "\\*", replacement = "", x = segmat)
  } else {
    segmat <- matrix(str2expression(segmat), nrow = ploidy + 1)
  }

  colnames(segmat) <- 0:(ploidy / 2)
  rownames(segmat) <- 0:ploidy

  return(segmat)
}


#' Special case of gsegmat() for diploids.
#'
#' @author David Gerard
#'
#' @noRd
gsegmat_diploid <- function() {
  segmat <- matrix(c(1.0, 0.0,
                     0.5, 0.5,
                     0.0, 1.0),
                   byrow = TRUE,
                   ncol = 2,
                   nrow = 3)
  colnames(segmat) <- 0:1
  rownames(segmat) <- 0:2
  return(segmat)
}

#' Special case of gsegmat() for tetraploids
#'
#' @author David Gerard
#'
#' @noRd
gsegmat_tetraploid <- function(alpha) {
  stopifnot(length(alpha) == 1)
  ## stopifnot(alpha >= 0, alpha <= 1)
  segmat <- matrix(c(1, 0, 0,
                     0.25 * (2 + alpha), 0.5 * (1 - alpha), 0.25 * alpha,
                     (1 + 2 * alpha) / 6, 4 * (1 - alpha) / 6, (1 + 2 * alpha) / 6,
                     0.25 * alpha, 0.5 * (1 - alpha), 0.25 * (2 + alpha),
                     0, 0, 1),
                   byrow = TRUE,
                   ncol = 3,
                   nrow = 5)
  colnames(segmat) <- 0:2
  rownames(segmat) <- 0:4
  return(segmat)
}

#' Special case of gsegmat() for tetraploids
#'
#' @author David Gerard
#'
#' @noRd
gsegmat_hexaploid <- function(alpha) {
  stopifnot(length(alpha) == 1)
  ## stopifnot(alpha >= 0, alpha <= 1)
  segmat <- matrix(c(1, 0, 0, 0,
                     (3 + alpha) / 6, (3 - 2 * alpha) / 6, alpha / 6, 0,
                     (3 + 3 * alpha) / 15, (9 - 5 * alpha) / 15, (3 + alpha) / 15, alpha / 15,
                     (1 + 3 * alpha) / 20, (9 - 3 * alpha) / 20, (9 - 3 * alpha) / 20, (1 + 3 * alpha) / 20,
                     alpha / 15, (3 + alpha) / 15, (9 - 5 * alpha) / 15, (3 + 3 * alpha) / 15,
                     0, alpha / 6, (3 - 2 * alpha) / 6, (3 + alpha) / 6,
                     0, 0, 0, 1),
                   byrow = TRUE,
                   ncol = 4,
                   nrow = 7)
  colnames(segmat) <- 0:3
  rownames(segmat) <- 0:6
  return(segmat)
}

#' Zygote dosage probabilities.
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
  ## stopifnot(alpha >= 0, sum(alpha) <= 1)

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
  ## stopifnot(alpha >= 0, sum(alpha) <= 1)

  segarray <- array(data = NA_real_, dim = rep(ploidy + 1, length.out = 3))
  for (G1 in 0:ploidy) {
    for (G2 in 0:ploidy) {
      segarray[G1 + 1, G2 + 1, ] <- zygdist(alpha = alpha,
                                            G1 = G1,
                                            G2 = G2,
                                            ploidy = ploidy)
    }
  }

  dimnames(segarray) <- list(offspring = 0:ploidy,
                             parent1 = 0:ploidy,
                             parent2 = 0:ploidy)

  return(segarray)
}
