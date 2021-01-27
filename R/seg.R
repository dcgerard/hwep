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
#' @param alpha The double reduction parameter(s). This should be a
#'     vector of length \code{floor(ploidy/4)} where \code{alpha[i]}
#'     is the probability of exactly \code{i} pairs of IBDR alleles
#'     being in the gamete. Note that \code{sum(alpha)} should be less than
#'     1, as \code{1 - sum(alpha)} is the probability of no double reduction.
#' @param G The dosage of the parent. Should be one of \code{seq(0, ploidy)}.
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

  ijmat <- cbind(utils::combn(x = 0:ibdr, m = 2),
                 matrix(rep(0:ibdr, each = 2), nrow = 2))
  jvec <- ijmat[1, ]
  ivec <- ijmat[2, ]
  alpha <- c(1 - sum(alpha), alpha)
  alphavec <- alpha[ivec + 1]

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
