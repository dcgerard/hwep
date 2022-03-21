###########################
## Utility functions
###########################

#' Greatest common divisor
#'
#' @param x An integer
#' @param y Anothe rinteger
#'
#' @author David Gerard
#'
#' @examples
#' gcd(105, 168)
#'
#' @noRd
gcd <- function(x, y) {
  r <- x %% y
  return(ifelse(r > 0, gcd(y, r), y))
}


#' Use log-sum-exp trick on vector
#'
#' Calculates \deqn{\log(\sum_i\exp(x_i))}
#'
#' @param x A numeric. We want to sum the exponent of x, then log.
#' @param na.rm A logical. Should we remove \code{NA}'s (\code{TRUE}) or not
#'     (\code{FALSE})? Defaults to \code{FALSE}.
#'
#' @author David Gerard
#'
#' @noRd
log_sum_exp <- function(x, na.rm = FALSE) {
  stopifnot(is.numeric(x))
  stopifnot(is.logical(na.rm))
  stopifnot(length(na.rm) == 1)

  x_max <- max(x, na.rm = na.rm)
  return(log(sum(exp(x - x_max), na.rm = na.rm)) + x_max)
}

logit <- function(x) {
  log(x) - log(1 - x)
}

expit <- function(x) {
  1 / (1 + exp(-x))
}


#' Function to aggregate groups.
#'
#' @param which_keep A logical vector. Should be \code{FALSE} if that
#'     group is being aggregated.
#'
#' @author David Gerard
#'
#' @noRd
aggfun <- function(which_keep) {
  stopifnot(is.logical(which_keep))
  ploidy <- length(which_keep) - 1
  numkeep <- sum(which_keep)
  if (numkeep >= ploidy) {
    H <- diag(ploidy + 1)
  } else {
    H <- matrix(0, nrow = numkeep + 1, ncol = ploidy + 1)
    j <- 1
    for (i in 0:ploidy) {
      if (which_keep[[i + 1]]) {
        H[j, i + 1] <- 1
        j <- j + 1
      } else {
        H[numkeep + 1, i + 1] <- 1
      }
    }
  }
  return(H)
}

#' Applies a transform from simplex to reals.
#'
#' Transform due to Betancourt (2012).
#'
#' @param q A point on the unit simplex.
#'
#' @references
#' \itemize{
#'   \item{Betancourt, M. (2012, May). Cruising the simplex:
#'         Hamiltonian Monte Carlo and the Dirichlet distribution.
#'         In AIP Conference Proceedings 31st (Vol. 1443, No. 1, pp. 157-164).
#'         American Institute of Physics. \doi{10.1063/1.3703631}}
#'   \item{\url{https://mc-stan.org/docs/2_18/reference-manual/simplex-transform-section.html}}
#' }
#'
#' @seealso \code{\link{real_to_simplex}()} for inverse.
#'
#' @author David Gerard
#'
#' @noRd
simplex_to_real <- function(q) {
  K <- length(q)
  stopifnot(K > 1)
  z <- q / (1 - c(0, cumsum(q)[-K]))
  y <- logit(z[-K]) - log(1 / (K - 1:(K-1)))
  return(y)
}

#' Applies a transform from reals to simplex
#'
#' Transform due to Betancourt (2012).
#'
#' @param y A point on the reals.
#'
#' @references
#' \itemize{
#'   \item{Betancourt, M. (2012, May). Cruising the simplex:
#'         Hamiltonian Monte Carlo and the Dirichlet distribution.
#'         In AIP Conference Proceedings 31st (Vol. 1443, No. 1, pp. 157-164).
#'         American Institute of Physics. \doi{10.1063/1.3703631}}
#'   \item{\url{https://mc-stan.org/docs/2_18/reference-manual/simplex-transform-section.html}}
#' }
#'
#' @seealso \code{\link{simplex_to_real}()} for inverse.
#'
#' @author David Gerard
#'
#' @noRd
real_to_simplex <- function(y) {
  stopifnot(length(y) > 0)
  K <- length(y) + 1
  z <- c(expit(y + log(1 / (K - 1:(K-1)))), 1)

  x <- rep(NA_real_, length.out = K)
  cumval <- 0
  for (i in seq_len(K)) {
    x[[i]] <- (1 - cumval) * z[[i]]
    cumval <- cumval + x[[i]]
  }

  return(x)
}



#' Choose which levels to aggregate
#'
#' @param x A non-negative numeric vector.
#' @param thresh A non-negative double. We aggregate when x < thresh
#'
#' @return A logical vector. It will place a \code{FALSE} wherever we will
#'     aggregate those cells. It aggregates small cells until all cells have
#'     at least the threshhold.
#'
#' @author David Gerard
#'
#' @noRd
choose_agg <- function(x, thresh = 0) {
  stopifnot(thresh >= 0, length(thresh) == 1)
  stopifnot(x >= 0)

  if (thresh == 0) {
    return(rep(TRUE, length.out = length(x)))
  }

  which_keep <- x >= thresh
  if (all(which_keep)) {
    return(rep(TRUE, length.out = length(x)))
  }

  aggvec <- c(x[which_keep], sum(x[!which_keep]))
  while(any(aggvec < thresh)) {
    which_keep[which_keep][which.min(x[which_keep])] <- FALSE
    aggvec <- c(x[which_keep], sum(x[!which_keep]))
  }

  return(which_keep)
}






