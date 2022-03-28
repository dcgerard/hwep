########################
## Methods to estimate double reduction in F1 populations
## with known parental dosages. Mostly as a comparison to
## HWE approach.
########################

#' Estimate Double Reduction in F1 Populations
#'
#' Estimates double reduction in F1 populations by maximum likelihood.
#'
#' @param nvec A vector containing the observed genotype counts,
#'     where \code{nvec[[i]]} is the number of individuals with genotype
#'     \code{i-1}. This should be of length \code{ploidy+1}.
#' @inheritParams zygdist
#'
#' @return A list with some or all of the following elements:
#' \describe{
#'   \item{\code{alpha}}{A vector of numerics of length
#'     \code{floor(ploidy / 4)}, the estimated double reduction rate.}
#'   \item{\code{llike}}{The final log-likelihood.}
#' }
#'
#' @author David Gerard
#'
#' @seealso \code{\link{zygdist}()} for calculating the probability of
#'    offpring genotypes given parental genotypes and the double reduction
#'    rate.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' size <- 100
#' qvec <- zygdist(alpha = 0.1, G1 = 2, G2 = 2, ploidy = 4)
#' nvec <- c(stats::rmultinom(n = 1, size = size, prob = qvec))
#' f1dr(nvec = nvec, G1 = 2, G2 = 2)
#'
f1dr <- function(nvec, G1, G2) {
  ploidy <- length(nvec) - 1
  stopifnot(ploidy > 2, ploidy %% 2 == 0)
  stopifnot(nvec >= 0)
  stopifnot(length(G1) == 1, length(G2) == 1,
            G1 >= 0, G1 <= ploidy,
            G2 >= 0, G2 <= ploidy)

  upper <- drbounds(ploidy = ploidy)
  ibdr <- length(upper)
  lower <- rep(sqrt(.Machine$double.eps), length.out = ibdr)
  alpha_init <- rep(0, length.out = ibdr)
  method <- ifelse(ibdr == 1, "Brent", "L-BFGS-B")

  if ((G1 == 0 & G2 == 0) | (G1 == 4 & G2 == 4)) {
    retlist <- list()
    retlist$alpha <- rep(NA_real_, length.out = ibdr)
    retlist$llike <- NA_real_
    return(retlist)
  }

  oout <- stats::optim(par = alpha_init,
                       fn = f1obj,
                       method = method,
                       lower = lower,
                       upper = upper,
                       control = list(fnscale = -1),
                       nvec = nvec,
                       G1 = G1,
                       G2 = G2)

  retlist <- list()
  retlist$alpha <- oout$par
  retlist$llike <- oout$value
  return(retlist)
}

#' Log-likelihood of nvec given double reduction in F1 populations
#'
#' Used in \code{f1dr()}.
#'
#' @inheritParams f1dr
#' @param alpha The current double reduction rate.
#'
#' @author David Gerard
#'
#' @noRd
f1obj <- function(alpha, nvec, G1, G2) {
  qvec <- zygdist(alpha = alpha, G1 = G1, G2 = G2, ploidy = length(nvec) - 1)
  return(stats::dmultinom(x = nvec, prob = qvec, log = TRUE))
}
