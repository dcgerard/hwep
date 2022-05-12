#################
## Bayes tests for random mating
#################

#' PMF of Dirichlet-multinomial distribution
#'
#' @param x The vector of counts.
#' @param alpha The vector of concentration parameters.
#' @param lg A logical. Should we log the density (TRUE) or not (FALSE)
#'
#' @author David Gerard
#'
#' @noRd
ddirmult <- function(x, alpha, lg = FALSE) {
  stopifnot(length(x) == length(alpha))
  stopifnot(all(alpha > 0))
  asum <- sum(alpha)
  n <- sum(x)
  ll <- lgamma(asum) + lgamma(n + 1) - lgamma(n + asum) +
    sum(lgamma(x + alpha)) - sum(lgamma(alpha)) - sum(lgamma(x + 1))
  if (!lg) {
    ll <- exp(ll)
  }
  return(ll)
}

#' Bayes test for random mating with known genotypes
#'
#' @param nvec A vector containing the observed genotype counts,
#'     where \code{nvec[[i]]} is the number of individuals with genotype
#'     \code{i-1}. This should be of length \code{ploidy+1}.
#' @param lg A logical. Should we return the log Bayes factor (\code{TRUE})
#'     or the Bayes factor (\code{FALSE})?
#'
#' @examples
#' set.seed(1)
#' ploidy <- 8
#'
#' ## Simulate under the null
#' p <- stats::runif(ploidy / 2 + 1)
#' p <- p / sum(p)
#' q <- stats::convolve(p, rev(p), type = "open")
#'
#' ## See BF increase
#' nvec <- c(stats::rmultinom(n = 1, size = 100, prob = q))
#' hwebrm(nvec = nvec)
#'
#' nvec <- c(stats::rmultinom(n = 1, size = 1000, prob = q))
#' hwebrm(nvec = nvec)
#'
#' nvec <- c(stats::rmultinom(n = 1, size = 10000, prob = q))
#' hwebrm(nvec = nvec)
#'
#' ## Simulate under the alternative
#' q <- stats::runif(ploidy + 1)
#' q <- q / sum(q)
#'
#' ## See BF decrease
#' nvec <- c(stats::rmultinom(n = 1, size = 100, prob = q))
#' hwebrm(nvec = nvec)
#'
#' nvec <- c(stats::rmultinom(n = 1, size = 1000, prob = q))
#' hwebrm(nvec = nvec)
#'
#' nvec <- c(stats::rmultinom(n = 1, size = 10000, prob = q))
#' hwebrm(nvec = nvec)
#'
#' @author David Gerard
#'
#' @export
hwebrm <- function(nvec, lg = TRUE) {
  ploidy <- length(nvec) - 1
  alpha <- rep(1, ploidy / 2 + 1)
  alpha <- alpha / sum(alpha)
  beta <- stats::convolve(alpha, rev(alpha), type = "open")
  beta <- beta / sum(beta)

  mnull <- gibbs_known(x = nvec, alpha = alpha, more = FALSE, lg = TRUE)$mx
  malt <- ddirmult(x = nvec, alpha = beta, lg = TRUE)

  lbf = mnull - malt

  if (!lg) {
    lbf <- exp(lbf)
  }

  return(lbf)
}
