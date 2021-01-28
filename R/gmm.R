###############################
## Functions for Generalized Method of Moments
###############################

#' Chi-square divergence between empirical proportions and updated proportions.
#'
#' @param nvec A vector of length ploidy + 1. \code{nvec[i]} is the observed
#'     number of individuals with dosage \code{i-1}.
#' @param alpha The double reduction parameter. Should be of length
#'     \code{floor(ploidy/4)}. These should sum to less than 1, since
#'     \code{1-sum(alpha)} is the assumed probability of no double
#'     reduction.
#' @param denom What should we use for the denominator? The expected
#'     genotype frequencies (\code{denom = "expected"}) or the
#'     observed genotype frequencies (\code{denom = "observed"})?
#'
#' @author David Gerard
#'
#' @noRd
chisqdiv <- function(nvec, alpha, denom = c("expected", "observed")) {
  ## Check parameters ----
  ploidy <- length(nvec) - 1
  stopifnot(ploidy %% 2 == 0)
  stopifnot(length(alpha) == floor(ploidy / 4))
  stopifnot(alpha >= 0, sum(alpha) <= 1)
  stopifnot(nvec >= 0)
  denom <- match.arg(denom)

  ## calculate divergence ----
  n <- sum(nvec)
  qhat <- nvec / n
  fq <- freqnext(freq = qhat, alpha = alpha)

  if (denom == "expected") {
    chisq <- n * sum((qhat - fq) ^ 2 / fq)
  } else if (denom == "observed") {
    chisq <- n * sum((qhat - fq) ^ 2 / qhat)
  }

  return(chisq)
}

#' Estimate double reduction parameter and test for HWE
#'
#' Implements a generalized method-of-moments approach to test for
#' Hardy-Weinberg equilibrium in the presence of double reduction. This
#' also provides an estimate of double reduction given that the population
#' is at equilibrium.
#'
#' @param nvec A vector containing the observed genotype counts,
#'     where \code{nvec[[i]]} is the number of individuals with genotype
#'     \code{i-1}. This should be of length \code{ploidy+1}.
#' @param denom What should we use for the denominator? The expected
#'     genotype frequencies (\code{denom = "expected"}) or the
#'     observed genotype frequencies (\code{denom = "observed"})?
#'
#' @return A list with the following elements
#' \describe{
#'   \item{\code{alpha}}{The estimated double reduction parameter(s).}
#'   \item{\code{chisq}}{The chi-square test statistic for testing
#'       against the null of equilibrium.}
#'   \item{\code{p.value}}{The p-value against the null of equilibrium.}
#' }
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' ## Diploid at various levels of deviation from HWE
#' hwefit(c(10, 20, 10))
#' hwefit(c(12, 16, 12))
#' hwefit(c(14, 12, 14))
#'
#' ## Tetraploid with strong deviation from HWE
#' nvec <- c(25, 0, 0, 0, 25)
#' hwefit(nvec)
#'
#' ## Hexaploid with exact frequencies at HWE
#' nvec <- round(hwefreq(p = 0.5, alpha = 0.4, ploidy = 6, niter = 100, tol = -Inf) * 50)
#' hwefit(nvec)
#'
hwefit <- function(nvec, denom = c("expected", "observed")) {
  ## Check parameters ----
  ploidy <- length(nvec) - 1
  stopifnot(ploidy %% 2 == 0)
  stopifnot(nvec >= 0)
  ibdr <- floor(ploidy / 4)
  denom <- match.arg(denom)

  if (ibdr == 0) {
    ## Diploid: Just return chi-square
    chisq <- chisqdiv(nvec = nvec, alpha = numeric(length = 0), denom = denom)
    pval <- stats::pchisq(q = chisq,
                          df = ploidy - ibdr - 1,
                          lower.tail = FALSE)
    retlist <- list(alpha = numeric(length = 0),
                    chisq = chisq,
                    p.value = pval)
  } else if (ibdr == 1) {
    ## Tetraploid or Hexaploid: Use Brent's method
    oout <- stats::optim(par = 0.1,
                         fn = chisqdiv,
                         method = "Brent",
                         lower = 0,
                         upper = 1,
                         nvec = nvec,
                         denom = denom)
    pval <- stats::pchisq(q = oout$value,
                          df = ploidy - ibdr - 1,
                          lower.tail = FALSE)
    retlist <- list(alpha = oout$par,
                    chisq = oout$value,
                    p.value = pval)
  } else {
    ## Higher Ploidy: Use L-BFGS-B on transformed space
  }

  return(retlist)
}
