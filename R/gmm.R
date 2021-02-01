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

#' Calculate \code{\link{chisqdiv}()} with real parameterization.
#'
#' Uses \code{\link{real_to_simplex}()} to convert \code{y} to
#' \code{alpha}. Then calculates chi-square divergence with
#' \code{\link{chisqdiv}()}. Transformation is that of
#' Betancourt (2012).
#'
#' @inheritParams chisqdiv
#' @param y A real vector of length at least 1.
#'
#' @author David Gerard
#'
#' @references
#' \itemize{
#'   \item{Betancourt, M. (2012, May). Cruising the simplex:
#'         Hamiltonian Monte Carlo and the Dirichlet distribution.
#'         In AIP Conference Proceedings 31st (Vol. 1443, No. 1, pp. 157-164).
#'         American Institute of Physics.
#'         \href{https://doi.org/10.1063/1.3703631}{doi:10.1063/1.3703631}}
#'   \item{\url{https://mc-stan.org/docs/2_18/reference-manual/simplex-transform-section.html}}
#' }
#'
#' @noRd
obj_reals <- function(y, nvec, denom = c("expected", "observed")) {
  ploidy <- length(nvec) - 1
  stopifnot(floor(ploidy / 4) == length(y))
  denom <- match.arg(denom)
  if (length(y) == 0) {
    return(chisqdiv(nvec = nvec, alpha = NULL, denom = denom))
  } else {
    alpha_full <- real_to_simplex(y)
    return(chisqdiv(nvec = nvec, alpha = alpha_full[-1], denom = denom))
  }
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
#'   \item{\code{alpha}}{The estimated double reduction parameter(s).
#'       In diploids, this value is \code{NULL}.}
#'   \item{\code{chisq_hwe}}{The chi-square test statistic for testing
#'       against the null of equilibrium.}
#'   \item{\code{df_hwe}}{The degrees of freedom associated with
#'       \code{chisq_hwe}.}
#'   \item{\code{p_hwe}}{The p-value against the null of equilibrium.}
#'   \item{\code{chisq_alpha}}{The chi-square test statistic for
#'       testing against the null of no double reduction (conditional
#'       on equilibrium). In diploids, this value is \code{NULL}.}
#'   \item{\code{df_alpha}}{The degrees of freedom associated with
#'       \code{chisq_alpha}.}
#'   \item{\code{p_alpha}}{The p-value against the null of no double
#'       reduction (conditional on equilibrium). In diploids, this value
#'       is \code{NULL}.}
#' }
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' ## Diploid at various levels of deviation from HWE
#' hwemom(c(10, 20, 10))
#' hwemom(c(12, 16, 12))
#' hwemom(c(14, 12, 14))
#'
#' ## Tetraploid with strong deviation from HWE
#' nvec <- c(25, 0, 0, 0, 25)
#' hwemom(nvec)
#'
#' ## Hexaploid with exact frequencies at HWE
#' nvec <- round(hwefreq(p = 0.5,
#'                       alpha = 0.4,
#'                       ploidy = 6,
#'                       niter = 100,
#'                       tol = -Inf) * 100)
#' hwemom(nvec)
#'
#' ## Octoploid case with exact frequencies at HWE
#' nvec <- round(hwefreq(p = 0.5,
#'                       alpha = c(0.4, 0.1),
#'                       ploidy = 8,
#'                       niter = 100,
#'                       tol = -Inf) * 100)
#' hwemom(nvec)
#'
hwemom <- function(nvec, denom = c("expected", "observed")) {
  ## Check parameters ----
  ploidy <- length(nvec) - 1
  stopifnot(ploidy %% 2 == 0)
  stopifnot(nvec >= 0)
  ibdr <- floor(ploidy / 4)
  denom <- match.arg(denom)

  if (ibdr == 0) {
    ## Diploid: Just return chi-square
    chisq <- chisqdiv(nvec = nvec, alpha = NULL, denom = denom)
    pval <- stats::pchisq(q = chisq,
                          df = ploidy - ibdr - 1,
                          lower.tail = FALSE)
    retlist <- list(alpha = NULL,
                    chisq_hwe = chisq,
                    df_hwe = ploidy - ibdr - 1,
                    p_hwe = pval,
                    chisq_alpha = NULL,
                    df_alpha = NULL,
                    p_alpha = NULL)
  } else if (ibdr == 1) {
    ## Tetraploid or Hexaploid: Use Brent's method
    oout <- stats::optim(par = 0.1,
                         fn = chisqdiv,
                         method = "Brent",
                         lower = 0,
                         upper = 1,
                         nvec = nvec,
                         denom = denom)
    pval_hwe <- stats::pchisq(q = oout$value,
                              df = ploidy - ibdr - 1,
                              lower.tail = FALSE)
    chisq_null <- chisqdiv(nvec = nvec,
                           alpha = rep(0, length.out = ibdr),
                           denom = denom)
    chisq_alpha <- 2 * (chisq_null - oout$value)
    pval_alpha <- stats::pchisq(q = chisq_alpha,
                                df = ibdr,
                                lower.tail = FALSE) / 2
    retlist <- list(alpha = oout$par,
                    chisq_hwe = oout$value,
                    df_hwe = ploidy - ibdr - 1,
                    p_hwe = pval_hwe,
                    chisq_alpha = chisq_alpha,
                    df_alpha = ibdr,
                    p_alpha = pval_alpha)
  } else {
    ## Higher Ploidy: Use BFGS on transformed space
    oout <- stats::optim(par = rep(0, length.out = ibdr),
                         fn = obj_reals,
                         method = "BFGS",
                         nvec = nvec,
                         denom = denom)
    alpha <- real_to_simplex(oout$par)[-1]
    pval_hwe <- stats::pchisq(q = oout$value,
                              df = ploidy - ibdr - 1,
                              lower.tail = FALSE)
    chisq_null <- chisqdiv(nvec = nvec,
                           alpha = rep(0, length.out = ibdr),
                           denom = denom)
    chisq_alpha <- 2 * (chisq_null - oout$value)
    pval_alpha <- stats::pchisq(q = chisq_alpha,
                                df = ibdr,
                                lower.tail = FALSE) / 2
    retlist <- list(alpha = alpha,
                    chisq_hwe = oout$value,
                    df_hwe = ploidy - ibdr - 1,
                    p_hwe = pval_hwe,
                    chisq_alpha = chisq_alpha,
                    df_alpha = ibdr,
                    p_alpha = pval_alpha)
  }

  return(retlist)
}
