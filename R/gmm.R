###############################
## Functions for Generalized Method of Moments
###############################

#' Upper bounds on rates of double reduction
#'
#' Calculates the upper bounds of the double reduction parameters
#' according to the complete equation segregation model. See
#' Huang et. al. (2019) for details.
#'
#' @param ploidy The ploidy of the species. Should be even and at least 4.
#'
#' @return A vector of length \code{floor(ploidy/4)}. Element \code{i} is
#'     the upper bound on the probability of \code{i} pairs of
#'     identical-by-double-reduction alleles being in an individual.
#'
#' @author David Gerard
#'
#' @export
#'
#' @references
#' \itemize{
#'   \item{Huang, K., Wang, T., Dunn, D. W., Zhang, P., Cao, X., Liu, R., & Li, B. (2019). Genotypic frequencies at equilibrium for polysomic inheritance under double-reduction. \emph{G3: Genes, Genomes, Genetics}, 9(5), 1693-1706. \doi{10.1534/g3.119.400132}}
#' }
#'
#' @examples
#' drbounds(4)
#' drbounds(6)
#' drbounds(8)
#' drbounds(10)
#' drbounds(12)
#' drbounds(14)
#' drbounds(16)
#'
drbounds <- function(ploidy) {
  stopifnot(ploidy > 2)
  stopifnot(ploidy %% 2 == 0)
  ibdr <- floor(ploidy / 4)
  alpha <- rep(NA_real_, length = ibdr)
  for (i in seq_along(alpha)) {
    jvec <- i:ibdr
    alpha[[i]] <-
      exp(
        log_sum_exp(
          (ploidy / 2 - 3 * jvec) * log(2) +
            lchoose(jvec, i) +
            lchoose(ploidy / 2, jvec) +
            lchoose(ploidy / 2 - jvec, ploidy / 2 - 2 * jvec) -
            lchoose(ploidy, ploidy / 2)
          )
      )
  }
  return(alpha)
}

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
#' @param ngen The number of generations of random mating of current
#'     empirical frequencies to compare against.
#'
#' @author David Gerard
#'
#' @noRd
chisqdiv <- function(nvec,
                     alpha,
                     denom = c("expected", "observed"),
                     ngen = 1) {
  ## Check parameters ----
  ploidy <- length(nvec) - 1
  stopifnot(ploidy %% 2 == 0, ploidy > 0)
  stopifnot(length(alpha) == floor(ploidy / 4))
  stopifnot(alpha >= 0, sum(alpha) <= 1)
  stopifnot(nvec >= 0)
  stopifnot(length(ngen) == 1, ngen >= 1)
  denom <- match.arg(denom)

  ## calculate divergence ----
  n <- sum(nvec)
  qhat <- nvec / n

  fq <- qhat
  for (i in seq_len(ngen)) {
    fq <- freqnext(freq = fq, alpha = alpha)
  }

  notzero <- qhat > 10^-6 | fq > 10^-6
  if (denom == "expected") {
    chisq <- n * sum((qhat[notzero] - fq[notzero]) ^ 2 / fq[notzero])
  } else if (denom == "observed") {
    chisq <- n * sum((qhat[notzero] - fq[notzero]) ^ 2 / qhat[notzero])
  }

  return(chisq)
}

neymandiv <- function(nvec, alpha, ngen = 1) {
  chisqdiv(nvec = nvec, alpha = alpha, denom = "observed", ngen = ngen)
}

pearsondiv <- function(nvec, alpha, ngen = 1) {
  chisqdiv(nvec = nvec, alpha = alpha, denom = "expected", ngen = ngen)
}


#' G-test statistic. I.e. likelihood ratio of multinomial
#'
#' @inheritParams chisqdiv
#'
#' @author David Gerard
#'
#' @noRd
gdiv <- function(nvec, alpha, ngen = 1) {
  ## Check parameters ----
  ploidy <- length(nvec) - 1
  stopifnot(ploidy %% 2 == 0, ploidy > 0)
  stopifnot(length(alpha) == floor(ploidy / 4))
  stopifnot(alpha >= 0, sum(alpha) <= 1)
  stopifnot(nvec >= 0)

  ## Iterate through ngen generations of random mating ----
  n <- sum(nvec)
  qhat <- nvec / n

  fq <- qhat
  for (i in seq_len(ngen)) {
    fq <- freqnext(freq = fq, alpha = alpha)
  }

  ## Calculate statistic
  ei <- fq * n
  notzero <- nvec > 0
  gstat <- 2 * sum(nvec[notzero] * log(nvec[notzero] / ei[notzero]))

  return(gstat)
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
#'         American Institute of Physics. \doi{10.1063/1.3703631}}
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
#' For ploides greater than 4, not all gametic frequencies are possible
#' given rates of double reduction. That is, if the model for meiosis
#' from Huang et. al. (2019) is incorrect, then this is not an appropriate
#' test for HWE.
#'
#' @param nvec A vector containing the observed genotype counts,
#'     where \code{nvec[[i]]} is the number of individuals with genotype
#'     \code{i-1}. This should be of length \code{ploidy+1}.
#' @param obj What chi-square function should we minimize?
#'     \describe{
#'       \item{\code{pearson}}{\deqn{\sum (o-e)^2/e}}
#'       \item{\code{g}}{\deqn{2\sum o \log(o/2)}}
#'       \item{\code{neyman}}{\deqn{\sum (o-e)^2/o}}
#'     }
#'     See Berkson (1980) for a description.
#' @param ngen The number of generations of random mating of current
#'     empirical frequencies to compare against.
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
#'   \item{\code{chisq_ndr}}{The chi-square test statistic for
#'       testing against the null of no double reduction (conditional
#'       on equilibrium). In diploids, this value is \code{NULL}.}
#'   \item{\code{df_ndr}}{The degrees of freedom associated with
#'       \code{chisq_ndr}.}
#'   \item{\code{p_ndr}}{The p-value against the null of no double
#'       reduction (conditional on equilibrium). In diploids, this value
#'       is \code{NULL}.}
#' }
#'
#' @author David Gerard
#'
#' @references
#' \itemize{
#'   \item{Agresti, A., & Coull, B. A. (1998). Approximate is better than "exact" for interval estimation of binomial proportions. The American Statistician, 52(2), 119-126. \doi{10.1080/00031305.1998.10480550}}
#'   \item{Berkson, J. (1980). Minimum chi-square, not maximum likelihood!. The Annals of Statistics, 8(3), 457-487. \doi{10.1214/aos/1176345003}}
#'   \item{Huang, K., Wang, T., Dunn, D. W., Zhang, P., Cao, X., Liu, R., & Li, B. (2019). Genotypic frequencies at equilibrium for polysomic inheritance under double-reduction. G3: Genes, Genomes, Genetics, 9(5), 1693-1706. \doi{10.1534/g3.119.400132}}
#' }
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
#' nvec <- round(hwefreq(r = 0.5,
#'                       alpha = 0.2,
#'                       ploidy = 6,
#'                       niter = 100,
#'                       tol = -Inf) * 100)
#' hwemom(nvec)
#'
#' ## Octoploid case with exact frequencies at HWE
#' nvec <- round(hwefreq(r = 0.5,
#'                       alpha = c(0.1, 0.01),
#'                       ploidy = 8,
#'                       niter = 100,
#'                       tol = -Inf) * 1000)
#' hwemom(nvec)
#'
#' @noRd
hwemom <- function(nvec,
                   obj = c("g", "pearson", "neyman"),
                   ngen = 1) {
  ## Check parameters ----
  ploidy <- length(nvec) - 1
  stopifnot(ploidy %% 2 == 0, ploidy > 0)
  stopifnot(nvec >= 0)
  stopifnot(is.vector(nvec))
  stopifnot(length(ngen) == 1, ngen >= 1)
  ibdr <- floor(ploidy / 4)
  obj <- match.arg(obj)

  ## Return NA's if too few counts ----
  if (sum(nvec > 0.5) < ibdr + 1) {
    retlist <- list()

    retlist <- list(chisq_hwe = NA_real_,
                    df_hwe = NA_real_,
                    p_hwe = NA_real_,
                    alpha = rep(NA_real_, times = ibdr),
                    chisq_ndr = NA_real_,
                    df_ndr = NA_real_,
                    p_ndr = NA_real_)

    retlist$r <- sum(0:ploidy * nvec) / (ploidy * sum(nvec))

    return(retlist)
  }

  ## Choose objective function ----
  if (obj == "pearson") {
    divfun <- pearsondiv
    if (ngen == 1) {
      grfun <- dpearsondiv_dalpha
    } else {
      grfun <- NULL
    }
  } else if (obj == "g") {
    divfun <- gdiv
    if (ngen == 1) {
      grfun <- dgdiv_dalpha
    } else {
      grfun <- NULL
    }
  } else if (obj == "neyman") {
    divfun <- neymandiv
    if (ngen == 1) {
      grfun <- dneymandiv_dalpha
    } else {
      grfun <- NULL
    }
  }

  ## tell folks to use other stuff ----
  if (ploidy == 4) {
    message("You should use `hwetetra()` for tetraploids")
  } else if (ploidy == 2) {
    message("Don't use this function. There are far better packages for diploids.")
  }

  minval <- sqrt(.Machine$double.eps) ## minimum DR to explore
  ## optimize ----
  if (ibdr == 0) {
    ## Diploid: Just return chi-square
    chisq <- divfun(nvec = nvec, alpha = NULL, ngen = ngen)
    pval <- stats::pchisq(q = chisq,
                          df = ploidy - ibdr - 1,
                          lower.tail = FALSE)
    phat <- sum(0:ploidy * nvec) / (ploidy * sum(nvec))
    retlist <- list(p = c(1 - phat, phat),
                    chisq_hwe = chisq,
                    df_hwe = ploidy - ibdr - 1,
                    p_hwe = pval,
                    alpha = NULL,
                    chisq_ndr = NULL,
                    df_ndr = NULL,
                    p_ndr = NULL)
  } else if (ibdr == 1) {
    ## Tetraploid or Hexaploid: Use Brent's method
    upper_alpha <- drbounds(ploidy = ploidy)
    oout <- stats::optim(par = minval,
                         fn = divfun,
                         method = "Brent",
                         lower = minval,
                         upper = upper_alpha,
                         nvec = nvec,
                         ngen = ngen)
    pval_hwe <- stats::pchisq(q = oout$value,
                              df = ploidy - ibdr - 1,
                              lower.tail = FALSE)
    chisq_null <- divfun(nvec = nvec,
                         alpha = rep(0, length.out = ibdr),
                         ngen = ngen)
    chisq_alpha <- (chisq_null - oout$value)
    pval_alpha <- stats::pchisq(q = chisq_alpha,
                                df = ibdr,
                                lower.tail = FALSE) / 2
    retlist <- list(chisq_hwe = oout$value,
                    df_hwe = ploidy - ibdr - 1,
                    p_hwe = pval_hwe,
                    alpha = oout$par,
                    chisq_ndr = chisq_alpha,
                    df_ndr = ibdr,
                    p_ndr = pval_alpha)
  } else {
    ## Higher Ploidy: Use L-BFGS-B using boundaries from CES model
    upper_alpha <- drbounds(ploidy = ploidy)
    oout <- stats::optim(par = rep(minval, ibdr),
                         fn = divfun,
                         gr = grfun,
                         method = "L-BFGS-B",
                         lower = rep(minval, ibdr),
                         upper = upper_alpha,
                         nvec = nvec,
                         ngen = ngen)
    alpha <- oout$par
    pval_hwe <- stats::pchisq(q = oout$value,
                              df = ploidy - ibdr - 1,
                              lower.tail = FALSE)
    chisq_null <- divfun(nvec = nvec,
                         alpha = rep(0, length.out = ibdr),
                         ngen = ngen)
    chisq_alpha <- (chisq_null - oout$value)
    pval_alpha <- stats::pchisq(q = chisq_alpha,
                                df = ibdr,
                                lower.tail = FALSE) / 2
    retlist <- list(chisq_hwe = oout$value,
                    df_hwe = ploidy - ibdr - 1,
                    p_hwe = pval_hwe,
                    alpha = alpha,
                    chisq_ndr = chisq_alpha,
                    df_ndr = ibdr,
                    p_ndr = pval_alpha)
  }

  ## get estimate of allele frequency ----
  retlist$r <- sum(0:ploidy * nvec) / (ploidy * sum(nvec))

  return(retlist)
}
