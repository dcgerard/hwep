############################
## Functions for likelihood inference
############################

#' Log-likelihood of gamete frequencies under random mating
#'
#' @param nvec Of a vector of counts. \code{nvec[[i]]} is the number of
#'     individuals that have genotype \code{i-1}. The ploidy is assumed
#'     to be \code{length(nvec)-1}.
#' @param pvec The vector of gamete frequencies. \code{pvec[[i]]} is the
#'     probability that a gamete will have dosage \code{i-1}. This should
#'     be of length \code{ploidy/2 + 1}.
#'
#' @author David Gerard
#'
#' @noRd
llike <- function(nvec, pvec) {
  ploidy <- length(nvec) - 1
  stopifnot(ploidy %% 2 == 0)
  stopifnot(length(pvec) == ploidy / 2 + 1)
  stopifnot(abs(sum(pvec) - 1) < 10^-6)
  stopifnot(nvec >= 0)

  stats::dmultinom(x = nvec,
                   prob = stats::convolve(pvec, rev(pvec), type = "open"),
                   log = TRUE)
}

#' Estimate gametic proportions under random mating
#'
#' This uses an EM algorithm.
#'
#' @param nvec Of a vector of counts. \code{nvec[[i]]} is the number of
#'     individuals that have genotype \code{i-1}. The ploidy is assumed
#'     to be \code{length(nvec)-1}.
#' @param tol The stopping criterion tolerance.
#' @param maxit The maximum number of EM iterations to run.
#'
#' @author David Gerard
#'
#' @examples
#' ## Gamete frequencies
#' pvec <- stats::runif(6)
#' pvec <- pvec / sum(pvec)
#'
#' ## Genotype frequencies
#' qvec <- stats::convolve(pvec, rev(pvec), type = "open")
#'
#' ## Generate data
#' nvec <- round(100000000 * qvec)
#'
#' ## Estimate pvec
#' rmem(nvec = nvec)
#' pvec
#'
#' @noRd
rmem <- function(nvec, tol = 10^-3, maxit = 100) {
  ploidy <- length(nvec) - 1
  stopifnot(ploidy %% 2 == 0)
  stopifnot(nvec >= 0)

  ## Initialize under assumption of random mating with alpha = 0 ----
  pvec <- stats::dbinom(x = 0:(ploidy / 2),
                        size = ploidy / 2,
                        prob = sum(nvec * 0:ploidy) / (sum(nvec) * ploidy))
  ll <- llike(nvec = nvec, pvec = pvec)

  ## Initialize parameters -----
  paramdf <- as.data.frame(which(upper.tri(matrix(nrow = ploidy / 2 + 1, ncol = ploidy / 2 + 1), diag = TRUE), arr.ind = TRUE))
  names(paramdf) <- c("i", "j")
  paramdf$geno <- paramdf$i + paramdf$j - 2
  paramdf$w <- NA_real_
  paramdf$xi <- NA_real_

  etavec <- rep(NA_real_, length = ploidy / 2 + 1)

  ## Run EM algorthm
  i <- 1
  err <- Inf
  while (i < maxit && err > tol) {
    llold <- ll
    ## One fixed point iteration ----
    paramdf$w <- pvec[paramdf$i] * pvec[paramdf$j]
    paramdf$w[paramdf$i != paramdf$j] <- 2 * paramdf$w[paramdf$i != paramdf$j]
    sumout <- by(data = paramdf, INDICES = paramdf$geno, FUN = function(x) sum(x$w), simplify = TRUE)
    paramdf$w <- paramdf$w / as.vector(sumout)[match(paramdf$geno, names(sumout))]
    paramdf$xi <- paramdf$w * nvec[paramdf$geno + 1]
    paramdf$xi[paramdf$i == paramdf$j] <- 2 * paramdf$xi[paramdf$i == paramdf$j]

    for (j in seq_len(ploidy / 2 + 1)) {
      etavec[[j]] <- sum(paramdf$xi[paramdf$i == j | paramdf$j == j])
    }

    pvec <- etavec / sum(etavec)

    ## calculate log-likelihood ----
    ll <- llike(nvec = nvec, pvec = pvec)

    if (ll - llold < -10^-6) {
      stop("rmem: log-likelihood is not increasing")
    }

    err <- ll - llold

    i <- i + 1
  }

  return(pvec)
}


#' Likelihood inference for random mating
#'
#' Estimates gamete genotype frequencies using a maximum likelihood approach
#' and runs a likelihood ratio test for random mating.
#'
#' Let \code{q} be the genotype frequencies. Let \code{p} be the gamete
#' frequencies. Then random mating occurs if
#' \code{q == stats::convolve(p, rev(p), type = "open")}. We test for
#' this hypothesis using likelihood inference, while estimating \code{p}.
#'
#' @inheritParams hwemom
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{p}}{The estimated gamete genotype frequencies. \code{p[[i}
#'       is the estimated frequency for gamete genotype \code{i-1}.}
#'   \item{\code{chisq_rm}}{The likelihood ratio test statistic for testing
#'       against the null of random mating.}
#'   \item{\code{df_rm}}{The degrees of freedom associated with
#'       \code{chisq_rm}.}
#'   \item{\code{p_rm}}{The p-value against the null of random mating.}
#' }
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' ## Randomly generate gamete frequencies
#' set.seed(1)
#' ploidy <- 10
#' pvec <- stats::runif(ploidy / 2 + 1)
#' pvec <- pvec / sum(pvec)
#'
#' ## Genotype frequencies from gamete frequencies under random mating
#' qvec <- stats::convolve(pvec, rev(pvec), type = "open")
#'
#' ## Generate data
#' nvec <- c(stats::rmultinom(n = 1, size = 100, prob = qvec))
#'
#' ## Run rmlike()
#' rmlike(nvec = nvec)
#'
rmlike <- function(nvec) {
  ploidy <- length(nvec) - 1
  stopifnot(ploidy %% 2 == 0)
  stopifnot(nvec >= 0)

  if (ploidy == 4) {
    message("You should use `hwetetra()` for tetraploids")
  } else if (ploidy == 2) {
    message("Don't use this function. There are far better packages for diploids.")
  }

  ## Estimate gametic frequencies ----
  hout <- rmem(nvec = nvec)
  names(hout) <- 0:(ploidy / 2)

  ## Get log-likelihood under null ----
  ll0 <- llike(nvec = nvec, pvec = hout)

  ## Get log-likelihood under alternative ----
  qmle <- nvec / sum(nvec)
  lla <- stats::dmultinom(x = nvec, prob = qmle, log = TRUE)

  ## Run likelihood ratio test against null of random mating ----
  chisq <- -2 * (ll0 - lla)
  pval <- stats::pchisq(q = chisq, df = ploidy / 2, lower.tail = FALSE)

  retlist <- list(p = hout,
                  chisq_rm = chisq,
                  df_rm = ploidy / 2,
                  p_rm = pval)

  return(retlist)
}








