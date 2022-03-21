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
#' @param addval The penalization applied to each genotype for the random
#'     mating hypothesis. This corresponds to the number of counts each
#'     genotype has \emph{a priori}.
#' @param which_keep A logical vector. The \code{FALSE}'s will be aggregated
#'     into one group.
#'
#' @author David Gerard
#'
#' @noRd
llike <- function(nvec, pvec, addval = 0, which_keep = NULL) {

  ## Check input ----
  ploidy <- length(nvec) - 1
  stopifnot(ploidy %% 2 == 0)
  stopifnot(length(pvec) == ploidy / 2 + 1)
  stopifnot(abs(sum(pvec) - 1) < 10^-6)
  stopifnot(nvec >= 0)
  stopifnot(addval >= 0, length(addval) == 1)
  if (is.null(which_keep)) {
    which_keep <- rep(TRUE, ploidy + 1)
  }

  ## Calculate frequencies ----
  q <- stats::convolve(pvec, rev(pvec), type = "open")
  if (!all(which_keep)) {
    nvec <- c(nvec[which_keep], sum(nvec[!which_keep]))
    q <- c(q[which_keep], sum(q[!which_keep]))
  }

  ## Calculate likelihood ----
  stats::dmultinom(x = nvec,
                   prob = q,
                   log = TRUE) +
    addval * sum(log(pvec[pvec > sqrt(.Machine$double.eps)]))
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
#' @param addval The penalization applied to each genotype for the random
#'     mating hypothesis. This corresponds to the number of counts each
#'     genotype has \emph{a priori}.
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
rmem <- function(nvec, tol = 10^-3, maxit = 100, addval = 1 / 100) {
  ploidy <- length(nvec) - 1
  stopifnot(ploidy %% 2 == 0)
  stopifnot(nvec >= 0)
  stopifnot(addval >= 0, length(addval) == 1)

  ## Initialize under assumption of random mating with alpha = 0 (small penalty) ----
  pvec <- stats::dbinom(x = 0:(ploidy / 2),
                        size = ploidy / 2,
                        prob = sum(nvec * 0:ploidy) / (sum(nvec) * ploidy)) + addval
  pvec <- pvec / sum(pvec)
  ll <- llike(nvec = nvec, pvec = pvec, addval = addval)

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

    ## Add penalty to etavec ----
    etavec <- etavec + addval

    ## normalize ----
    pvec <- etavec / sum(etavec)

    ## calculate (penalized) log-likelihood ----
    ll <- llike(nvec = nvec, pvec = pvec, addval = addval)

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
#' @param nvec A vector containing the observed genotype counts,
#'     where \code{nvec[[i]]} is the number of individuals with genotype
#'     \code{i-1}. This should be of length \code{ploidy+1}.
#' @param thresh All groups with counts less than \code{nvec} will
#'     be aggregated together.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{p}}{The estimated gamete genotype frequencies. \code{p[[i]]}
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
rmlike <- function(nvec, thresh = 1) {
  ploidy <- length(nvec) - 1
  stopifnot(ploidy %% 2 == 0)
  stopifnot(nvec >= 0)
  stopifnot(is.vector(nvec))
  stopifnot(length(thresh) == 1, thresh >= 0)

  if (ploidy == 2) {
    message("Don't use this function. There are far better packages for diploids.")
  }

  ## Choose which groups to aggregate ----
  which_keep <- nvec >= thresh
  if (sum(which_keep) <= ploidy / 2) {
    which_keep <- rep(TRUE, ploidy + 1)
    which_keep[order(nvec, decreasing = FALSE)[1:(ploidy / 2)]] <- FALSE
  }
  if (sum(which_keep) >= ploidy) {
    ## aggregating one group = no aggregation at all.
    which_keep <- rep(TRUE, ploidy + 1)
  }

  ## Estimate gametic frequencies ----
  hout <- rmem(nvec = nvec)
  names(hout) <- 0:(ploidy / 2)

  ## Get log-likelihood under null ----
  ll0 <- llike(nvec = nvec, pvec = hout, which_keep = which_keep)

  ## Get log-likelihood under alternative ----
  qmle <- nvec / sum(nvec)
  if (!all(which_keep)) {
    nmle <- c(nvec[which_keep], sum(nvec[!which_keep]))
    qmle <- c(qmle[which_keep], sum(qmle[!which_keep]))
  } else {
    nmle <- nvec
  }
  lla <- stats::dmultinom(x = nmle, prob = qmle, log = TRUE)

  ## Get degrees of freedom ---
  if (!all(which_keep)) {
    df_rm <- sum(which_keep) - ploidy / 2 + sum(hout < 0.001)
  } else {
    df_rm <- ploidy / 2
  }

  ## Run likelihood ratio test against null of random mating ----
  chisq <- -2 * (ll0 - lla)
  pval <- stats::pchisq(q = chisq, df = df_rm, lower.tail = FALSE)

  if (df_rm <= 0) {
    pval <- NA_real_
  }

  retlist <- list(p = hout,
                  chisq_rm = chisq,
                  df_rm = df_rm,
                  p_rm = pval)

  return(retlist)
}
