##################################
## U-statistic approach to HWE and DR
##################################

#' Objective used in \code{\link{hweustat}()}
#'
#' @param nvec A vector of length ploidy + 1. \code{nvec[i]} is the observed
#'     number of individuals with dosage \code{i-1}.
#' @param alpha The double reduction parameter. Should be of length
#'     \code{floor(ploidy/4)}. These should sum to less than 1, since
#'     \code{1-sum(alpha)} is the assumed probability of no double
#'     reduction.
#' @param omega A (ploidy+1) by (ploidy+1) positive semidefinite matrix.
#'     This is an estimate of the \emph{inverse} covariance of the elements
#'     of the U-statistic.
#' @param which_keep A logical vector. Which genotypes do we use in the
#'     objective function evaluation?
#'
#' @author David Gerard
#'
#' @noRd
uobj <- function(nvec, alpha, omega = NULL, which_keep = NULL) {
  ## check parameters ----
  ploidy <- length(nvec) - 1
  stopifnot(ploidy %% 2 == 0, ploidy > 0)
  stopifnot(length(alpha) == floor(ploidy / 4))
  stopifnot(nvec >= 0)
  if (is.null(which_keep)) {
    which_keep <- rep(TRUE, ploidy + 1)
  }
  stopifnot(length(which_keep) == ploidy + 1)
  stopifnot(is.logical(which_keep))
  numkeep <- sum(which_keep)
  if (!is.null(omega)) {
    stopifnot(dim(omega) == c(sum(numkeep), sum(numkeep)))
  }


  ## Calculate objective ----
  n <- sum(nvec)
  qhat <- nvec / n
  fq <- freqnext(freq = qhat, alpha = alpha, check = FALSE)

  if (is.null(omega)) {
    diff <- (qhat - fq)[which_keep]
    return(sum(diff^2) * n)
  } else {
    diff <- matrix((qhat - fq)[which_keep], ncol = 1)
    return(c(t(diff) %*% omega %*% diff) * n)
  }
}

#' Generalized inverse of a symmetric p.d. matrix
#'
#' @author David Gerard
#'
#' @noRd
ginv <- function(omega) {
  stopifnot(is.matrix(omega))
  stopifnot(nrow(omega) == ncol(omega))

  eout <- eigen(omega, symmetric = TRUE)

  stopifnot(eout$values > -sqrt(.Machine$double.eps))

  which_pos <- eout$values > sqrt(.Machine$double.eps)

  f <- c(1 / eout$values[which_pos], rep(1, length.out = sum(!which_pos)))

  return(eout$vectors %*% diag(f) %*% t(eout$vectors))
}

#' \emph{Inverse} covariance estimate of elements of ustat
#'
#'
#' @inheritParams uobj
#'
#' @return \emph{Generalized} inverse covariance estimate.
#'
#' @author David Gerard
#'
#' @noRd
ucov <- function(nvec, alpha, which_keep = NULL) {
  ## Check parameters ----
  ploidy <- length(nvec) - 1
  stopifnot(ploidy %% 2 == 0, ploidy > 0)
  stopifnot(length(alpha) == floor(ploidy / 4))
  stopifnot(nvec >= 0)
  if (is.null(which_keep)) {
    which_keep <- rep(TRUE, ploidy + 1)
  }
  stopifnot(length(which_keep) == ploidy + 1)
  stopifnot(is.logical(which_keep))

  ## Calculate covariance ----
  n <- sum(nvec)
  qhat <- nvec / n
  fq <- freqnext(freq = qhat, alpha = alpha, check = FALSE)
  A <- zsegarray(alpha = alpha, ploidy = ploidy)

  omega <- diag(qhat) - tcrossprod(qhat, fq) - tcrossprod(fq, qhat)

  indset <- 0:ploidy
  for (p in indset) {
    for (r in indset) {
      for (ell in indset) {
        for (m in indset) {
          omega[p + 1, r + 1] <- omega[p + 1, r + 1] +
            A[ell + 1, m + 1, p + 1] * A[ell + 1, m + 1, r + 1] * qhat[[ell + 1]] * qhat[[m + 1]]
        }
      }
    }
  }

  ## return generalized inverse ----
  ominv <- ginv(omega[which_keep, which_keep, drop = FALSE])

  return(ominv)
}

#' U-process minimizer approach to equilibrium testing and double reduction
#' estimation
#'
#' Estimates double reduction and tests for equilibrium while accounting
#' for double reduction. It does this using an approach called
#' "U-process minimization", where we minimize a function of a U-statistic
#' that should be 0 at equilibrium given the true double reduction rate.
#'
#' This is a two-step estimator, where we first obtain a consistent
#' estimate of the double reduction parameter, use this to
#' estimate the covariance of estimators, then use this to obtain
#' our final estimate of the double reduction parameter.
#'
#' @param nvec A vector containing the observed genotype counts,
#'     where \code{nvec[[i]]} is the number of individuals with genotype
#'     \code{i-1}. This should be of length \code{ploidy+1}.
#' @param thresh_mult The threshhold for ignoring the genotype. We keep
#'     genotypes such that \code{max(nvec) / nvec <= thresh_mult}.
#'     Setting this to \code{Inf} uses all genotypes.
#' @param thresh_tot The threshhold for ignoring the genotype. We keep
#'     genotypes such that \code{nvec >= thresh_tot}.
#'     Setting this to \code{0} uses all genotypes.
#'
#' @return A list with some or all of the following elements:
#' \describe{
#'   \item{\code{alpha}}{The estimated double reduction parameter(s).
#'       In diploids, this value is \code{NULL}.}
#'   \item{\code{chisq_hwe}}{The chi-square test statistic for testing
#'       against the null of equilibrium.}
#'   \item{\code{df_hwe}}{The degrees of freedom associated with
#'       \code{chisq_hwe}.}
#'   \item{\code{p_hwe}}{The p-value against the null of equilibrium.}
#' }
#'
#' @author David Gerard
#'
#' @examples
#' set.seed(99)
#' ploidy <- 6
#' size <- 100
#' r <- 0.5
#' alpha <- 0.1
#' qvec <- hwefreq(r = r, alpha = alpha, ploidy = ploidy)
#' nvec <- c(rmultinom(n = 1, size = size, prob = qvec))
#' hweustat(nvec = nvec)
#'
#' @export
hweustat <- function(nvec, thresh_mult = 100, thresh_tot = 10) {
  ploidy <- length(nvec) - 1
  stopifnot(ploidy %% 2 == 0, ploidy >= 4)
  ibdr <- floor(ploidy / 4)
  stopifnot(thresh_mult >= 0, length(thresh_mult) == 1)
  stopifnot(thresh_tot >= 0, length(thresh_tot) == 1)
  minval <- sqrt(.Machine$double.eps)

  which_keep <- ((max(nvec) / nvec) <= thresh_mult) & (nvec >= thresh_tot)

  if (sum(which_keep) < ibdr + 1) {
    return(
      list(
        alpha = NA_real_,
        chisq_hwe = NA_real_,
        df_hwe = NA_real_,
        p_hwe = NA_real_
      )
    )
  }

  if (ibdr == 1) {
    upper_alpha <- drbounds(ploidy = ploidy)

    ## step 1 ----
    oout <- stats::optim(par = minval,
                         fn = uobj,
                         method = "Brent",
                         lower = minval,
                         upper = upper_alpha,
                         nvec = nvec,
                         omega = NULL,
                         which_keep = which_keep)

    ## Calculate covariance ----
    alpha_tilde <- oout$par
    omega <- ucov(nvec = nvec, alpha = alpha_tilde, which_keep = which_keep)

    ## step 2 ----
    oout <- stats::optim(par = minval,
                         fn = uobj,
                         method = "Brent",
                         lower = minval,
                         upper = upper_alpha,
                         nvec = nvec,
                         omega = omega,
                         which_keep = which_keep)
  } else {
    upper_alpha <- drbounds(ploidy = ploidy)

    ## step 1 ----
    oout <- stats::optim(par = rep(minval, ibdr),
                         fn = uobj,
                         method = "L-BFGS-B",
                         lower = rep(minval, ibdr),
                         upper = upper_alpha,
                         nvec = nvec,
                         omega = NULL,
                         which_keep = which_keep)

    ## Calculate covariance ----
    alpha_tilde <- oout$par
    omega <- ucov(nvec = nvec, alpha = alpha_tilde, which_keep = which_keep)

    ## step 2 ----
    oout <- stats::optim(par = rep(minval, ibdr),
                         fn = uobj,
                         method = "L-BFGS-B",
                         lower = rep(minval, ibdr),
                         upper = upper_alpha,
                         nvec = nvec,
                         omega = omega,
                         which_keep = which_keep)
  }

  ## Calculate degrees of freedom ----
  df_hwe <- sum(which_keep) -
    ibdr -
    sum(eigen(omega)$values - 1 < sqrt(.Machine$double.eps)) -
    1

  ## return ----
  retlist <- list(
    alpha = oout$par,
    chisq_hwe = oout$value,
    df_hwe = df_hwe
  )
  retlist$p_hwe <- stats::pchisq(q = retlist$chisq_hwe,
                                 df = retlist$df_hwe,
                                 lower.tail = FALSE)
  retlist$p_hwe[retlist$df_hwe <= 0] <- NA_real_
  return(retlist)
}
