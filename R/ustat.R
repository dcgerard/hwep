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
#'
#' @author David Gerard
#'
#' @export
uobj <- function(nvec, alpha, omega = NULL) {
  ## check parameters ----
  ploidy <- length(nvec) - 1
  stopifnot(ploidy %% 2 == 0, ploidy > 0)
  stopifnot(length(alpha) == floor(ploidy / 4))
  stopifnot(alpha >= 0, sum(alpha) <= 1)
  stopifnot(nvec >= 0)
  if (!is.null(omega)) {
    stopifnot(dim(omega) == c(ploidy + 1, ploidy + 1))
  }

  ## Calculate objective ----
  n <- sum(nvec)
  qhat <- nvec / n
  fq <- freqnext(freq = qhat, alpha = alpha)

  if (is.null(omega)) {
    return(sum((qhat - fq)^2) * n)
  } else {
    diff <- matrix(qhat - fq, ncol = 1)
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
ucov <- function(nvec, alpha) {
  ## Check parameters ----
  ploidy <- length(nvec) - 1
  stopifnot(ploidy %% 2 == 0, ploidy > 0)
  stopifnot(length(alpha) == floor(ploidy / 4))
  stopifnot(alpha >= 0, sum(alpha) <= 1)
  stopifnot(nvec >= 0)

  ## Calculate covariance ----
  n <- sum(nvec)
  qhat <- nvec / n
  fq <- freqnext(freq = qhat, alpha = alpha)
  A <- zsegarray(alpha = alpha, ploidy = ploidy)

  omega <- diag(qhat) - tcrossprod(qhat, fq) - tcrossprod(fq, qhat)
  for (p in 0:ploidy) {
    for (r in 0:ploidy) {
      for (ell in 0:ploidy) {
        for (m in 0:ploidy) {
          omega[p + 1, r + 1] <- omega[p + 1, r + 1] +
            A[ell + 1, m + 1, p + 1] * A[ell + 1, m + 1, r + 1] * qhat[[ell + 1]] * qhat[[m + 1]]
        }
      }
    }
  }

  ## return generalized inverse ----
  return(ginv(omega))
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
#' @export
hweustat <- function(nvec) {
  ploidy <- length(nvec) - 1
  stopifnot(ploidy %% 2 == 0, ploidy >= 4)
  ibdr <- floor(ploidy / 4)
  minval <- sqrt(.Machine$double.eps)

  if (ibdr == 1) {
    upper_alpha <- drbounds(ploidy = ploidy)

    ## step 1 ----
    oout <- stats::optim(par = minval,
                         fn = uobj,
                         method = "Brent",
                         lower = minval,
                         upper = upper_alpha,
                         nvec = nvec,
                         omega = NULL)

    ## Calculate covariance ----
    alpha_tilde <- oout$par
    omega <- ucov(nvec = nvec, alpha = alpha_tilde)

    ## step 2 ----
    oout <- stats::optim(par = minval,
                         fn = uobj,
                         method = "Brent",
                         lower = minval,
                         upper = upper_alpha,
                         nvec = nvec,
                         omega = omega)
  } else {
    ## step 1 ----
    upper_alpha <- drbounds(ploidy = ploidy)
    oout <- stats::optim(par = rep(minval, ibdr),
                         fn = uobj,
                         method = "L-BFGS-B",
                         lower = rep(minval, ibdr),
                         upper = upper_alpha,
                         nvec = nvec,
                         omega = NULL)

    ## Calculate covariance ----
    alpha_tilde <- oout$par
    omega <- ucov(nvec = nvec, alpha = alpha_tilde)

    ## step 2 ----
    upper_alpha <- drbounds(ploidy = ploidy)
    oout <- stats::optim(par = rep(minval, ibdr),
                         fn = uobj,
                         method = "L-BFGS-B",
                         lower = rep(minval, ibdr),
                         upper = upper_alpha,
                         nvec = nvec,
                         omega = omega)
  }


  ## return ----
  retlist <- list(
    alpha = oout$par,
    chisq_hwe = oout$value,
    df_hwe = ploidy - ibdr - 1
  )
  retlist$p_hwe <- stats::pchisq(q = retlist$chisq_hwe,
                                 df = retlist$df_hwe,
                                 lower.tail = FALSE)

  return(retlist)
}
