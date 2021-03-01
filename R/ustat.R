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
#' @param which_keep A lotical vector of length \code{ploidy + 1}. All of the
#'     \code{FALSE}'s will be combined with each other. They form the last
#'     group.
#' @param omega A \code{sum(which_keep)+1} by \code{sum(which_keep)+1}
#'     positive semidefinite matrix.
#'     This is an estimate of the \emph{generalized inverse} covariance
#'     of the elements of the U-statistic. The last row and column
#'     correspond to all aggregated groups determined by \code{which_keep}.
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
  stopifnot(is.logical(which_keep), length(which_keep) == ploidy + 1)
  numkeep <- sum(which_keep)
  if (!is.null(omega)) {
    if (numkeep == ploidy + 1) {
      stopifnot(dim(omega) == c(ploidy + 1, ploidy + 1))
    } else {
      stopifnot(dim(omega) == c(numkeep + 1, numkeep + 1))
    }
  }

  ## Calculate objective ----
  n <- sum(nvec)
  qhat <- nvec / n
  fq <- freqnext(freq = qhat, alpha = alpha, check = FALSE)

  if (numkeep != ploidy + 1) {
    qhat <- c(qhat[which_keep], sum(qhat[!which_keep]))
    fq   <- c(fq[which_keep], sum(fq[!which_keep]))
  }

  if (is.null(omega)) {
    diff <- (qhat - fq)
    return(sum(diff^2) * n)
  } else {
    diff <- matrix((qhat - fq), ncol = 1)
    return(c(t(diff) %*% omega %*% diff) * n)
  }
}

#' Generalized inverse of a symmetric p.d. matrix
#'
#' This is the Moore-penrose inverse of a \emph{symmetric positive definite}
#' matrix.
#'
#' @param omega A matrix to invert. Must be symmetric and positive definite.
#' @param tol The tolerance on the eigenvalues to discard.
#'
#' @author David Gerard
#'
#' @noRd
ginv <- function(omega, tol = sqrt(.Machine$double.eps)) {
  stopifnot(is.matrix(omega))
  stopifnot(nrow(omega) == ncol(omega))

  eout <- eigen(omega, symmetric = TRUE)

  stopifnot(eout$values > -tol)

  which_pos <- eout$values > tol

  f <- c(1 / eout$values[which_pos], rep(0, length.out = sum(!which_pos)))

  return(list(mat = eout$vectors %*% diag(f) %*% t(eout$vectors),
              rank = sum(which_pos)))
}


#' Same as \code{\link{ginv}()}, but explicitly specify rank
#'
#' @param Q The asymptotic covariance matrix, or an estimate of such matrix.
#' @param df The degrees of freedom to use.
#'
#' @author David Gerard
#'
#' @noRd
projme <- function(Q, df = nrow(Q) - 1) {
  stopifnot(is.matrix(Q))
  K <- nrow(Q)
  stopifnot(nrow(Q) == ncol(Q))
  stopifnot(df >= 0, df <= K)

  eout <- eigen(Q, symmetric = TRUE)
  f <- c(1 / eout$values[1:df], rep(0, K - df))
  return(mat = eout$vectors %*% diag(f) %*% t(eout$vectors))
}

#' Covariance estimate of elements of ustat
#'
#' This is an empirical one, that doesn't work as well as ucov in practice
#'
#' @inheritParams uobj
#'
#' @return Covariance estimate.
#'
#' @author David Gerard
#'
#' @noRd
ucov_empirical <- function(nvec, alpha) {
  ## Check parameters ----
  ploidy <- length(nvec) - 1
  stopifnot(ploidy %% 2 == 0, ploidy > 0)
  stopifnot(length(alpha) == floor(ploidy / 4))
  stopifnot(nvec >= 0)

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

  return(omega)
}


#' \emph{Inverse} covariance estimate of elements of ustat
#'
#' This is the theoretical one, based on U-statistic theory.
#'
#' @inheritParams uobj
#'
#' @return Covariance estimate.
#'
#' @author David Gerard
#'
#' @noRd
ucov <- function(nvec, alpha) {
  ## Check parameters ----
  ploidy <- length(nvec) - 1
  stopifnot(ploidy %% 2 == 0, ploidy > 0)
  stopifnot(length(alpha) == floor(ploidy / 4))
  stopifnot(nvec >= 0)

  ## Calculate covariance ----
  n <- sum(nvec)
  qhat <- nvec / n
  fq <- freqnext(freq = qhat, alpha = alpha, check = FALSE)
  A <- aperm(zsegarray(alpha = alpha, ploidy = ploidy), c(3, 1, 2)) ## put offspring in first index

  Qmat <- diag(fq) - tcrossprod(fq)

  bread <- diag(ploidy + 1) -
    2 * tensr::mat(
      A = tensr::amprod(A = A, M = matrix(fq, nrow = 1), k = 2),
      k = 3
      )

  omega <- t(bread) %*% Qmat %*% bread

  return(omega)
}

#' \emph{Inverse} covariance estimate of elements of ustat
#'
#' This is the naive one, based only on the multinomial distribution
#' and ignoring that we are estimating q in f(q). Does not work well.
#'
#' @inheritParams uobj
#'
#' @return Covariance estimate.
#'
#' @author David Gerard
#'
#' @noRd
ucov_naive  <- function(nvec, alpha) {
  ## Check parameters ----
  ploidy <- length(nvec) - 1
  stopifnot(ploidy %% 2 == 0, ploidy > 0)
  stopifnot(length(alpha) == floor(ploidy / 4))
  stopifnot(nvec >= 0)

  n <- sum(nvec)
  qhat <- nvec / n
  fq <- freqnext(freq = qhat, alpha = alpha, check = FALSE)

  Qmat <- diag(fq) - outer(X = fq, Y = fq, FUN = `*`)

  return(Qmat)
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
#' @param thresh The threshhold for ignoring the genotype. We keep
#'     genotypes such that \code{nvec >= thresh}.
#'     Setting this to \code{0} uses all genotypes.
#' @param effdf A logical. Should we use the ad-hoc
#'     "effective degrees of freedom" (\code{TRUE}) or not (\code{FALSE})?
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
hweustat <- function(nvec,
                     thresh = 1,
                     effdf = TRUE) {
  ploidy <- length(nvec) - 1
  stopifnot(ploidy %% 2 == 0, ploidy >= 4)
  ibdr <- floor(ploidy / 4)
  stopifnot(thresh >= 0, length(thresh) == 1)
  stopifnot(is.logical(effdf), length(effdf) == 1)
  minval <- sqrt(.Machine$double.eps)

  ## Choose which groups to aggregate ----
  which_keep <- nvec >= thresh
  if (sum(!which_keep) >= 1) {
    which_keep[which_keep][which.min(nvec[which_keep])] <- FALSE
  }

  ## Return early if too few groups ----
  if (sum(which_keep) - ibdr <= 0) {
    return(
      list(
        alpha = rep(NA_real_, ibdr),
        chisq_hwe = NA_real_,
        df_hwe = NA_real_,
        p_hwe = NA_real_
      )
    )
  }

  ## Create aggregation matrix ----
  numkeep <- sum(which_keep)
  H <- aggfun(which_keep = which_keep)
  if (numkeep >= ploidy) {
    ## aggregating one group = no aggregation at all.
    which_keep <- rep(TRUE, ploidy + 1)
    numkeep <- ploidy + 1
  }

  ## Upper bounds on alpha ----
  upper_alpha <- drbounds(ploidy = ploidy)

  ## optimization method ----
  ometh <- ifelse(ibdr == 1, "Brent", "L-BFGS-B")

  ## Run two-step procedure ---------------------------------------------------

  ## step 1 ----
  oout <- stats::optim(par = rep(minval, ibdr),
                       fn = uobj,
                       method = ometh,
                       lower = rep(minval, ibdr),
                       upper = upper_alpha,
                       nvec = nvec,
                       omega = NULL,
                       which_keep = which_keep)

  ## Calculate covariance ----
  alpha_tilde <- oout$par
  covmat <- ucov(nvec = nvec, alpha = alpha_tilde)
  projout <- ginv(H %*% covmat %*% t(H))
  omega <- projout$mat

  ## step 2 ----
  oout <- stats::optim(par = rep(minval, ibdr),
                       fn = uobj,
                       method = ometh,
                       lower = rep(minval, ibdr),
                       upper = upper_alpha,
                       nvec = nvec,
                       omega = omega,
                       which_keep = which_keep)

  ## Get alpha ----
  alpha <- oout$par

  ## Calculate chi-square statistic ---
  chisq_hwe <- oout$value

  ## Calculate degrees of freedom ----
  if (effdf) {
    dfadd <- sum((abs(alpha - minval) < minval) | (abs(alpha - upper_alpha) < minval))
    ecounts <- freqnext(freq = nvec / sum(nvec), alpha = alpha) * sum(nvec)
    if (any(!which_keep)) {
      ecounts <- c(ecounts[which_keep], sum(ecounts[!which_keep]))
    }
    dfadd <- dfadd - sum(ecounts < 0.1)
  } else {
    dfadd <- 0
  }
  df_hwe <- projout$rank - ibdr + dfadd

  ## Run test ----
  p_hwe <- stats::pchisq(q = chisq_hwe, df = df_hwe, lower.tail = FALSE)

  ## return ----
  retlist <- list(
    alpha = alpha,
    chisq_hwe = chisq_hwe,
    df_hwe = df_hwe,
    p_hwe = p_hwe
  )

  retlist$p_hwe[retlist$df_hwe == 0] <- NA_real_

  return(retlist)
}
