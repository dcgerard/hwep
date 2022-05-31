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

#' Default concentration hyperparameters for the dirichlet priors
#'
#' This is used in \code{\link{rmbayes}()} and \code{\link{rmbayesgl}()}.
#'
#' @param ploidy The ploidy of the species
#'
#' @return A list of length two. Element \code{alpha} are the concentration
#'     parameters of the gamete frequencies under the null of random mating.
#'     Element \code{beta} are the conentration parameters of the genotype
#'     frequencies under the alternative of non-random mating.
#'
#' @author David Gerard
#'
#' @noRd
conc_default <- function(ploidy) {
  alpha <- rep(1, ploidy / 2 + 1)
  beta <- pmin(floor(0:ploidy / 2), floor((ploidy - 0:ploidy) / 2)) + 1
  beta <- beta * (ploidy + 1) / sum(beta)
  return(list(alpha = alpha, beta = beta))
}

#' Probability of pair counts given total numbers in each category.
#'
#' There are \code{sum(y)} total objects of \code{length(y)} categories and
#' we randomly pair them up. The counts of pairs are in \code{A}. This
#' function will calculate the probability of \code{A} given \code{y}
#' as calculated by Levene (1949).
#'
#' @param A An upper-triangular matrix. \code{A[i,j]} is the number of pairs
#'     of categories \code{i} and \end{j}.
#' @param y A vector. \code{y[i]} are the number of objects in category
#'     \code{i}.
#' @param lg A logical. Should we return the log probability \code{TRUE} or
#'     not (\code{FALSE})?
#'
#' @references
#' \itemize{
#'   \item{H. Levene. On a matching problem arising in genetics. The Annals of Mathematical Statistics, 20 (1):91–94, 1949. doi: 10.1214/aoms/1177730093.}
#' }
#'
#' @author David Gerard
#'
#' @examples
#' A <- matrix(c(1, 2, 1,
#'               0, 3, 1,
#'               0, 0, 4),
#'             ncol = 3,
#'             byrow = TRUE)
#' y <- c(5, 9, 10)
#'
#' @noRd
dpairs <- function(A, y, lg = FALSE) {
  stopifnot(sum(y) == 2 * sum(A))
  stopifnot(length(y) == nrow(A), length(y) == ncol(A))
  stopifnot(A[lower.tri(A)] == 0)

  ## check for compatibility
  for (i in 1:length(y)) {
    if (y[[i]] != sum(A[i,]) + sum(A[, i])) {
      if (lg) {
        return(-Inf)
      } else {
        return(0)
      }
    }
  }

  ## Use formula from Levene (1949).
  n <- sum(y) / 2
  retval <- lfactorial(n) +
    sum(lfactorial(y)) -
    lfactorial(2 * n) -
    sum(lfactorial(A[upper.tri(A, diag = TRUE)])) +
    sum(A[upper.tri(A)]) * log(2)

  if (!lg) {
    retval <- exp(retval)
  }

  return(retval)
}

#' Get gamete counts from pairs of gametes matrix
#'
#' @param A An upper-triangular matrix. \code{A[i,j]} is the number of pairs
#'     of categories \code{i} and \end{j}.
#'
#' @return A vector. Element \code{i} is the number of gametes with dosage
#'     \code{i - 1}.
#'
#' @author David Gerard
#'
#' @noRd
gam_from_pairs <- function(A) {
  stopifnot(nrow(A) == ncol(A))
  stopifnot(A[lower.tri(A)] == 0)

  ngam <- nrow(A)

  y <- rep(NA_real_, length.out = ngam)
  for (i in seq_len(ngam)) {
    y[[i]] <- sum(A[i,]) + sum(A[, i])
  }

  return(y)
}

#' Marginal likelihood for tetraploids under random mating
#'
#' @param x A vector of length 5, containing the number of individuals at
#'     each dosage.
#' @param alpha A vector of length 3, containing the concentration
#'     hyperparameters for the gamete frequencies.
#'
#' @author David Gerard
#'
#' @noRd
tetra_rm_marg <- function(x, alpha, lg = FALSE) {
  stopifnot(length(x) == 5)
  stopifnot(length(alpha) == 3)
  stopifnot(x >= 0)
  stopifnot(alpha > 0)

  A <- matrix(0, nrow = 3, ncol = 3)
  A[1, 1] <- x[[1]]
  A[1, 2] <- x[[2]]
  A[2, 3] <- x[[4]]
  A[3, 3] <- x[[5]]

  mx <- -Inf
  for (i in 0:x[[3]]) {
    A[2, 2] <- i
    A[1, 3] <- x[[3]] - i
    y <- gam_from_pairs(A)
    cval <- dpairs(A = A, y = y, lg = TRUE) + ddirmult(x = y, alpha = alpha, lg = TRUE)
    mx <- log_sum_exp_2_cpp(mx, cval)
  }


  if (!lg) {
    mx <- exp(mx)
  }

  return(mx)
}

#' Marginal likelihood for hexaploids under random mating
#'
#' @param x A vector of length 7, containing the number of individuals at
#'     each dosage.
#' @param alpha A vector of length 4, containing the concentration
#'     hyperparameters for the gamete frequencies.
#'
#' @author David Gerard
#'
#' @noRd
hexa_rm_marg <- function(x, alpha, lg = FALSE) {
  stopifnot(length(x) == 7)
  stopifnot(length(alpha) == 4)
  stopifnot(x >= 0)
  stopifnot(alpha > 0)

  A <- matrix(0, nrow = 4, ncol = 4)
  A[1, 1] <- x[[1]]
  A[1, 2] <- x[[2]]
  A[3, 4] <- x[[6]]
  A[4, 4] <- x[[7]]


  mx <- -Inf
  for (i in 0:x[[3]]) {
    A[2, 2] <- i
    A[1, 3] <- x[[3]] - i
    for (j in 0:x[[4]]) {
      A[2, 3] <- j
      A[1, 4] <- x[[4]] - j
      for(k in 0:x[[5]]) {
        A[3, 3] <- k
        A[2, 4] <- x[[5]] - k
        y <- gam_from_pairs(A)
        cval <- dpairs(A = A, y = y, lg = TRUE) + ddirmult(x = y, alpha = alpha, lg = TRUE)
        mx <- log_sum_exp_2_cpp(mx, cval)
      }
    }
  }


  if (!lg) {
    mx <- exp(mx)
  }

  return(mx)
}

#' Bayes test for random mating with known genotypes
#'
#' @param nvec A vector containing the observed genotype counts,
#'     where \code{nvec[[i]]} is the number of individuals with genotype
#'     \code{i-1}. This should be of length \code{ploidy+1}.
#' @param lg A logical. Should we return the log Bayes factor (\code{TRUE})
#'     or the Bayes factor (\code{FALSE})?
#' @param alpha The concentration hyperparameters of the gamete frequencies
#'     under the null of random mating. Should be length ploidy/2 + 1.
#' @param beta The concentration hyperparameters of the genotype frequencies
#'     under the alternative of no random mating. Should be length ploidy + 1.
#' @param nburn The number of iterations in the Gibbs sampler to burn-in.
#' @param niter The number of sampling iterations in the Gibbs sampler.
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
#' rmbayes(nvec = nvec)
#'
#' nvec <- c(stats::rmultinom(n = 1, size = 1000, prob = q))
#' rmbayes(nvec = nvec)
#'
#' nvec <- c(stats::rmultinom(n = 1, size = 10000, prob = q))
#' rmbayes(nvec = nvec)
#'
#' ## Simulate under the alternative
#' q <- stats::runif(ploidy + 1)
#' q <- q / sum(q)
#'
#' ## See BF decrease
#' nvec <- c(stats::rmultinom(n = 1, size = 100, prob = q))
#' rmbayes(nvec = nvec)
#'
#' nvec <- c(stats::rmultinom(n = 1, size = 1000, prob = q))
#' rmbayes(nvec = nvec)
#'
#' nvec <- c(stats::rmultinom(n = 1, size = 10000, prob = q))
#' rmbayes(nvec = nvec)
#'
#' @author David Gerard
#'
#' @export
rmbayes <- function(nvec,
                    lg = TRUE,
                    alpha = NULL,
                    beta = NULL,
                    nburn = 1000,
                    niter = 10000) {
  ploidy <- length(nvec) - 1

  ## Default concentration parameters ----
  if (xor(is.null(alpha), is.null(beta))) {
    warning(paste0("You need to specify both alpha and beta to use custom hyperparameters.\n",
                   "If you only specify one then default values are used."))
  }

  if (is.null(alpha) || is.null(beta)) {
    clist = conc_default(ploidy = ploidy)
    alpha <- clist$alpha
    beta <- clist$beta
  } else {
    stopifnot(length(alpha) == ploidy / 2 + 1)
    stopifnot(length(beta) == ploidy + 1)
    stopifnot(alpha > 0)
    stopifnot(beta > 0)
  }

  ## Get marginal likelihoods under null and alternative
  if (ploidy == 4) {
    mnull <- tetra_rm_marg(x = nvec, alpha = alpha, lg = TRUE)
  } else if ((ploidy == 6) && sum(nvec) <= 50) {
    mnull <- hexa_rm_marg(x = nvec, alpha = alpha, lg = TRUE)
  } else {
    mnull <- gibbs_known(x = nvec, alpha = alpha, more = FALSE, lg = TRUE, B = niter, T = nburn)$mx
  }
  malt <- ddirmult(x = nvec, alpha = beta, lg = TRUE)

  ## Bayes factor ----
  lbf = mnull - malt

  if (!lg) {
    lbf <- exp(lbf)
  }

  return(lbf)
}


#' Marginal likelihood under alternative when using genotype likelihoods
#'
#' This is the R implementation of plq
#'
#' @param gl Genotype log-likelihoods. Rows index individuals and columns
#'     index genotypes.
#' @param beta The concentration parameters under the alternative.
#' @param lg A logical. Should we return the log (TRUE) or not (FALSE)?
#'
#' @author David Gerard
#'
#' @seealso plq
#'
#' @noRd
gl_alt_marg <- function(gl, beta, lg = TRUE) {
  lqvec <- log(beta / sum(beta))
  mx <- sum(apply(gl + rep(lqvec, each = nrow(gl)), 1, log_sum_exp))
  if (!lg) {
    mx <- exp(mx)
  }
  return(mx)
}

#' Bayes test for random mating using genotype log-likelihoods
#'
#' @param gl A matrix of genotype log-likelihoods. The rows index the
#'     individuals and the columns index the genotypes. So \code{gl[i,k]}
#'     is the genotype log-likelihood for individual \code{i} at
#'     dosage \code{k-1}. We assume the \emph{natural} log is used (base e).
#' @inheritParams rmbayes
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
#' \dontrun{
#' ## See BF increase
#' nvec <- c(stats::rmultinom(n = 1, size = 100, prob = q))
#' gl <- simgl(nvec = nvec)
#' rmbayesgl(gl = gl)
#'
#' nvec <- c(stats::rmultinom(n = 1, size = 1000, prob = q))
#' gl <- simgl(nvec = nvec)
#' rmbayesgl(gl = gl)
#'
#' ## Simulate under the alternative
#' q <- stats::runif(ploidy + 1)
#' q <- q / sum(q)
#'
#' ## See BF decrease
#' nvec <- c(stats::rmultinom(n = 1, size = 100, prob = q))
#' gl <- simgl(nvec = nvec)
#' rmbayesgl(gl = gl)
#'
#' nvec <- c(stats::rmultinom(n = 1, size = 1000, prob = q))
#' gl <- simgl(nvec = nvec)
#' rmbayesgl(gl = gl)
#'
#' }
#'
#' @author David Gerard
#'
#' @export
rmbayesgl <- function(gl,
                      lg = TRUE,
                      alpha = NULL,
                      beta = NULL,
                      nburn = 1000,
                      niter = 10000) {
  ploidy <- ncol(gl) - 1
  n <- nrow(gl)

  ## Default concentration parameters ----
  if (xor(is.null(alpha), is.null(beta))) {
    warning(paste0("You need to specify both alpha and beta to use custom hyperparameters.\n",
                   "If you only specify one then default values are used."))
  }

  if (is.null(alpha) || is.null(beta)) {
    clist = conc_default(ploidy = ploidy)
    alpha <- clist$alpha
    beta <- clist$beta
  } else {
    stopifnot(length(alpha) == ploidy / 2 + 1)
    stopifnot(length(beta) == ploidy + 1)
    stopifnot(alpha > 0)
    stopifnot(beta > 0)
  }

  ## Get marginal likelihoods under null and alternative
  mnull <- gibbs_gl(gl = gl, alpha = alpha, B = niter, T = nburn, more = FALSE, lg = TRUE)$mx
  malt <- gibbs_gl_alt(gl = gl, beta = beta, B = niter, T = nburn, more = FALSE, lg = TRUE)$mx

  ## Bayes factor ----
  lbf = mnull - malt

  if (!lg) {
    lbf <- exp(lbf)
  }

  return(lbf)
}

#' Stupid simulator for genotype log-likelihoods
#'
#' Gets at the idea of genotype uncertainty, but don't use this for simulations.
#' I just created this to help me debug.
#'
#' @param nvec The genotype counts. \code{nvec[k]} contains the number of
#'     individuals with genotype \code{k-1}.
#' @param sig The standard deviation. Higher means more uncertainty. You should
#'     think of this as in units of dosages.
#'
#' @author David Gerard
#'
#' @export
simgl <- function(nvec, sig = 0.1) {
  ploidy <- length(nvec) - 1
  z <- unlist(mapply(FUN = rep, x = 0:ploidy, each = nvec)) / ploidy
  zsim <- stats::rnorm(n = length(z), mean = z, sd = sig / ploidy)

  f <- function(x, y) {
    stats::dnorm(x = x, mean = y, sd = sig, log = TRUE)
  }

  gl <- outer(X = zsim, Y = (0:ploidy) / ploidy, FUN = f)

  return(gl)
}

