#################
## Bayes tests for random mating
#################

#' PMF of Dirichlet-multinomial distribution
#'
#' @param x The vector of counts.
#' @param alpha The vector of concentration parameters.
#' @param lg A logical. Should we log the density (\code{TRUE}) or not
#'     (\code{FALSE})?
#'
#' @author David Gerard
#'
#' @examples
#' ddirmult(c(1, 2, 3), c(1, 1, 1))
#' ddirmult(c(2, 2, 2), c(1, 1, 1))
#'
#' @export
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

#' Default values from beta given alpha was specified.
#'
#' This makes the Bayes factor testing for random mating equal to 1 for
#' samples of size 1.
#'
#' @author David Gerard
#'
#' @noRd
beta_from_alpha <- function(alpha) {
  ploidy <- 2 * (length(alpha) - 1)
  beta <- rep(0, length.out = ploidy + 1)

  for (i in 0:ploidy) {
    il <- max(0, i - ploidy / 2)
    ih <- floor(i / 2)
    for (j in il:ih) {
      y <- rep(0, ploidy / 2 + 1)
      y[[j + 1]] <- 1
      y[[i - j + 1]] <- y[[i - j + 1]] + 1
      beta[[i + 1]] <- beta[[i + 1]] + ddirmult(x = y, alpha = alpha, lg = FALSE)
    }
  }

  beta <- beta * (ploidy + 1)

  return(beta)
}

#' Default concentration hyperparameters for the dirichlet priors
#'
#' This is used in \code{\link{rmbayes}()} and \code{\link{rmbayesgl}()}.
#'
#' @param ploidy The ploidy of the species
#'
#' @return A list of length two. Element \code{alpha} are the concentration
#'     parameters of the gamete frequencies under the null of random mating.
#'     Element \code{beta} are the concentration parameters of the genotype
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

#' Default concentration hyperparameters for the dirichlet priors for allopolyploids.
#'
#' Thi sis used in \code{\link{rmbayes}()} and \code{\link{rmbayesgl}()}.
#'
#' @param ploidy The ploidy of the species
#'
#' @return A list of length two. Element \code{alpha} are the concentration
#'     parameters of the gamete frequencies under the null of random mating.
#'     Element \code{beta} are the concentration parameters of the genotype
#'     frequencies under the alternative of non-random mating.
#'
#' @author David Gerard
#'
#' @noRd
conc_allo <- function(ploidy) {
  alpha <- exp(lchoose(ploidy / 2, 0:(ploidy / 2)) + log(ploidy + 2) - (ploidy / 2 + 1) * log(2))
  beta <- beta_from_alpha(alpha = alpha)
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
#' @param type If \code{alpha} is \code{NULL}, then the default priors depend
#'     on if you have autopolyploids (\code{"auto"}) or allopolyploids
#'     (\code{"allo"}).
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
#' @references
#' \itemize{
#'   \item{Gerard D (2022). "Bayesian tests for random mating in autopolyploids." \emph{bioRxiv}. \doi{10.1101/2022.08.11.503635}.}
#' }
#'
#' @export
rmbayes <- function(nvec,
                    lg = TRUE,
                    alpha = NULL,
                    beta = NULL,
                    nburn = 10000,
                    niter = 10000,
                    type = c("auto", "allo")) {
  ploidy <- length(nvec) - 1
  nvec <- round(nvec)
  type <- match.arg(type)

  ## Default concentration parameters ----
  if (is.null(alpha) && !is.null(beta)) {
    warning(paste0("You cannot specify beta and not specify alpha",
                   "Default values are being used."))
  }

  if (is.null(beta) && !is.null(alpha)) {
    beta <- beta_from_alpha(alpha = alpha)
  } else if (is.null(alpha) && type == "auto") {
    clist <- conc_default(ploidy = ploidy)
    alpha <- clist$alpha
    beta <- clist$beta
  } else if (is.null(alpha) && type == "allo") {
    clist <- conc_allo(ploidy = ploidy)
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
#' @param chains Number of MCMC chains. Almost always 1 is enough, but we
#'     use 2 as a default to be conservative.
#' @param cores Number of cores to use.
#' @param iter Total number of iterations.
#' @param warmup Number of those iterations used in the burnin.
#' @param method Should we use Stan (\code{"stan"}) or Gibbs sampling
#'     (\code{"gibbs"})?
#' @param ... Control arguments passed to \code{\link[rstan]{sampling}()}.
#' @inheritParams rmbayes
#'
#' @examples
#'
#' \dontrun{
#' set.seed(1)
#' ploidy <- 4
#'
#' ## Simulate under the null ----
#' p <- stats::runif(ploidy / 2 + 1)
#' p <- p / sum(p)
#' q <- stats::convolve(p, rev(p), type = "open")
#'
#' ## See BF increases
#' nvec <- c(stats::rmultinom(n = 1, size = 10, prob = q))
#' gl <- simgl(nvec = nvec)
#' rmbayesgl(gl = gl)
#' rmbayesgl(gl = gl, method = "gibbs")
#'
#' nvec <- c(stats::rmultinom(n = 1, size = 100, prob = q))
#' gl <- simgl(nvec = nvec)
#' rmbayesgl(gl = gl)
#' rmbayesgl(gl = gl, method = "gibbs")
#'
#' nvec <- c(stats::rmultinom(n = 1, size = 1000, prob = q))
#' gl <- simgl(nvec = nvec)
#' rmbayesgl(gl = gl)
#' rmbayesgl(gl = gl, method = "gibbs")
#'
#' ## Simulate under the alternative ----
#' q <- stats::runif(ploidy + 1)
#' q <- q / sum(q)
#'
#' ## See BF decreases
#' nvec <- c(stats::rmultinom(n = 1, size = 10, prob = q))
#' gl <- simgl(nvec = nvec)
#' rmbayesgl(gl = gl)
#' rmbayesgl(gl = gl, method = "gibbs")
#'
#' nvec <- c(stats::rmultinom(n = 1, size = 100, prob = q))
#' gl <- simgl(nvec = nvec)
#' rmbayesgl(gl = gl)
#' rmbayesgl(gl = gl, method = "gibbs")
#'
#' nvec <- c(stats::rmultinom(n = 1, size = 1000, prob = q))
#' gl <- simgl(nvec = nvec)
#' rmbayesgl(gl = gl)
#' rmbayesgl(gl = gl, method = "gibbs")
#'
#' }
#'
#' @author David Gerard
#'
#' @references
#' \itemize{
#'   \item{Gerard D (2022). "Bayesian tests for random mating in autopolyploids." \emph{bioRxiv}. \doi{10.1101/2022.08.11.503635}.}
#' }
#'
#' @export
rmbayesgl <- function(gl,
                      method = c("stan", "gibbs"),
                      lg = TRUE,
                      alpha = NULL,
                      beta = NULL,
                      type = c("auto", "allo"),
                      chains = 2,
                      cores = 1,
                      iter = 2000,
                      warmup = floor(iter / 2),
                      ...) {
  method <- match.arg(method)
  type <- match.arg(type)

  ## Remove rows with missing data ----
  which_row_na <- apply(is.na(gl), 1, any)
  gl <- gl[!which_row_na, , drop = FALSE]

  ## dimensions
  ploidy <- ncol(gl) - 1
  n <- nrow(gl)

  ## parallel processing
  oldpar <- options(mc.cores = cores)
  on.exit(options(oldpar))

  ## Default concentration parameters ----
  if (is.null(alpha) && !is.null(beta)) {
    warning(paste0("You cannot specify beta and not specify alpha",
                   "Default values are being used."))
  }

  if (is.null(beta) && !is.null(alpha)) {
    beta <- beta_from_alpha(alpha = alpha)
  } else if (is.null(alpha) && type == "auto") {
    clist <- conc_default(ploidy = ploidy)
    alpha <- clist$alpha
    beta <- clist$beta
  } else if (is.null(alpha) && type == "allo") {
    clist <- conc_allo(ploidy = ploidy)
    alpha <- clist$alpha
    beta <- clist$beta
  } else {
    stopifnot(length(alpha) == ploidy / 2 + 1)
    stopifnot(length(beta) == ploidy + 1)
    stopifnot(alpha > 0)
    stopifnot(beta > 0)
  }

  if (method == "stan") {
    ## Use stan to get Bayes factors.
    dat_alt <- list(gl = gl,
                    K = ploidy,
                    N = n,
                    beta = beta)
    dat_null <- list(gl = gl,
                     K = ploidy,
                     N = n,
                     alpha = alpha,
                     khalf = ploidy / 2 + 1)
    rmout_alt <- rstan::sampling(object = stanmodels$gl_alt,
                                 data = dat_alt,
                                 verbose = FALSE,
                                 show_messages = FALSE,
                                 chains = chains,
                                 iter = iter,
                                 warmup = warmup,
                                 ...)
    rmout_null <- rstan::sampling(object = stanmodels$gl_null,
                                  data = dat_null,
                                  verbose = FALSE,
                                  show_messages = FALSE,
                                  chains = chains,
                                  iter = iter,
                                  warmup = warmup,
                                  ...)
    balt <- bridgesampling::bridge_sampler(rmout_alt, verbose = FALSE, silent = TRUE)
    bnull <- bridgesampling::bridge_sampler(rmout_null, verbose = FALSE, silent = TRUE)
    lbf <- bridgesampling::bayes_factor(x1 = bnull, x2 = balt, log = TRUE)[[1]]
  } else if (method == "gibbs") {
    mnull <- gibbs_gl(gl = gl,
                      alpha = alpha,
                      B = iter - warmup,
                      T = warmup,
                      more = FALSE,
                      lg = TRUE)$mx
    malt <- gibbs_gl_alt(gl = gl,
                         beta = beta,
                         B = iter - warmup,
                         T = warmup,
                         more = FALSE,
                         lg = TRUE)$mx
    lbf = mnull - malt
  }

  if (!lg) {
    lbf <- exp(lbf)
  }

  return(lbf)
}

#' Simulator for genotype likelihoods.
#'
#' Uses the {updog} R package for simulating read counts and generating
#' genotype log-likelihoods.
#'
#' @param nvec The genotype counts. \code{nvec[k]} contains the number of
#'     individuals with genotype \code{k-1}.
#' @param rdepth The read depth. Lower means more uncertain.
#' @param od The overdispersion parameter. Higher means more uncertain.
#' @param bias The allele bias parameter. Further from 1 means more bias.
#'     Must greater than 0.
#' @param seq The sequencing error rate. Higher means more uncertain.
#' @param ret The return type. Should we just return the genotype
#'     likelihoods (\code{"gl"}), just the genotype posteriors
#'     (\code{"gp"}), or the entire updog output (\code{"all"})
#' @param est A logical. Estimate the updog likelihood parameters while
#'     genotype (\code{TRUE}) or fix them at the true values (\code{FALSE})?
#'     More realistic simulations would set this to \code{TRUE}, but it makes
#'     the method much slower.
#' @param ... Additional arguments to pass to
#'     \code{\link[updog]{flexdog_full}()}.
#'
#' @return By default, a matrix. The genotype (natural) log likelihoods.
#'     The rows index the individuals and the columns index the dosage. So
#'     \code{gl[i,j]} is the genotype log-likelihood for individual i
#'     at dosage j - 1.
#'
#' @author David Gerard
#'
#' @examples
#' set.seed(1)
#' simgl(c(1, 2, 1, 0, 0), model = "norm", est = TRUE)
#' simgl(c(1, 2, 1, 0, 0), model = "norm", est = FALSE)
#'
#' @export
simgl <- function(nvec,
                  rdepth = 10,
                  od = 0.01,
                  bias = 1,
                  seq = 0.01,
                  ret = c("gl", "gp", "all"),
                  est = FALSE,
                  ...) {
  stopifnot(is.logical(est), length(est) == 1)
  ploidy <- length(nvec) - 1
  nind <- sum(nvec)
  ret <- match.arg(ret)
  z <- unlist(mapply(FUN = rep, x = 0:ploidy, each = nvec))

  sizevec <- rep(rdepth, length.out = nind)
  refvec <- updog::rflexdog(sizevec = sizevec,
                            geno = z,
                            ploidy = ploidy,
                            seq = seq,
                            bias = bias,
                            od = od)

  if (est) {
    fout <- updog::flexdog_full(refvec = refvec,
                                sizevec = sizevec,
                                ploidy = ploidy,
                                verbose = FALSE,
                                seq = seq,
                                bias = bias,
                                od = od,
                                ...)
  } else {
    fout <- updog::flexdog_full(refvec = refvec,
                                sizevec = sizevec,
                                ploidy = ploidy,
                                verbose = FALSE,
                                seq = seq,
                                bias = bias,
                                od = od,
                                update_bias = FALSE,
                                update_od = FALSE,
                                update_seq = FALSE,
                                ...)
  }

  if (ret == "gl") {
    retval <- fout$genologlike
  } else if (ret == "gp") {
    retval <- fout$postmat
  } else if (ret == "all") {
    retval <- fout
  }

  return(retval)
}

#' R Implementation of gibbs_gl()
#'
#' @inheritParams gibbs_gl
#'
#' @author David Gerard
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
#' nvec <- c(stats::rmultinom(n = 1, size = 100, prob = q))
#' gl <- simgl(nvec = nvec)
#' alpha <- rep(1, 5)
#' gibbs_gl_r(gl = gl, alpha = alpha)
#'
#'
#' @noRd
gibbs_gl_r <- function(gl,
                       alpha,
                       nsamp = 10000,
                       nburn = 1000,
                       more = FALSE,
                       lg = FALSE) {
  ploidy <- ncol(gl) - 1
  nind <- nrow(gl)

  pvec <- rdirichlet1(alpha)
  qvec <- stats::convolve(pvec, rev(pvec), type = "open")
  ptilde <- pvec
  qtilde <- qvec
  lp <- ddirichlet(x = ptilde, alpha = alpha, lg = TRUE) + gl_alt_marg(gl = gl, beta = qtilde, lg = TRUE)
  postval <- -Inf

  for (i in seq_len(nburn + nsamp)) {
    ## sample genotypes ----
    postmat <- gl + rep(log(qvec), each = nind)
    postmat <- exp(postmat - apply(postmat, 1, log_sum_exp))
    zvec <- apply(postmat, 1, function(x) sample(x = 0:ploidy, size = 1, replace = FALSE, prob = x))

    ## calculate genotype counts ----
    xvec <- table(factor(zvec, levels = 0:ploidy))

    ## sample y ----
    yvec <- samp_gametes(x = xvec, p = pvec)

    ## sample p ----
    pvec <- rdirichlet1(alpha = yvec + alpha)

    ## calculate q ---
    qvec <- stats::convolve(pvec, rev(pvec), type = "open")

    if (i <= nburn) {
      lp_new <- ddirichlet(x = pvec, alpha = alpha, lg = TRUE) + gl_alt_marg(gl = gl, beta = qvec, lg = TRUE)
      if (lp_new > lp) {
        lp <- lp_new
        ptilde <- pvec
        qtilde <- qvec
      }
    } else {
      postval <- log_sum_exp_2_cpp(postval, ddirichlet(x = ptilde, yvec + alpha, lg = TRUE))
    }
  }
  postval <- postval - log(nsamp)

  ## return ----
  retlist <- list()
  retlist$mx <- lp - postval

  if (!lg) {
    retlist$mx <- exp(retlist$mx)
  }

  return(retlist)
}

#' Implementation of gibbs_gl_alt() in R
#'
#' @inheritParams gibbs_gl_alt()
#'
#' @author David Gerard
#'
#' @noRd
gibbs_gl_alt_r <- function(gl, beta, nsamp = 10000, nburn = 1000, more = FALSE, lg = FALSE) {
  nind <- nrow(gl)
  ploidy <- ncol(gl) - 1
  qvec <- rdirichlet1(alpha = beta)
  qtilde <- qvec
  lp <- ddirichlet(x = qtilde, alpha = beta, lg = TRUE) + gl_alt_marg(gl = gl, beta = qtilde, lg = TRUE)
  postval <- -Inf

  for (i in seq_len(nburn + nsamp)) {
    ## sample genotypes ----
    postmat <- gl + rep(log(qvec), each = nind)
    postmat <- exp(postmat - apply(postmat, 1, log_sum_exp))
    zvec <- apply(postmat, 1, function(x) sample(x = 0:ploidy, size = 1, replace = FALSE, prob = x))

    ## calculate genotype counts ----
    xvec <- table(factor(zvec, levels = 0:ploidy))

    ## sample q ----
    qvec <- rdirichlet1(xvec + beta)

    if (i <= nburn) {
      lp_new <- ddirichlet(x = qvec, alpha = beta, lg = TRUE) + gl_alt_marg(gl = gl, beta = qvec, lg = TRUE)
      if (lp_new > lp) {
        lp <- lp_new
        qtilde <- qvec
      }
    } else {
      postval <- log_sum_exp_2_cpp(postval, ddirichlet(x = qtilde, alpha = xvec + beta, lg = TRUE))
    }
  }
  postval <- postval - log(nsamp)

  ## return ----
  retlist <- list()
  retlist$mx <- lp - postval

  if (!lg) {
    retlist$mx <- exp(retlist$mx)
  }

  return(retlist)
}

##############################################################################
## F1 / S1 test
##############################################################################

#' Bayes test for F1/S1 genotype frequencies using genotype likelihoods
#'
#' Uses \code{\link[updog]{get_q_array}()} from the updog R package to
#' calculate segregation probabilities (assuming no double reduction) and
#' tests that offspring genotypes follow this distribution.
#'
#' @inheritParams rmbayesgl
#' @param method Should we test for F1 proportions (\code{"f1"}) or
#'     S1 proportions (\code{"s1"})?
#' @param p1gl A vector of genotype log-likelihoods for parent 1.
#'     \code{p1gl[k]} is the log-likelihood of parent 1's data given
#'     their genotype is \code{k}.
#' @param p2gl A vector of genotype log-likelihoods for parent 2.
#'     \code{p2gl[k]} is the log-likelihood of parent 2's data given
#'     their genotype is \code{k}.
#'
#' @author David Gerard
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' ploidy <- 4
#'
#' ## Simulate under the null ----
#' q <- updog::get_q_array(ploidy = 4)[3, 3, ]
#'
#' ## See BF increases
#' nvec <- c(stats::rmultinom(n = 1, size = 10, prob = q))
#' gl <- simgl(nvec = nvec)
#' menbayesgl(gl = gl, method = "f1")
#'
#' nvec <- c(stats::rmultinom(n = 1, size = 100, prob = q))
#' gl <- simgl(nvec = nvec)
#' menbayesgl(gl = gl, method = "f1")
#'
#' nvec <- c(stats::rmultinom(n = 1, size = 1000, prob = q))
#' gl <- simgl(nvec = nvec)
#' menbayesgl(gl = gl, method = "f1")
#'
#' ## Simulate under the alternative ----
#' q <- stats::runif(ploidy + 1)
#' q <- q / sum(q)
#'
#' ## See BF decreases
#' nvec <- c(stats::rmultinom(n = 1, size = 10, prob = q))
#' gl <- simgl(nvec = nvec)
#' menbayesgl(gl = gl, method = "f1")
#'
#' nvec <- c(stats::rmultinom(n = 1, size = 100, prob = q))
#' gl <- simgl(nvec = nvec)
#' menbayesgl(gl = gl, method = "f1")
#'
#' nvec <- c(stats::rmultinom(n = 1, size = 1000, prob = q))
#' gl <- simgl(nvec = nvec)
#' menbayesgl(gl = gl, method = "f1")
#'
#' }
#'
#' @references
#' \itemize{
#'   \item{Gerard D (2022). "Bayesian tests for random mating in autopolyploids." \emph{bioRxiv}. \doi{10.1101/2022.08.11.503635}.}
#' }
#'
#' @export
menbayesgl <- function(gl,
                       method = c("f1", "s1"),
                       p1gl = NULL,
                       p2gl = NULL,
                       lg = TRUE,
                       beta = NULL,
                       chains = 2,
                       cores = 1,
                       iter = 2000,
                       warmup = floor(iter / 2),
                       ...) {
  ## Remove rows with missing data ----
  which_row_na <- apply(is.na(gl), 1, any)
  gl <- gl[!which_row_na, , drop = FALSE]

  ## Check input ----
  ploidy <- ncol(gl) - 1
  nind <- nrow(gl)
  method <- match.arg(method)

  if (is.null(p1gl)) {
    p1gl <- rep(0, ploidy + 1)
  }
  if (is.null(p2gl) && method == "f1") {
    p2gl <- rep(0, ploidy + 1)
  }
  if (!is.null(p2gl) && method == "s1") {
    warning("p2gl is not used because method = 's1'")
  }
  if (is.null(beta)) {
    beta <- rep(1, ploidy + 1)
  }

  ## Get alternative marginal likelihood
  dat_alt <- list(gl = gl,
                  K = ploidy,
                  N = nind,
                  beta = beta)
  rmout_alt <- rstan::sampling(object = stanmodels$gl_alt,
                               data = dat_alt,
                               verbose = FALSE,
                               show_messages = FALSE,
                               chains = chains,
                               iter = iter,
                               warmup = warmup,
                               ...)
  balt <- bridgesampling::bridge_sampler(rmout_alt, verbose = FALSE, silent = TRUE)

  ## Get null marginal likelihood
  qarray <- updog::get_q_array(ploidy = ploidy)

  if (method == "s1") {
    p1post <- p1gl - log_sum_exp_cpp(p1gl)
    bnull <- -Inf
    for (k in 0:ploidy) {
      bnull <- log_sum_exp_2_cpp(bnull, p1post[[k + 1]] + plq(gl = gl, beta = qarray[k + 1, k + 1, ], lg = TRUE))
    }
  } else if (method == "f1") {
    p1post <- p1gl - log_sum_exp(p1gl)
    p2post <- p2gl - log_sum_exp(p2gl)
    bnull <- -Inf
    for (k1 in 0:ploidy) {
      for (k2 in 0:ploidy) {
        bnull <- log_sum_exp_2_cpp(bnull, p1post[[k1 + 1]] + p2post[[k2 + 1]] + plq(gl = gl, beta = qarray[k1 + 1, k2 + 1, ], lg = TRUE))
      }
    }
  }

  lbf <- bnull - balt[[1]]

  if (!lg) {
    lbf <- exp(lbf)
  }

  return(lbf)
}
