########################
## Bootstrap procedure
########################


#' Bootstrap procedure to test for equilibrium
#'
#' Iteratively resample individuals/genotypes, calculating the U-statistic
#' for each resample, and use these resamples to test against the null
#' of no equilibrium.
#'
#' @param n One of two forms
#' \describe{
#'   \item{A vector of length ploidy + 1}{Element i is the number of
#'         individuals with genotype i.}
#'   \item{A matrix with nsamp rows and ploidy+1 columns}{Element (i, j)
#'         is the posterior probability that individual i has ploidy j-1.}
#' }
#' @param nboot The number of bootstrap samples to run.
#' @param more A logical. Should we return the bootstrap replicates
#'     (\code{FALSE}) or just the p-value, with 95% confidence interval
#'     of the p-value (\code{TRUE}).
#'
#' @author David Gerard
#'
#' @return A list with some or all of the following elements
#' \describe{
#'   \item{\code{p_hwe}}{The bootstrap p-value against the null of equilibrium.}
#'   \item{\code{p_ci}}{The 95% confidence interval of p_hwe.}
#'   \item{\code{alpha_boot}}{The bootstrap samples of the double reduction
#'       parameter.}
#'   \item{\code{u_boot}}{The bootstrap samples of the U-statistic.}
#' }
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' ploidy <- 6
#' size <- 100
#' r <- 0.5
#' alpha <- 0.1
#' qvec <- hwefreq(r = r, alpha = alpha, ploidy = ploidy)
#' nvec <- c(rmultinom(n = 1, size = size, prob = qvec))
#' bout <- hweboot(n = nvec, more = TRUE, nboot = 1000)
#' bout$p_hwe
#' bout$p_ci
#' hist(bout$test_boot)
#' abline(v = bout$test_stat, lty = 2, col = 2)
#'
hweboot <- function(n, nboot = 2000, more = FALSE) {
  ## Check inputs ----
  if (is.vector(n)) {
    ploidy <- length(n) - 1
    nind <- sum(n)
    type <- "geno"
  } else if (is.matrix(n)) {
    ploidy <- ncol(n) - 1
    nind <- nrow(n)
    type <- "genolike"
  } else {
    stop("hweboot: n is neither a vector nor a matrix")
  }
  stopifnot(ploidy %% 2 == 0, ploidy >= 4)
  ibdr <- floor(ploidy / 4)
  stopifnot(nboot >= 1, length(nboot) == 1)
  stopifnot(is.logical(more), length(more) == 1)

  ## Parameters for optimization ----
  minval <- sqrt(.Machine$double.eps)
  upper_alpha <- drbounds(ploidy = ploidy)
  ometh <- ifelse(ibdr == 1, "Brent", "L-BFGS-B")

  ## Estimate parameters using original data ----
  if (type == "geno") {
    nvec <- n
  } else if (type == "genolike") {
    ## use posterior modes
    ## nvec <- as.vector(table(c(apply(X = n, MARGIN = 1, FUN = which.max) - 1, 0:ploidy)) - 1)
    ## use mean counts
    nvec <- colSums(n)
    stopifnot(abs(sum(nvec) - nind) < 10^-4)
  }
  qhat <- nvec / nind

  oout <- stats::optim(par = rep(minval, ibdr),
                         fn = uobj,
                         method = ometh,
                         lower = rep(minval, ibdr),
                         upper = upper_alpha,
                         nvec = nvec,
                         omega = NULL,
                         which_keep = NULL)

  alphahat <- oout$par
  fq <- freqnext(freq = qhat, alpha = oout$par, check = FALSE)
  u_stat <- (qhat - fq) * sqrt(nind)
  test_stat <- sum(u_stat ^ 2)

  ## Run bootstrap ----
  alpha_boot <- matrix(NA_real_, nrow = nboot, ncol = ibdr)
  u_boot <- matrix(NA_real_, nrow = nboot, ncol = ploidy + 1)
  for (b in seq_len(nboot)) {
    ## Resample
    if (type == "geno") {
      n_boot <- c(stats::rmultinom(n = 1, size = nind, prob = qhat))
    } else if (type == "genolike") {
      ind_boot <- sample(x = seq_len(nind), size = nind, replace = TRUE)
      genovec <- apply(X = n[ind_boot, ], MARGIN = 1, FUN = function(x) sample(x = 0:ploidy, size = 1, prob = x))
      n_boot <- as.vector(table(c(genovec, 0:ploidy)) - 1)
    }

    ## Estimate alpha
    oout <- stats::optim(par = rep(minval, ibdr),
                         fn = uobj,
                         method = ometh,
                         lower = rep(minval, ibdr),
                         upper = upper_alpha,
                         nvec = n_boot,
                         omega = NULL,
                         which_keep = NULL)

    alpha_boot[b, ] <- oout$par

    ## Calculate u-stat
    qhat_boot <- n_boot / nind
    fq_boot <- freqnext(freq = qhat_boot, alpha = oout$par, check = FALSE)
    u_boot[b, ] <- (qhat_boot - fq_boot) * sqrt(nind)
  }

  ## Calculate test_boot
  u_mean <- colMeans(u_boot)
  diffmat <- t(t(u_boot) - u_mean)
  test_boot <- apply(X = diffmat, MARGIN = 1, FUN = function(x) sum(x ^ 2))

  ## p-values
  p_hwe <- mean(test_boot > test_stat)
  p_ci <- stats::binom.test(x = sum(test_boot > test_stat), n = nboot)$conf.int

  ## Return
  retlist <- list(alpha = alphahat,
                  p_hwe = p_hwe,
                  p_ci = p_ci)

  if (more) {
    retlist$u_stat <- u_stat
    retlist$test_stat <- test_stat
    retlist$alpha_boot <- alpha_boot
    retlist$u_boot <- u_boot
    retlist$test_boot <- test_boot
    retlist$u_mean <- u_mean
  }

  return(retlist)
}
