####################################
## Likelihood inference for tetraploids
####################################

#' Gamete frequencies from double reduction and allele frequency for tetraploids
#'
#' @param a The rate of double reduction
#' @param r The major allele frequency
#'
#' @return A vector of length 3. \code{p[[i]]} is the probability a gamete has
#' dosage \code{i-1}
#'
#' @author David Gerard
#'
#' @noRd
gam_from_a <- function(a, r) {
  stopifnot(length(a) == 1, length(r) == 1)
  stopifnot(r >= 0, r <= 1)
  c(p0 = (1 - 2 * r * (1 - a) / (2 + a)) * (1 - r),
    p1 = 4 * (1 - a) * r * (1 - r) / (2 + a),
    p2 = (3 * a / (2 + a) + 2 * r * (1 - a) / (2 + a)) * r)
}

#' Double reduction and allele frequency from gamete frequencies
#'
#' @param p The gamete frequencies.  \code{p[[i]]} is the probability a gamete
#'     has dosage \code{i-1}
#'
#' @return A vector of length 2. The first element is the double reduction
#' rate, the second element is the allele frequency.
#'
#' @author David Gerard
#'
#' @noRd
a_from_gam <- function(p) {
  stopifnot(length(p) == 3)
  stopifnot(abs(sum(p) - 1) < 10^-6)
  c(
    a = ((1 - 2 * p[[2]]) - (p[[3]] - p[[1]]) ^ 2) /
      ((1 + p[[2]]) - (p[[3]] - p[[1]]) ^ 2),
    r = (p[[2]] + 2 * p[[3]]) / 2
   )
}

#' Log-likelihood of tetraploid model under double reduction and allele freq
#'
#' @param a The double reduction parameter.
#' @param r The allele frequency.
#' @inheritParams hwetetra
#'
#' @noRd
lltetra <- function(a, r, nvec) {
  stopifnot(length(a) == 1, length(r) == 1)
  stopifnot(length(nvec) == 5)

  p <- gam_from_a(a = a, r = r)
  q <- stats::convolve(p, rev(p), type = "open")

  stats::dmultinom(x = nvec, prob = q, log = TRUE)
}

#' HWE likelihood inference for tetraploids
#'
#' @param nvec A vector containing the observed genotype counts,
#'     where \code{nvec[[i]]} is the number of individuals with genotype
#'     \code{i-1}. This should be of length 5 (since we have tetraploids).
#' @param upperdr The upper bound on the double reduction rate. Defaults
#'     to 1/6, which is the maximum value under the complete equational
#'     segregation model (Mather, 1935).
#' @param addval The penalization applied to each genotype for the random
#'     mating hypothesis. This corresponds to the number of counts each
#'     genotype has \emph{a priori}. Defaults to \code{1 / 100}.
#'
#' @return A list with the following elements
#' \describe{
#'  \item{\code{alpha}}{The estimated double reduction rate.}
#'  \item{\code{r}}{The estimated allele frequency.}
#'  \item{\code{q_u}}{The genotype frequencies under the unrestricted model.}
#'  \item{\code{ll_u}}{The log-likelihood under the urestricted model.}
#'  \item{\code{q_rm}}{The genotype frequencies under the random mating model.}
#'  \item{\code{ll_rm}}{The log-likelihood under the random mating model.}
#'  \item{\code{q_hwe}}{The genotype frequencies under the small deviations
#'      from HWE model.}
#'  \item{\code{ll_hwe}}{The log-likelihood under the small deviations from HWE
#'      model.}
#'  \item{\code{q_ndr}}{The genotype frequencies under the no double reduction
#'      at HWE model.}
#'  \item{\code{ll_ndr}}{The log-likelihood under the no double reduction at
#'      HWE model.}
#'  \item{\code{chisq_rm}}{The chi-square statistic against the null of random
#'      mating, with an alternative of unrestricted.}
#'  \item{\code{df_rm}}{The degrees of freedom corresponding to \code{chisq_rm}.}
#'  \item{\code{p_rm}}{The p-value corresponding to \code{chisq_rm}.}
#'  \item{\code{chisq_hwe}}{The chi-square statistic against the null of only
#'      small deviations from HWE, with an alternative of arbitrary
#'      genotype frequencies.}
#'  \item{\code{df_hwe}}{The degrees of freedom corresponding to \code{chisq_hwe}.}
#'  \item{\code{p_hwe}}{The p-value corresponding to \code{chisq_hwe}.}
#'  \item{\code{chisq_ndr}}{The chi-square statistic against the null of no
#'      double reduction at HWE, with an alternative of only small deviations
#'      from HWE.}
#'  \item{\code{df_ndr}}{The degrees of freedom corresponding to \code{chisq_ndr}.}
#'  \item{\code{p_ndr}}{The p-value corresponding to \code{chisq_ndr}.}
#' }
#' @author David Gerard
#'
#' @references
#' \itemize{
#'   \item{Mather, K. (1935). Reductional and equational separation of the chromosomes in bivalents and multivalents. Journal of Genetics, 30(1), 53-78. \doi{10.1007/BF02982205}}
#' }
#'
#' @export
#'
#' @examples
#' ## Generate data under HWE
#' set.seed(1)
#' alpha <- 0.1
#' r <- 0.4
#' q <- hwefreq(r = r, alpha = alpha, ploidy = 4)
#' nvec <- c(stats::rmultinom(n = 1, size = 500, prob = q))
#' hwetetra(nvec = nvec)
#'
#'
hwetetra <- function(nvec, upperdr = 1/6, addval = 1 / 100) {
  stopifnot(length(nvec) == 5)
  stopifnot(nvec >= 0)
  stopifnot(is.vector(nvec))
  stopifnot(addval >= 0, length(addval) == 1)

  ## Unconstrained ---------------------------
  q_u <- nvec / sum(nvec)
  ll_u <- stats::dmultinom(x = nvec, prob = q_u, log = TRUE)

  ## Random mating ---------------------------
  p_rm <- rmem(nvec = nvec, addval = addval)
  q_rm <- stats::convolve(p_rm, rev(p_rm), type = "open")
  ll_rm <- stats::dmultinom(x = nvec, prob = q_rm, log = TRUE)

  ## Small deviation from HWE ----------------
  r <- sum(0:4 * nvec) / (sum(nvec) * 4)
  oout <- stats::optim(par = 0.1,
                       fn = lltetra,
                       method = "Brent",
                       lower = 0,
                       upper = upperdr,
                       nvec = nvec,
                       r = r,
                       control = list(fnscale = -1))
  alpha <- oout$par
  p_hwe <- gam_from_a(a = alpha, r = r)
  names(p_hwe) <- NULL
  q_hwe <- stats::convolve(p_hwe, rev(p_hwe), type = "open")
  ll_hwe <- oout$value

  ## No dr ------------------------------------
  p_ndr <- c((1 - r) ^ 2, 2 * r * (1-r), r ^ 2)
  q_ndr <- stats::convolve(p_ndr, rev(p_ndr), type = "open")
  ll_ndr <- stats::dmultinom(x = nvec, prob = q_ndr, log = TRUE)

  ## q names
  names(q_u) <- 0:4
  names(q_rm) <- 0:4
  names(q_hwe) <- 0:4
  names(q_ndr) <- 0:4

  ## Chi square tests -------------------------
  chisq_rm <- -2 * (ll_rm - ll_u)
  pval_rm <- stats::pchisq(q = chisq_rm, df = 2, lower.tail = FALSE)

  chisq_hwe <- -2 * (ll_hwe - ll_u)
  pval_hwe <- stats::pchisq(q = chisq_hwe, df = 2, lower.tail = FALSE)

  chisq_ndr <- -2 * (ll_ndr - ll_hwe)
  pval_ndr <- stats::pchisq(q = chisq_ndr, df = 1, lower.tail = FALSE)

  ## Return everything -----------------------
  retlist <- list(
    alpha = alpha,
    r = r,
    q_u = q_u,
    ll_u = ll_u,
    q_rm = q_rm,
    ll_rm = ll_rm,
    q_hwe = q_hwe,
    ll_hwe = ll_hwe,
    q_ndr = q_ndr,
    ll_ndr = ll_ndr,
    chisq_rm = chisq_rm,
    df_rm = 2,
    p_rm = pval_rm,
    chisq_hwe = chisq_hwe,
    df_hwe = 2,
    p_hwe = pval_hwe,
    chisq_ndr = chisq_ndr,
    df_ndr = 1,
    p_ndr = pval_ndr)
  return(retlist)
}
