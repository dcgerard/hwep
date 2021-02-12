#################
## Methods when there is no double reduction
#################


#' Test for HWE in autopolyploids under the assumption of no double reduction
#'
#' We run a likelihood ratio test against the null of no HWE.
#'
#' @param nvec A vector containing the observed genotype counts,
#'     where \code{nvec[[i]]} is the number of individuals with genotype
#'     \code{i-1}. This should be of length \code{ploidy+1}.
#'
#' @return A list with some or all of the following elements
#' \describe{
#'   \item{\code{r}}{The estimated allele frequency.}
#'   \item{\code{chisq_hwe}}{The chi-square statistic against the null of
#'       equilibrium given no double reduction.}
#'   \item{\code{df_hwe}}{The degrees of freedom associated with
#'       \code{chisq_hwe}.}
#'   \item{\code{p_hwe}}{The p-value against the null of
#'       equilibrium given no double reduction.}
#' }
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' set.seed(10)
#' qvec <- c(0.2, 0.3, 0.4, 0.1)
#' nvec <- c(stats::rmultinom(n = 1, size = 100, prob = qvec))
#' hwenodr(nvec = nvec)
#'
hwenodr <- function(nvec) {
  ## Check parameters ----
  ploidy <- length(nvec) - 1
  stopifnot(ploidy > 1)
  stopifnot(nvec >= 0)
  stopifnot(is.vector(nvec))

  ## Estimate parameters ----
  n <- sum(nvec)
  qhat <- nvec / n
  rhat <- sum(0:ploidy * nvec) / (n * ploidy)
  qnull <- stats::dbinom(x = 0:ploidy, size = ploidy, prob = rhat)

  ## Run tests ----
  llnull <- stats::dmultinom(x = nvec, prob = qnull, log = TRUE)
  llalt <- stats::dmultinom(x = nvec, prob = qhat, log = TRUE)
  chisq_hwe <- -2 * (llnull - llalt)
  df_hwe <- ploidy - 1
  p_hwe <- stats::pchisq(q = chisq_hwe, df = df_hwe, lower.tail = FALSE)

  retlist <- list(chisq_hwe = chisq_hwe,
                  df_hwe = df_hwe,
                  p_hwe = p_hwe)

  return(retlist)
}
