#################################
## Code to fit across multiple loci
#################################

#' Equilibrium and random mating estimation and testing for many loci.
#'
#' Estimates and tests for either equilibrium or random mating  across many loci
#' using \code{\link{hwelike}()}, \code{\link{hweustat}()},
#' \code{\link{rmlike}()}, \code{\link{hwenodr}()}, or \code{\link{hweboot}()}.
#'
#' We provide parallelization support through the \link[future]{future}
#' package.
#'
#' @param nmat A matrix of counts. The rows index the loci and the columns
#'     index the genotypes. So \code{nmat[i, j]} is the number of individuals
#'     that have genotype \code{j-1} at locus \code{i}. The ploidy is
#'     assumed to be \code{ncol(nmat)-1}.
#' @param type The method to use:
#' \describe{
#'   \item{\code{"ustat"}}{U-statistic approach to test for equilibrium and
#'       estimate double reduction rates given equilibrium. The default.
#'       See \code{\link{hweustat}()}.}
#'   \item{\code{"mle"}}{Maximum likelihood estimation and testing. Only
#'       supported for ploidies less than or equal to 10.
#'       See \code{\link{hwelike}()}.}
#'   \item{\code{"rm"}}{Testing random mating, and estimating gamete
#'       frequencies given random mating. See \code{\link{rmlike}()}.}
#'   \item{\code{"nodr"}}{Testing equilibrium given no double reduction.
#'       See \code{\link{hwenodr}()}.}
#'   \item{\code{"boot"}}{Bootstrap approach to test for equilibrium.
#'       See \code{\link{hweboot}()}.}
#' }
#' @param effdf A logical. Should we use the effective degrees of freedom?
#'     Only applicable if \code{type = "mle"} or \code{type = "ustat"}.
#' @param thresh A non-negative numeric. The threshold for aggregating
#'     genotypes. Only applicable if \code{type = "mle"},
#'     \code{type = "ustat"}, or \code{type = "rm"}.
#' @param nboot The number of bootstrap iterations to use if
#'     \code{type = "boot"}.
#'
#' @return A data frame. The columns of which can are described in
#'     \code{\link{hwelike}()}, \code{\link{hweustat}()},
#'     \code{\link{rmlike}()}, or \code{\link{hwenodr}()}.
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' ## Generate random data
#' set.seed(5)
#' ploidy <- 4
#' nloc <- 100
#' size <- 1000
#' r <- 0.25
#' alpha <- 1/12
#' qvec <- hwefreq(r = r, alpha = alpha, ploidy = ploidy)
#' nmat <- t(rmultinom(n = nloc, size = size, prob = qvec))
#'
#' ## Run the analysis in parallel on the local computer with two workers
#' future::plan(future::multisession, workers = 2)
#' hout <- hwefit(nmat = nmat,
#'                type = "ustat")
#'
#' ## Shut down parallel workers
#' future::plan("sequential")
#'
#' ## Show that p-values are uniform
#' obs <- sort(hout$p_hwe)
#' plot(x = ppoints(n = length(obs)),
#'      y = obs,
#'      xlab = "theoretical",
#'      ylab = "observed",
#'      main = "qqplot")
#' abline(0, 1, col = 2, lty = 2)
#' mean(hout$p_hwe < 0.05, na.rm = TRUE)
#'
#' ## Consistent estimate for alpha
#' alpha
#' mean(hout$alpha1)
#'
hwefit <- function(nmat,
                   type = c("ustat", "mle", "rm", "nodr", "boot"),
                   effdf = TRUE,
                   thresh = 3,
                   nboot = 2000) {
  ## Check parameters ----
  stopifnot(is.matrix(nmat))
  ploidy <- ncol(nmat) - 1
  stopifnot(ploidy %% 2 == 0, ploidy > 2)
  stopifnot(is.logical(effdf), length(effdf) == 1)
  stopifnot(length(thresh) == 1, thresh >= 0)
  stopifnot(length(nboot) == 1, nboot >= 1)
  type <- match.arg(type)

  if (type == "mle") {
    stopifnot(ploidy <= 10)
  }

  ## Register doFutures() ----
  oldDoPar <- registerDoFuture()
  on.exit(with(oldDoPar, foreach::setDoPar(fun=fun, data=data, info=info)), add = TRUE)

  ## Run foreach ----
  nvec <- NULL
  outdf <- foreach(nvec = iterators::iter(nmat, by = "row"),
                   .combine = rbind) %dorng% {
                     if (type == "mle") {
                       hout <- hwelike(nvec = c(nvec),
                                       thresh = thresh,
                                       effdf = effdf)
                     } else if (type == "ustat") {
                       hout <- hweustat(nvec = c(nvec),
                                        thresh = thresh,
                                        effdf = effdf)
                     } else if (type == "rm") {
                       hout <- rmlike(nvec = c(nvec),
                                      thresh = thresh)
                     } else if (type == "nodr") {
                       hout <- hwenodr(nvec = c(nvec))
                     } else if (type == "boot") {
                       hout <- hweboot(n = c(nvec), nboot = nboot, more = FALSE)
                     }
                     unlist(hout)
                   }

  ## for consistency with higher ploidy levels ----
  if (ploidy %in% c(4, 6)) {
    colnames(outdf)[colnames(outdf) == "alpha"] <- "alpha1"
  }

  return(as.data.frame(outdf))
}

#' Equilibrium and random mating estimation and testing for many loci using
#' posterior genotype probabilities.
#'
#' Right now, only the bootstrap method in \code{\link{hweboot}()} is
#' implemented. We hope in the future to offer more options.
#'
#' @param genoarray A three-dimensional array of posterior genotype probabilties.
#'     The first dimension indexes
#'     the loci, the second dimension indexes the individuals, and the
#'     third dimension indexes the genotypes. So \code{genoarray[i,j,k]}
#'     contains the posterior probability that individual \code{j} at locus
#'     \code{i} has genotype \code{k-1}.
#' @param nboot The number of bootstrap samples to apply to
#'     \code{\link{hweboot}()}.
#'
#' @author David Gerard
#'
#' @noRd
hwefit_genolike <- function(genoarray, nboot = 2000) {
  stopifnot(is.array(genoarray))

}
