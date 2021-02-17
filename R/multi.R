#################################
## Code to fit across multiple loci
#################################

#' HWE and random mating estimation and testing for many loci.
#'
#' Estimates and tests for either random mating or HWE across many loci
#' using \code{\link{hwetetra}()}, \code{\link{hweustat}()},
#' \code{\link{hwemom}()}, \code{\link{rmlike}()}, or \code{\link{hwenodr}()}.
#'
#' We provide parallelization support through the \link[future]{future}
#' package.
#'
#' @param nmat A matrix of counts. The rows index the loci and the columns
#'     index the genotypes. So \code{nmat[i, j]} is the number of individuals
#'     that have genotype \code{j-1} at locus \code{i}. The ploidy is
#'     assumed to be \code{ncol(nmat)-1}.
#' @param type Should we test for random mating (\code{type = "rm"})
#'     equilibrium using a U-statistic approach (\code{type = "ustat"}),
#'     equilibrium using an ad-hoc generalized method of moments approach
#'     (\code{type = "gmm"}), or equilibrium assuming no double reduction
#'     (\code{type = "nodr"})? This is only applicable for ploidies greater
#'     than 4 (unless \code{overwrite = TRUE}).
#' @param overwrite A logical. The default is to run \code{hwetetra()} if
#'     you have tetraploids, regardless of the selection of \code{type}. Set
#'     this to \code{TRUE} to overwrite this default.
#' @param ... Any other parameters to send to \code{\link{hweustat}()},
#'     \code{\link{hwetetra}()}, \code{\link{hwemom}()},
#'     \code{\link{rmlike}()}, or \code{\link{hwenodr}()}.
#'
#' @return A data frame. The columns of which can are described in
#'     \code{\link{hwetetra}()}, \code{\link{hweustat}()},
#'     \code{\link{hwemom}()}, \code{\link{rmlike}()},
#'     or \code{\link{hwenodr}()}.
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' ## Generate random data
#' set.seed(5)
#' ploidy <- 6
#' nloc <- 1000
#' size <- 1000
#' r <- 0.1
#' alpha <- 0
#' qvec <- hwefreq(r = r, alpha = alpha, ploidy = ploidy)
#' nmat <- t(rmultinom(n = nloc, size = size, prob = qvec))
#'
#' ## Run the analysis in parallel on the local computer with two workers
#' future::plan(future::multisession, workers = 6)
#' hout <- hwefit(nmat = nmat, type = "ustat", thresh_mult = 100, thresh_tot = 10)
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
#' obs <- sort(hout$chisq_hwe)
#' plot(x = qchisq(ppoints(n = length(obs)), df = 2),
#'      y = obs,
#'      xlab = "theoretical",
#'      ylab = "observed",
#'      main = "qqplot")
#' abline(0, 1, col = 2, lty = 2)
#' hist(obs, breaks = 30)
#'
#' ## Consistent estimate for alpha
#' mean(hout$alpha)
#'
hwefit <- function(nmat, type = c("ustat", "rm", "nodr", "gmm"), overwrite = FALSE, ...) {
  stopifnot(is.matrix(nmat))
  ploidy <- ncol(nmat) - 1
  stopifnot(ploidy %% 2 == 0, ploidy > 2)
  stopifnot(length(overwrite) == 1, is.logical(overwrite))
  type <- match.arg(type)

  ## Choose appropriate function ----
  if (ploidy == 4 & !overwrite) {
    fun <- hwetetra
  } else if (type == "ustat") {
    fun <- hweustat
  } else if (type == "gmm") {
    fun <- hwemom
  } else if (type == "rm") {
    fun <- rmlike
  } else if (type == "nodr") {
    fun <- hwenodr
  }

  ## Register doFutures()
  oldDoPar <- registerDoFuture()
  on.exit(with(oldDoPar, foreach::setDoPar(fun=fun, data=data, info=info)), add = TRUE)

  ## Run foreach
  nvec <- NULL
  outdf <- foreach(nvec = iterators::iter(nmat, by = "row"),
                   .combine = rbind) %dopar% {
                     hout <- fun(nvec = c(nvec), ...)
                     unlist(hout)
                     }

  return(as.data.frame(outdf))
}
