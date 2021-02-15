#################################
## Code to fit across multiple loci
#################################

#' HWE and random mating estimation and testing for many loci.
#'
#' Estimates and tests for either random mating or HWE across many loci
#' using \code{\link{hwetetra}()}, \code{\link{hwemom}()},
#' \code{\link{rmlike}()}, or \code{\link{hwenodr}()}.
#'
#' We provide parallelization support through the \link[future]{future}
#' package.
#'
#' @param nmat A matrix of counts. The rows index the loci and the columns
#'     index the genotypes. So \code{nmat[i, j]} is the number of individuals
#'     that have genotype \code{j-1} at locus \code{i}. The ploidy is
#'     assumed to be \code{ncol(nmat)-1}.
#' @param type Should we test for random mating (\code{type = "rm"})
#'     or equilibrium (\code{type = "hwe"}). For tetraploids, both
#'     tests will be run, so this is only applicable for ploidies greater
#'     than 4.
#' @param overwrite A logical. The default is to run \code{hwetetra()} if
#'     you have tetraploids, regardless of the selection of \code{type}. Set
#'     this to \code{TRUE} to overwrite this default.
#' @param ... Any other parameters to send to \code{\link{hwetetra}()},
#'     \code{\link{hwemom}()}, \code{\link{rmlike}()}, or
#'     \code{\link{hwenodr}()}.
#'
#' @return A data frame. The columns of which can are described in
#'     \code{\link{hwetetra}()}, \code{\link{hwemom}()},
#'     \code{\link{rmlike}()}, or \code{\link{hwenodr}()}.
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' ## Generate random data where there is no double reduction at HWE
#' set.seed(2)
#' ploidy <- 8
#' nloc <- 100
#' size <- 10000
#' r <- runif(nloc)
#' probmat <- t(sapply(r, dbinom, x = 0:ploidy, size = ploidy))
#' nmat <- t(apply(X = probmat, MARGIN = 1, FUN = rmultinom, n = 1, size = size))
#'
#' ## Run the analysis in parallel on the local computer with two workers
#' future::plan(future::multisession, workers = 2)
#' hout <- hwefit(nmat = nmat, obj = "pearson")
#'
#' ## Shut down parallel workers
#' future::plan("sequential")
#'
#' ## Show that test statistic follows theoretical distribution
#' plot(x = qchisq(ppoints(n = nrow(hout)), df = hout$df_hwe),
#'      y = sort(hout$chisq_hwe),
#'      xlab = "theoretical",
#'      ylab = "observed",
#'      main = "qqplot")
#' abline(0, 1, col = 2, lty = 2)
#'
hwefit <- function(nmat, type = c("hwe", "rm", "nodr"), overwrite = FALSE, ...) {
  stopifnot(is.matrix(nmat))
  ploidy <- ncol(nmat) - 1
  stopifnot(ploidy %% 2 == 0, ploidy > 2)
  stopifnot(length(overwrite) == 1, is.logical(overwrite))
  type <- match.arg(type)

  ## Choose appropriate function ----
  if (ploidy == 4 & !overwrite) {
    fun <- hwetetra
  } else if (type == "hwe") {
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
