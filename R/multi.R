#################################
## Code to fit across multiple loci
#################################

#' HWE and random mating estimation and testing for many loci.
#'
#' Estimates and tests for either random mating or HWE across many loci
#' using \code{\link{hwetetra}()}, \code{\link{hwemom}()},
#' or \code{\link{rmlike}()}.
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
#' @param ... Any other parameters to send to \code{\link{hwetetra}()},
#'     \code{\link{hwemom}()}, or \code{\link{rmlike}()}.
#'
#' @return A data frame. The columns of which can are described in
#'     \code{\link{hwetetra}()}, \code{\link{hwemom}()},
#'     or \code{\link{rmlike}()}.
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' ## Generate random data where there is no double reduction at HWE
#' set.seed(1)
#' ploidy <- 4
#' nloc <- 100
#' size <- 1000
#' r <- runif(nloc)
#' probmat <- t(sapply(r, dbinom, x = 0:ploidy, size = ploidy))
#' nmat <- t(apply(X = probmat, MARGIN = 1, FUN = rmultinom, n = 1, size = size))
#'
#' ## Run the analysis in parallel on the local computer with two workers
#' future::plan(future::multisession, workers = 2)
#' hout <- hwefit(nmat = nmat)
#'
#' ## Shut down parallel workers
#' future::plan("sequential")
#'
hwefit <- function(nmat, type = c("hwe", "rm"), ...) {
  stopifnot(is.matrix(nmat))
  ploidy <- ncol(nmat) - 1
  stopifnot(ploidy %% 2 == 0, ploidy > 2)
  type <- match.arg(type)

  ## Choose appropriate function ----
  if (ploidy == 4) {
    fun <- hwetetra
  } else if (type == "hwe") {
    fun <- hwemom
  } else {
    fun <- rmlike
  }

  ## Register doFutures()
  oldDoPar <- registerDoFuture()
  on.exit(with(oldDoPar, foreach::setDoPar(fun=fun, data=data, info=info)), add = TRUE)

  ## Run foreach
  nvec <- NULL
  outdf <- foreach (nvec = iterators::iter(nmat, by = "row"),
                    .combine = rbind) %dopar% {
                      hout <- fun(nvec = c(nvec), ...)
                      unlist(hout)
                    }

  return(as.data.frame(outdf))
}
