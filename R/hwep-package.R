#' Hardy-Weinberg Equilibrium in Polyploids
#'
#' Inference concerning equilibrium and random mating in autopolyploids.
#' Methods are available to test for equilibrium and random mating at
#' any even ploidy level (>2) in the presence of double reduction. For
#' autopolyploid populations in equilibrium, methods are available to
#' estimate the degree of double reduction. We also provide functions
#' to calculate genotype frequencies at equilibrium, or after one or
#' several rounds of random mating, given rates of double reduction.
#' This material is based upon work supported by the
#' National Science Foundation under Grant No. 2132247. The opinions,
#' findings, and conclusions or recommendations expressed are those of
#' the author and do not necessarily reflect the views of the National
#' Science Foundation. For details of these methods, see
#' Gerard (2021) \doi{10.1101/2021.09.24.461731}.
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{hwefit}()}}{Fit either \code{\link{hwelike}()},
#'       \code{\link{rmlike}()}, \code{\link{hweustat}()}, or
#'       \code{\link{hwenodr}()} across many loci.
#'       Parallelization is supported through the \link{future} package.}
#'   \item{\code{\link{hwelike}()}}{Likelihood inference for equilibrium.
#'       This function estimates the rate of double reduction given
#'       equilibrium, and tests for at most small deviations from
#'       equilibrium.}
#'   \item{\code{\link{rmlike}()}}{Likelihood inference for random
#'       mating in polyploids. This function tests for random mating
#'       and estimates gametic frequencies given random mating.
#'       This function does not assume a model for meiosis.}
#'   \item{\code{\link{hweustat}()}}{U-statistic approach for equilibrium
#'       and double reduction. This function tests for equilibrium given
#'       double reduction rates and estimates these rates given equilibrium.}
#'   \item{\code{\link{hwenodr}()}}{Implements a likelihood ratio test
#'       that tests for equilibrium in autopolyploids given no double reduction.}
#'   \item{\code{\link{hweboot}()}}{Implements a bootstrap approach to test
#'       for equilibrium which is more appropriate for small samples and
#'       uncertain genotypes.}
#' }
#'
#' @section Other Functions:
#' \describe{
#'   \item{\code{\link{dgamete}()}}{Gamete dosage probability given
#'       parental dosage.}
#'   \item{\code{\link{drbounds}()}}{Upper bounds on the rates of double
#'       reduction given the complete equational segregation model.}
#'   \item{\code{\link{freqnext}()}}{Update genotype frequencies after one
#'       generation of random mating.}
#'   \item{\code{\link{gsegmat}()}}{Gamete dosage probabilities for all
#'       possible parental dosages.}
#'   \item{\code{\link{hwefreq}()}}{Generate equilibrium
#'       genotype frequencies.}
#'   \item{\code{\link{p_from_alpha}()}}{Obtain gamete frequencies from
#'       the major allele frequency and double reduction rates.}
#'   \item{\code{\link{zsegarray}()}}{All zygote dosage distributions
#'       given all possible parental dosages.}
#'   \item{\code{\link{zygdist}()}}{Zygote dosage distribution given
#'       one pair of parental dosages.}
#' }
#'
#' @section Citation:
#' If you find the methods in this package useful, please run the following
#' in R for citation information: \code{citation("hwep")}
#'
#' @importFrom foreach %dopar%
#' @importFrom doRNG %dorng%
#' @import doFuture
#' @importFrom foreach foreach
#' @importFrom future plan
#' @import Rcpp
#' @import methods
#' @importFrom rstan sampling
#' @importFrom RcppParallel defaultNumThreads
#'
#' @docType package
#' @name hwep-package
#' @aliases hwep
#'
#' @useDynLib hwep, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @author David Gerard
NULL
