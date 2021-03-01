#' Hardy-Weinberg Equilibrium in Polyploids
#'
#' Inference concerning Hardy-Weinberg equilibrium (HWE) in
#' polyploids. Methods are available to test for HWE at any
#' ploidy level in the presence of double reduction. For
#' autopolyploid populations in HWE, methods are available
#' to estimate the degree of double reduction. We also
#' provide functions to calculate genotype frequencies at
#' equilibrium given rates of double reduction.
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{hwefit}()}}{Fit either \code{\link{hwelike}()},
#'       \code{\link{rmlike}()}, or \code{\link{hweustat}()} across
#'       many loci. Parallelization is supported through the
#'       \link{future} package.}
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
#'       that tests for HWE in autopolyploids given no double reduction.}
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
#'   \item{\code{\link{hwefreq}()}}{Generate Hardy-Weinberg equilibrium
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
#' @import doFuture
#' @importFrom foreach foreach
#'
#' @docType package
#' @name hwep-package
#' @aliases hwep
#'
#' @author David Gerard
NULL
