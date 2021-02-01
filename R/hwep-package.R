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
#' @section Functions:
#' The main functions are:
#' \describe{
#'   \item{\code{\link{hwelike}()}}{Likelihood inference for HWE in polyploids.
#'       This function tests for HWE and estimates gametic frequencies
#'       given HWE. This function does not assume a model for meiosis.}
#'   \item{\code{\link{hwemom}()}}{Generalized method of moments inference
#'       for double reduction. This function tests for equilibrium given
#'       double reduction rates and estimates these rates given equilibrium.
#'       Note that for ploidies greater than four, not all gametic
#'       frequencies are possible given this model for meiosis. That is,
#'       the test might reject not because the population is not in HWE,
#'       but because the model for meiosis is incorrect. So
#'       \code{\link{hwelike}()} is more general if you just want to
#'       test for HWE and do not care about double reduction.}
#' }
#'
#' Other functions include:
#' \describe{
#'   \item{\code{\link{dgamete}()}}{Gamete dosage probability given
#'       parental dosage.}
#'   \item{\code{\link{freqnext}()}}{Update genotype frequencies after one
#'       generation of random mating.}
#'   \item{\code{\link{hwefreq}()}}{Generate Hardy-Weinberg equilibrium
#'       genotype frequencies.}
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
#' @docType package
#' @name hwep-package
#' @aliases hwep
#'
#' @author David Gerard
NULL
