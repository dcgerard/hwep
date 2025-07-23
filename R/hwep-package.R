#' @keywords internal
#'
#' @import doFuture
#' @import Rcpp
#' @import methods
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom doRNG %dorng%
#' @importFrom future plan
#' @importFrom rstan sampling
#' @importFrom rstantools use_rstan
#' @importFrom RcppParallel defaultNumThreads
#' @importFrom Rcpp sourceCpp
#' @useDynLib hwep, .registration = TRUE
#'
"_PACKAGE"
