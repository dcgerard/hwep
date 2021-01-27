###########################
## Utility functions
###########################

#' Use log-sum-exp trick on vector
#'
#' Calculates \deqn{\log(\sum_i\exp(x_i))}
#'
#' @param x A numeric. We want to sum the exponent of x, then log.
#' @param na.rm A logical. Should we remove \code{NA}'s (\code{TRUE}) or not
#'     (\code{FALSE})? Defaults to \code{FALSE}.
#'
#' @author David Gerard
#'
#' @noRd
log_sum_exp <- function(x, na.rm = FALSE) {
  stopifnot(is.numeric(x))
  stopifnot(is.logical(na.rm))
  stopifnot(length(na.rm) == 1)

  x_max <- max(x, na.rm = na.rm)
  return(log(sum(exp(x - x_max), na.rm = na.rm)) + x_max)
}
