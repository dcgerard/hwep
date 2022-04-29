#' Get simultaneous confidence bands for a uniform QQ-plot
#'
#' This will provide 100(1-a)% simultaneous confidence bands for a
#' sample of size \code{n}. It does this by the "tail-sensitive" approach
#' of Aldor-Noiman et al (2013), which uses simulated uniform vectors. The
#' number of simulations is controlled by \code{nsamp}.
#'
#' The procedure used is described in Aldor-Noiman et al (2013). But note
#' that they have a mistake in their paper. Step (e) of their algorithm on
#' page 254 should be the CDF of the Beta distribution, not the quantile
#' function.
#'
#' @param n Sample size.
#' @param nsamp Number of simulation repetitions.
#' @param a The significance level.
#'
#' @return A list of length 3. The \code{$lower} and \code{$upper} confidence
#' limits at uniform quantiles \code{$q}.
#'
#' @references
#' \itemize{
#'   \item{Aldor-Noiman, S., Brown, L. D., Buja, A., Rolke, W., & Stine, R. A. (2013). The power to see: A new graphical test of normality. The American Statistician, 67(4), 249-260.}
#' }
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' ts <- ts_bands(100)
#'
#' graphics::plot(x = ts$q,
#'                y = ts$upper,
#'                type = "l",
#'                xlim = c(0, 1),
#'                ylim = c(0, 1),
#'                xlab = "Theoretical Quantiles",
#'                ylab = "Empirical Quantiles")
#' graphics::lines(x = ts$q, y = ts$lower)
#' graphics::lines(x = ts$q, y = ts$q, lty = 2)
#'
ts_bands <- function(n, nsamp = 1000, a = 0.05) {
  alpha <- seq_len(n)
  beta <- n + 1 - seq_len(n)
  simmat <- matrix(stats::runif(n * nsamp), ncol = nsamp)
  simmat <- apply(simmat, 2, sort)
  amat <- apply(simmat, 2, function(x) stats::pbeta(q = x, shape1 = alpha, shape2 = beta))
  cvec <- apply(amat, 2, function(x) 2 * min(c(x, 1 - x)))
  gamma <- stats::quantile(cvec, probs = a)
  upper <- stats::qbeta(p = 1 - gamma / 2, shape1 = alpha, shape2 = beta)
  lower <- stats::qbeta(p = gamma / 2, shape1 = alpha, shape2 = beta)
  return(list(q = stats::ppoints(n), lower = lower, upper = upper))
}

#' QQ-plot for p-values
#'
#' This will create a QQ-plot for p-values, comparing them to a uniform
#' distribution. We make our plot on the -log10 scale. We calculate
#' simultaneous confidence bands by the Tail Sensitive approach of
#' Aldor-Noiman et al (2013).
#'
#' @param pvals A vector of p-values.
#' @param method Should we use base plotting or ggplot2 (if installed)?
#' @param band_type Should we use the method of Aldor-Noiman et al (2013) or
#'     pointwise based on beta? Pointwise is not recommended since there is
#'     strong dependence between order statistics, and if one is beyond
#'     the pointwise bands, then likely lots are also beyond them.
#' @param conf_level Confidence level for the bands.
#'
#' @author David Gerard
#'
#' @references
#' \itemize{
#'   \item{Aldor-Noiman, S., Brown, L. D., Buja, A., Rolke, W., & Stine, R. A. (2013). The power to see: A new graphical test of normality. The American Statistician, 67(4), 249-260.}
#' }
#'
#' @seealso
#' \itemize{
#'   \item{The \code{qqPlot()} function from the car package.}
#' }
#'
#' @examples
#' set.seed(1)
#' pvals <- runif(100)
#' qqpvalue(pvals, band_type = "ts", method = "base")
#'
#' \dontrun{
#' qqpvalue(pvals, band_type = "ts", method = "ggplot2")
#' }
#'
#' @export
qqpvalue <- function(pvals,
                     method = c("ggplot2", "base"),
                     band_type = c("ts", "pointwise"),
                     conf_level = 0.95) {
  method <- match.arg(method)
  band_type <- match.arg(band_type)
  stopifnot(conf_level >= 0, conf_level <= 1)
  if (method == "ggplot2" &
      !requireNamespace(package = "ggplot2", quietly = TRUE)) {
    message("ggplot2 not installed, using base")
    method <- "base"
  }

  pvals <- sort(pvals)
  theo <- stats::ppoints(n = length(pvals))

  pvals_nl10 <- -log10(pvals)
  theo_nl10 <- -log10(theo)


  if (band_type == "ts") {
    bounds <- ts_bands(n = length(pvals), a = 1 - conf_level)
    lower <- bounds$lower
    upper <- bounds$upper
  } else if (band_type == "pointwise") {
    shape1 <- seq_len(length(pvals))
    shape2 <- length(pvals) + 1 - seq_len(length(pvals))
    lower <- stats::qbeta(p = (1 - conf_level)/2, shape1 = shape1, shape2 = shape2)
    upper <- stats::qbeta(p = 1 - (1 - conf_level)/2, shape1 = shape1, shape2 = shape2)
    bounds <- list(lower = lower, upper = upper)
  }

  maxlab <- ceiling(max(c(pvals_nl10, theo_nl10, -log10(bounds$lower))))
  breakvec <- 10^-(0:maxlab)
  lim <- c(min(c(pvals, theo, bounds$lower)), max(c(pvals, theo, bounds$upper)))

  max_limval <- max(c(pvals_nl10, theo_nl10))

  if(requireNamespace(package = "ggplot2", quietly = TRUE) &
     method == "ggplot2") {
    pl <- ggplot2::ggplot(data.frame(theo = theo_nl10,
                                     pvals = pvals_nl10,
                                     lower = -log10(lower),
                                     upper = -log10(upper))) +
      ggplot2::geom_point(mapping = ggplot2::aes(x = theo, y = pvals)) +
      ggplot2::geom_abline() +
      ggplot2::geom_line(mapping = ggplot2::aes(x = theo, y = lower), lty = 2) +
      ggplot2::geom_line(mapping = ggplot2::aes(x = theo, y = upper), lty = 2) +
      ggplot2::scale_x_continuous(name = "Theoretical Quantiles",
                                  breaks = 0:maxlab,
                                  labels = breakvec,
                                  minor_breaks = NULL) +
      ggplot2::scale_y_continuous(name = "Observed P-values",
                                  breaks = 0:maxlab,
                                  labels = breakvec,
                                  minor_breaks = NULL) +
      ggplot2::coord_cartesian(xlim = c(0, max_limval),
                               ylim = c(0, max_limval)) +
      ggplot2::theme_bw()
    print(pl)
    return(invisible(pvals))
  } else {
    oldpar <- graphics::par(pch = 16,
                            mar = c(3, 3, 0.5, 0.5),
                            mgp = c(1.8, 0.4, 0),
                            tcl = -0.25)
    on.exit(graphics::par(oldpar), add = TRUE)
    graphics::plot(c(),
                   xlim = c(0, max_limval),
                   ylim = c(0, max_limval),
                   axes = FALSE,
                   xlab = "Theoretical Quantiles",
                   ylab = "Observed P-values")
    graphics::box(lwd = 0.5)
    graphics::axis(side = 1, at = -log10(breakvec), labels = breakvec)
    graphics::axis(side = 2, at = -log10(breakvec), labels = breakvec)
    graphics::abline(h = -log10(breakvec), lwd = 1.5, col = "#B3B3B340")
    graphics::abline(v = -log10(breakvec), lwd = 1.5, col = "#B3B3B340")
    graphics::abline(a = 0, b = 1)
    graphics::points(x = theo_nl10, y = pvals_nl10, cex = 0.8)
    graphics::lines(x = theo_nl10, y = -log10(lower), lty = 2)
    graphics::lines(x = theo_nl10, y = -log10(upper), lty = 2)
    return(invisible(pvals))
  }
}
