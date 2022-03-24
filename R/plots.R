#' Critical values for the KS test
#'
#' Critical values taken from Birnbaum (1952). Spline interpolation chosen
#' for uncalculated sample sizes up to n = 100. Use 1.36/sqrt(n) for sample
#' sizes greater than 100. See, e.g. Aldor-Noiman et al (2013).
#'
#' @param n Sample size.
#'
#' @references
#' \itemize{
#'   \item{Aldor-Noiman, S., Brown, L. D., Buja, A., Rolke, W., & Stine, R. A. (2013). The power to see: A new graphical test of normality. The American Statistician, 67(4), 249-260.}
#'   \item{Birnbaum, Z. W. (1952). Numerical tabulation of the distribution of Kolmogorov's statistic for finite sample size. Journal of the American Statistical Association, 47(259), 425-441.}
#' }
#'
#' @noRd
#'
#' @examples
#' ks_crit(10)
#' ks_crit(1000)
#'
ks_crit <- function(n) {
  crit <- data.frame(val = c(0.8419,
                     0.7076,
                     0.6239,
                     0.5633,
                     0.4087,
                     0.3375,
                     0.2939,
                     0.2639,
                     0.2417,
                     0.2101,
                     0.1884,
                     0.1723,
                     0.1597,
                     0.1496,
                     0.1412,
                     0.1340),
             n = c(2,
                   3,
                   4,
                   5,
                   10,
                   15,
                   20,
                   25,
                   30,
                   40,
                   50,
                   60,
                   70,
                   80,
                   90,
                   100)
  )

  if (n <= 100) {
    ss <- stats::splinefun(x = crit$n, y = crit$val)
    cval <- ss(n)
  } else {
    cval <- 1.36 / sqrt(n)
  }
  return(cval)
}


#' Get bands for qqplot
#'
#' Procedure is described in Aldor-Noiman et al (2013). But note that they
#' have a mistake in their paper. Step (e) of their algorithm on page 254
#' should be the CDF of the Beta distribution, not the quantile function.
#'
#' @param n sample size
#' @param nsamp Number of simulation reps
#' @param a The significance level.
#'
#' @references
#' \itemize{
#'   \item{Aldor-Noiman, S., Brown, L. D., Buja, A., Rolke, W., & Stine, R. A. (2013). The power to see: A new graphical test of normality. The American Statistician, 67(4), 249-260.}
#' }
#'
#' @author David Gerard
#'
#' @noRd
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
  return(list(lower = lower, upper = upper))
}

#' QQ-plot for p-values
#'
#' This will create a QQ-plot for p-values, comparing them to a uniform
#' distribution. We make our plot on the -log10 scale. We calculate
#' the confidence bands by the Tail Sensative statistic of
#' Aldor-Noiman et al (2013).
#'
#' @param pvals A vector of p-values.
#' @param method Should we use base plotting or ggplot2 (if installed)?
#' @param band_type Should we use the method of Aldor-Noiman et al (2013) or
#'     pointwise based on beta? Pointwise is not recommended since there is
#'     strong dependence between order statistics and if one is beyond
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
#' qqpvalue(pvals, band_type = "pointwise", method = "base")
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
      (!requireNamespace(package = "ggplot2", quietly = TRUE) |
       !requireNamespace(package = "scales", quietly = TRUE))) {
    message("ggplot2 or scales not installed, using base")
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

  breakvec <- 10^-(0:ceiling(max(c(pvals_nl10, theo_nl10, -log10(bounds$lower)))))
  lim <- c(min(c(pvals, theo, bounds$lower)), max(c(pvals, theo, bounds$upper)))

  if(requireNamespace(package = "ggplot2", quietly = TRUE) &
     requireNamespace(package = "scales", quietly = TRUE) &
     method == "ggplot2") {
    nlog10 <- scales::trans_new(name = "nlog10",
                      transform = function(x) -log10(x),
                      inverse = function(x) 10 ^ -x,
                      domain = c(0, 1)
                      )

    pl <- ggplot2::ggplot(data.frame(theo = theo,
                                     pvals = pvals,
                                     lower = bounds$lower,
                                     upper = bounds$upper)) +
      ggplot2::coord_trans(x = nlog10, y = nlog10) +
      ggplot2::geom_point(mapping = ggplot2::aes(x = theo, y = pvals)) +
      ggplot2::geom_abline() +
      ggplot2::geom_line(mapping = ggplot2::aes(x = theo, y = lower)) +
      ggplot2::geom_line(mapping = ggplot2::aes(x = theo, y = upper)) +
      ggplot2::scale_x_continuous(name = "Theoretical Quantiles",
                                  breaks = breakvec,
                                  limits = lim) +
      ggplot2::scale_y_continuous(name = "Observed P-values",
                                  breaks = breakvec,
                                  limits = lim) +
      ggplot2::theme_bw()
    print(pl)
    return(invisible(pvals))
  } else {
    oldpar <- graphics::par(pch = 16,
                            mar = c(3, 3, 2, 1),
                            mgp = c(1.8, 0.4, 0),
                            tcl = -0.25)
    on.exit(graphics::par(oldpar), add = TRUE)
    graphics::plot(c(),
                   xlim = rev(-log10(lim)),
                   ylim = rev(-log10(lim)),
                   axes = FALSE,
                   xlab = "Theoretical Quantiles",
                   ylab = "Observed P-values")
    graphics::axis(side = 1, at = -log10(breakvec), labels = breakvec)
    graphics::axis(side = 2, at = -log10(breakvec), labels = breakvec)
    graphics::abline(a = 0, b = 1)
    graphics::points(x = theo_nl10, y = pvals_nl10)
    graphics::lines(x = theo_nl10, y = -log10(lower))
    graphics::lines(x = theo_nl10, y = -log10(upper))
    return(invisible(pvals))
  }
}
