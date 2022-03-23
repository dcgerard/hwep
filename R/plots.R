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
    ss <- spline(x = crit$n, y = crit$val)
    ss <- stats::splinefun(x = crit$n, y = crit$val)
    cval <- ss(n)
  } else {
    cval <- 1.36 / sqrt(n)
  }
  return(cval)
}

#' QQ-plot for p-values
#'
#' This will create a QQ-plot for p-values, comparing them to a uniform
#' distribution. We make our plot on the -log10 scale. We calculate
#' the confidence bands by inverting the KS test. See, e.g.
#' Aldor-Noiman et al (2013).
#'
#' @param pvals A vector of p-values.
#' @param method Should we use base plotting or ggplot2 (if installed)?
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
#' qqpvalue(pvals, method = "base")
#'
#' @export
qqpvalue <- function(pvals,
                     method = c("ggplot2", "base")) {
  method <- match.arg(method)
  if (method == "ggplot2" &
      (!requireNamespace(package = "ggplot2", quietly = TRUE) |
       !requireNamespace(package = "scales", quietly = TRUE))) {
    warning("ggplot2 or scales not installed, using base")
    method <- "base"
  }

  pvals <- sort(pvals)
  theo <- stats::ppoints(n = length(pvals))

  pvals_nl10 <- -log10(pvals)
  theo_nl10 <- -log10(theo)

  breakvec <- 10^-(0:ceiling(max(c(pvals_nl10, theo_nl10))))
  lim <- c(min(c(pvals, theo)), max(c(pvals, theo)))

  if(requireNamespace(package = "ggplot2", quietly = TRUE) &
     requireNamespace(package = "scales", quietly = TRUE) &
     method == "ggplot2") {
    nlog10 <- scales::trans_new(name = "nlog10",
                      transform = function(x) -log10(x),
                      inverse = function(x) 10 ^ -x,
                      domain = c(0, 1)
                      )

    pl <- ggplot2::ggplot(data.frame(theo = theo,
                                     pvals = pvals)) +
      ggplot2::coord_trans(x = nlog10, y = nlog10) +
      ggplot2::geom_point(mapping = ggplot2::aes(x = theo, y = pvals)) +
      ggplot2::geom_abline() +
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
    return(invisible(pvals))
  }
}
