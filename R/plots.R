#' QQ-plot for p-values
#'
#' This will create a QQ-plot for p-values, comparing them to a unifrom
#' distribution. We make our plot on the -log10 scale. We calculate
#' the confidence envelope using the method of Fox (2016), pp 37--41.
#'
#' @param pvals A vector of p-values.
#' @param method Should we use base plotting or ggplot2 (if installed)?
#'
#' @author David Gerard
#'
#' @references
#' \itemize{
#'   \item{Fox, J. (2016) \emph{Applied Regression Analysis and Generalized Linear Models}, Third Edition. Sage.}
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
#' qqpvalue(pvals, method = "ggplot2")
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
                                     pvals = pvals),
                          mapping = ggplot2::aes(x = theo, y = pvals)) +
      ggplot2::coord_trans(x = nlog10, y = nlog10) +
      ggplot2::geom_point() +
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
