########################
## Theoretical genotype frequencies from Huang et al (2019)
########################

################
## The following are functions from Supplementary Material S3 of Huang et al (2019)
## https://doi.org/10.1534/g3.119.400132
## That I have converted from C code to R code
################

GFG4_ii <- function(a1, pi) {
	(pi*(a1*(3.0-2.0*pi)+2.0*pi))/(2.0+a1)
}

GFG4_ij <- function(a1, pi, pj) {
	-((4.0*(-1.0+a1)*pi*pj)/(2.0+a1))
}

GFG6_iii <- function(a1, pi) {
	pi2 <- pi*pi
	a12 <- a1*a1
	(pi*(81.0*pi2-27.0*a1*pi*(-5.0+4.0*pi)+a12*(20.0-45.0*pi+27.0*pi2)))/(81.0+27.0*a1+2.0*a12)
}

GFG6_iij <- function(a1, pi, pj) {
	a12 <- a1*a1
	(9.0*(-3.0+a1)*pi*(-9.0*pi+a1*(-5.0+9.0*pi))*pj)/(81.0+27.0*a1+2.0*a12)
}

GFG8_iiii <- function(a1, a2, pi) {
	pi2 <- pi*pi
	pi3 <- pi2*pi
	a12 <- a1*a1
	a13 <- a12*a1
	a14 <- a13*a1
	a22 <- a2*a2
	a23 <- a22*a2
	a24 <- a23*a2
	C1 <- (-35.0+210.0*pi-336.0*pi2+160.0*pi3)
	-(1.0/((34.0+a1+a2)*(8.0+a1+2.0*a2)^2*(a1+2.0*(6.0+a2))))*pi*(a14*C1+a13*(-2.0*(140.0-301.0*pi+112.0*pi2+80.0*pi3)+7.0*a2*C1)+2.0*a12*(-4.0*pi*(1519.0-2632.0*pi+1264.0*pi2)+9.0*a22*C1-2.0*a2*(455.0-945.0*pi+336.0*pi2+240.0*pi3))+8.0*(-3264.0*pi3+32.0*a2*pi*(-105.0-147.0*pi+181.0*pi2)-2.0*a23*(175.0-343.0*pi+112.0*pi2+80.0*pi3)+a24*C1-4.0*a22*(210.0+609.0*pi-1316.0*pi2+632.0*pi3))+4.0*a1*(32.0*pi2*(-357.0+283.0*pi)+5.0*a23*C1-2.0*a22*(490.0-987.0*pi+336.0*pi2+240.0*pi3)-8.0*a2*(105.0+1064.0*pi-1974.0*pi2+948.0*pi3)))
}

GFG8_iiij <- function(a1, a2, pi, pj)
{
	pi2 <- pi*pi
	a12 <- a1*a1
	a13 <- a12*a1
	a14 <- a13*a1
	a22 <- a2*a2
	a23 <- a22*a2
	C1 <- (7.0-42.0*pi+40.0*pi2)
	-(1.0/((34.0+a1+a2)*(8.0+a1+2.0*a2)^2*(a1+2.0*(6.0+a2))))*8.0*pi*(2.0*a14*C1+a13*(-7.0-56.0*pi-80.0*pi2+14.0*a2*C1)+16.0*(-1.0+a2)*(a22*(28.0-70.0*pi)+4.0*a2*(147.0-158.0*pi)*pi+816.0*pi2+a23*C1)+4.0*a12*(9.0*a22*C1-2.0*a2*(-7.0+42.0*pi+60.0*pi2)-2.0*(119.0-658.0*pi+632.0*pi2))+4.0*a1*(8.0*pi*(-357.0+566.0*pi)+a22*(77.0-168.0*pi-240.0*pi2)+10.0*a23*C1-4.0*a2*(133.0-987.0*pi+948.0*pi2)))*pj
}

GFG8_iijj <- function(a1, a2, pi, pj) {
	a12 <- a1*a1
	a13 <- a12*a1
	a14 <- a13*a1
	a22 <- a2*a2
	a23 <- a22*a2
	a24 <- a23*a2
	C1 <- (49.0-84.0*pj+12.0*pi*(-7.0+20.0*pj))
	C2 <- (7.0+30.0*pj)
	C3 <- (-329.0+948.0*pj)
	-(1.0/((34.0+a1+a2)*(8.0+a1+2.0*a2)^2*(a1+2.0*(6.0+a2))))*4.0*pi*pj*(a14*C1+a13*(329.0-56.0*pj-8.0*pi*C2+7.0*a2*C1)+2.0*a12*(-1134.0+2632.0*pi+2632.0*pj-7584.0*pi*pj+9.0*a22*C1-a2*(-833.0+168.0*pj+24.0*pi*C2))+8.0*(-4896.0*pi*pj+24.0*a2*(-70.0-49.0*pi-49.0*pj+362.0*pi*pj)+a24*C1-a23*(7.0*(-25.0+8.0*pj)+8.0*pi*C2)-2.0*a22*(497.0-658.0*pj+2.0*pi*C3))+4.0*a1*(5.0*a23*C1-a22*(7.0*(-97.0+24.0*pj)+24.0*pi*C2)+24.0*(-119.0*pj+pi*(-119.0+566.0*pj))-4.0*a2*(532.0-987.0*pj+3.0*pi*C3)))
}

GFG10_iiiii <- function(a1, a2, pi) {
	pi2 <- pi*pi
	pi3 <- pi2*pi
	pi4 <- pi3*pi
	a12 <- a1*a1
	a13 <- a12*a1
	a14 <- a13*a1
	a15 <- a14*a1
	a22 <- a2*a2
	a23 <- a22*a2
	a24 <- a23*a2
	a25 <- a24*a2
	C1 <- (a1+2.0*(-5.0+a2))
	C2 <- (20.0+a2)
	(pi*(3125.0*pi4*C1*(-1.0+a1+a2)*(25000.0+a1*(-11050.0+9.0*a1*(75.0+7.0*a1))-21300.0*a2+a1*(2650.0+357.0*a1)*a2+8.0*(325.0+84.0*a1)*a22+420.0*a23)-750.0*pi*(45.0*a13*(-992.0+a1*(43.0+4.0*a1))+2.0*a1*(-133000.0+a1*(-104590.0+a1*(8039.0+780.0*a1)))*a2+2.0*(-255500.0+a1*(-154580.0+a1*(24727.0+2670.0*a1)))*a22+40.0*(-3468.0+5.0*a1*(334.0+45.0*a1))*a23+8.0*(4181.0+930.0*a1)*a24+2400.0*a25)+504.0*(a1+2.0*a2)*(25.0+2.0*a1+4.0*a2)*(9.0*a13+42.0*a12*a2+30.0*a22*C2+7.0*a1*a2*(50.0+9.0*a2))-18750.0*pi3*C1*(27.0*a14+6.0*a13*(59.0+30.0*a2)+10.0*(-1.0+a2)*a2*(-625.0+a2*(173.0+18.0*a2))+a12*(-4290.0+7.0*a2*(256.0+63.0*a2))+a1*(7500.0+a2*(-12395.0+9.0*a2*(327.0+52.0*a2))))+375.0*pi2*(1161.0*a15+234.0*a14*(34.0+43.0*a2)+9.0*a13*(-36655.0+a2*(7172.0+3827.0*a2))+2.0*a1*a2*(1510375.0+3.0*a2*(-430805.0+42236.0*a2+7998.0*a22))+10.0*a12*(129375.0+a2*(-162814.0+15.0*a2*(1289.0+387.0*a2)))+20.0*a2*(109375.0+a2*(50850.0+a2*(-64813.0+6126.0*a2+774.0*a22))))))/(2.0*(125.0+a1+a2)*(25.0+2.0*a1+4.0*a2)^2*(50.0+3.0*a1+6.0*a2)*(3.0*a1+5.0*C2))
}

GFG10_iiiij <- function(a1, a2, pi, pj) {
	pi2 <- pi*pi
	pi3 <- pi2*pi
	a12 <- a1*a1
	a13 <- a12*a1
	a14 <- a13*a1
	a15 <- a14*a1
	a22 <- a2*a2
	a23 <- a22*a2
	C1 <- (a1+2.0*(-5.0+a2))
	-(25.0*pi*(-625.0*pi3*C1*(-1.0+a1+a2)*(25000.0+a1*(-11050.0+9.0*a1*(75.0+7.0*a1))-21300.0*a2+a1*(2650.0+357.0*a1)*a2+8.0*(325.0+84.0*a1)*a22+420.0*a23)+84.0*(a1+2.0*a2)*(25.0+2.0*a1+4.0*a2)*(9.0*a13+3.0*a12*(-75.0+14.0*a2)+a1*a2*(-415.0+63.0*a2)+10.0*a2*(-125.0+a2*(-4.0+3.0*a2)))+2250.0*pi2*C1*(27.0*a14+6.0*a13*(59.0+30.0*a2)+10.0*(-1.0+a2)*a2*(-625.0+a2*(173.0+18.0*a2))+a12*(-4290.0+7.0*a2*(256.0+63.0*a2))+a1*(7500.0+a2*(-12395.0+9.0*a2*(327.0+52.0*a2))))-15.0*pi*(1593.0*a15+6.0*a14*(1165.0+2301.0*a2)+a13*(-442875.0+a2*(62840.0+47259.0*a2))+10.0*a12*(189375.0+a2*(-215540.0+a2*(20443.0+7965.0*a2)))+2.0*a1*a2*(2181875.0+a2*(-1683475.0+2.0*a2*(71830.0+16461.0*a2)))+20.0*a2*(109375.0+a2*(67250.0+a2*(-82765.0+2.0*a2*(3695.0+531.0*a2))))))*pj)/(2.0*(125.0+a1+a2)*(25.0+2.0*a1+4.0*a2)^2*(50.0+3.0*a1+6.0*a2)*(3.0*a1+5.0*(20.0+a2)))
}

GFG10_iiijj <- function(a1, a2, pi, pj) {
	pi2 <- pi*pi
	a12 <- a1*a1
	a13 <- a12*a1
	a14 <- a13*a1
	a15 <- a14*a1
	a22 <- a2*a2
	a23 <- a22*a2
	C1 <- (a1+2.0*(-5.0+a2))
	C2 <- (27.0*a14+6.0*a13*(59.0+30.0*a2)+10.0*(-1.0+a2)*a2*(-625.0+a2*(173.0+18.0*a2))+a12*(-4290.0+7.0*a2*(256.0+63.0*a2))+a1*(7500.0+a2*(-12395.0+9.0*a2*(327.0+52.0*a2))))
	C3 <- (20.0+a2)
	(25.0*pj*(3.0*pi*(-125.0*pi2*C1*C2+3.0*(-216.0*a15-9.0*a14*(475.0+208.0*a2)-160.0*a22*C3*(-175.0+2.0*a2*(-23.0+9.0*a2))-10.0*a12*a2*(-23605.0+2.0*a2*(4477.0+540.0*a2))-6.0*a13*(-8025.0+a2*(5345.0+1068.0*a2))+4.0*a1*a2*(74375.0+a2*(91775.0-2.0*a2*(13765.0+1116.0*a2))))+5.0*pi*(729.0*a15+a14*(8922.0+6318.0*a2)+a13*(-216915.0+a2*(66256.0+21627.0*a2))+10.0*a12*(69375.0+a2*(-110088.0+a2*(18227.0+3645.0*a2)))+2.0*a1*a2*(838875.0+a2*(-901355.0+2.0*a2*(54878.0+7533.0*a2)))+20.0*a2*(109375.0+a2*(34450.0+a2*(-46861.0+4862.0*a2+486.0*a22)))))+5.0*pi*(648.0*a15+9.0*a14*(-161.0+624.0*a2)+125.0*pi2*C1*(-1.0+a1+a2)*(25000.0+a1*(-11050.0+9.0*a1*(75.0+7.0*a1))-21300.0*a2+a1*(2650.0+357.0*a1)*a2+8.0*(325.0+84.0*a1)*a22+420.0*a23)+480.0*(-1.0+a2)*a22*(-1025.0+a2*(97.0+18.0*a2))+6.0*a13*(-28245.0+a2*(-427.0+3204.0*a2))+12.0*a1*a2*(167875.0+a2*(-97765.0+4238.0*a2+2232.0*a22))+30.0*a12*(30000.0+a2*(-26363.0+2.0*a2*(277.0+540.0*a2)))-225.0*pi*C1*C2)*pj))/((125.0+a1+a2)*(25.0+2.0*a1+4.0*a2)^2*(50.0+3.0*a1+6.0*a2)*(3.0*a1+5.0*C3))
}


#' Theoretical frequencies at equilibrium.
#'
#' These return gamete and genotype frequencies as calculated in
#' Huang et al (2019). Only supported for ploidies less than or equal to 10.
#'
#' @param alpha A numeric vector containing the double reduction parameter(s).
#'     This should be a
#'     vector of length \code{floor(ploidy/4)} where \code{alpha[i]}
#'     is the probability of exactly \code{i} pairs of IBDR alleles
#'     being in the gamete. Note that \code{sum(alpha)} should be less than
#'     1, as \code{1 - sum(alpha)} is the probability of no double reduction.
#' @param ploidy The ploidy of the species. This should be an even positive
#'     integer.
#' @param r The allele frequency of the reference allele.
#'
#' @references
#' \itemize{
#'   \item{Huang, K., Wang, T., Dunn, D. W., Zhang, P., Cao, X., Liu, R., & Li, B. (2019). Genotypic frequencies at equilibrium for polysomic inheritance under double-reduction. G3: Genes, Genomes, Genetics, 9(5), 1693-1706. \doi{10.1534/g3.119.400132}}
#' }
#'
#' @return A list with the following elements
#' \describe{
#'   \item{\code{p}}{The gamete frequencies at equilibrium.}
#'   \item{\code{q}}{The genotype frequencies at equilibrium.}
#' }
#'
#' @author David Gerard
#'
#' @noRd
theofreq <- function(alpha, r, ploidy) {
  stopifnot(length(r) == 1L, length(ploidy) == 1L)
  stopifnot(ploidy %% 2 == 0)
  stopifnot(ploidy > 1, ploidy < 11)
  stopifnot(length(alpha) == floor(ploidy / 4))
  stopifnot(sum(alpha) <= 1)

  stopifnot(alpha > -sqrt(.Machine$double.eps))

  alpha[alpha < 0] <- 0

  if (ploidy == 2) {
    hout <- c(r, 1-r)
  } else if (ploidy == 4) {
    hout <- c(
      GFG4_ii(a1 = alpha, pi = 1 - r),
      GFG4_ij(a1 = alpha, pi = 1 - r, pj = r),
      GFG4_ii(a1 = alpha, pi = r)
    )
  } else if (ploidy == 6) {
    hout <- c(
      GFG6_iii(a1 = alpha, pi = 1 - r),
      GFG6_iij(a1 = alpha, pi = 1 - r, pj = r),
      GFG6_iij(a1 = alpha, pi = r, pj = 1 - r),
      GFG6_iii(a1 = alpha, pi = r)
    )
  } else if (ploidy == 8) {
    hout <- c(
      GFG8_iiii(a1 = alpha[[1]], a2 = alpha[[2]], pi = 1 - r),
      GFG8_iiij(a1 = alpha[[1]], a2 = alpha[[2]], pi = 1 - r, pj = r),
      GFG8_iijj(a1 = alpha[[1]], a2 = alpha[[2]], pi = 1 - r, pj = r),
      GFG8_iiij(a1 = alpha[[1]], a2 = alpha[[2]], pi = r, pj = 1 - r),
      GFG8_iiii(a1 = alpha[[1]], a2 = alpha[[2]], pi = r)
    )
  } else if (ploidy == 10) {
    hout <- c(
      GFG10_iiiii(a1 = alpha[[1]], a2 = alpha[[2]], pi = 1 - r),
      GFG10_iiiij(a1 = alpha[[1]], a2 = alpha[[2]], pi = 1 - r, pj = r),
      GFG10_iiijj(a1 = alpha[[1]], a2 = alpha[[2]], pi = 1 - r, pj = r),
      GFG10_iiijj(a1 = alpha[[1]], a2 = alpha[[2]], pi = r, pj = 1 - r),
      GFG10_iiiij(a1 = alpha[[1]], a2 = alpha[[2]], pi = r, pj = 1 - r),
      GFG10_iiiii(a1 = alpha[[1]], a2 = alpha[[2]], pi = r)
    )
  }

  ## Fix numerical issues
  stopifnot(hout > -sqrt(.Machine$double.eps))
  hout[hout < 0] <- 0
  hout <- hout / sum(hout)
  q <- stats::convolve(x = hout, y = rev(hout), type = "open")
  q[q < 0] <- 0
  q <- q / sum(q)

  retlist <- list(q = q,
                  p = hout)

  return(retlist)
}

#' log-likelihood used in hwelike
#'
#' @param par first element is r, rest are alpha
#' @param nvec the counts
#' @param which_keep A logical vector the same length as nvec, indicating
#'     which genotypes to aggregate.
#'
#' @author David Gerard
#'
#' @noRd
like_obj <- function(par, nvec, which_keep = NULL) {
  if (is.null(which_keep)) {
    which_keep <- rep(TRUE, length(nvec))
  }

  ploidy <- length(nvec) - 1
  r <- par[[1]]
  alpha <- par[-1]
  fq <- theofreq(alpha = alpha, r = r, ploidy = ploidy)
  q <- fq$q

  if (!all(which_keep)) {
    nvec <- c(nvec[which_keep], sum(nvec[!which_keep]))
    q <- c(q[which_keep], sum(q[!which_keep]))
  }

  return(stats::dmultinom(x = nvec, prob = q, log = TRUE))
}


#' log-likelihood used in hwelike
#'
#' @param alpha double reduction parameter
#' @param r allele frequency
#' @param nvec the counts
#' @param which_keep A logical vector the same length as nvec, indicating
#'     which genotypes to aggregate.
#'
#' @author David Gerard
#'
#' @noRd
like_obj2 <- function(alpha, r, nvec, which_keep = NULL) {
  like_obj(par = c(r, alpha), nvec = nvec, which_keep = which_keep)
}

#' Maximum likelihood approach for equilibrium testing and double reduction
#' estimation.
#'
#' Genotype frequencies from Huang et al (2019) are used to implement a
#' likelihood procedure to estimate double reduction rates and to test
#' for equilibrium while accounting for double reduction. This approach
#' is only implemented for ploidies 4, 6, 8, and 10.
#'
#' @inheritParams hweustat
#'
#' @author David Gerard
#'
#' @references
#' \itemize{
#'   \item{Huang, K., Wang, T., Dunn, D. W., Zhang, P., Cao, X., Liu, R., & Li, B. (2019). Genotypic frequencies at equilibrium for polysomic inheritance under double-reduction. G3: Genes, Genomes, Genetics, 9(5), 1693-1706. \doi{10.1534/g3.119.400132}}
#' }
#'
#' @return A list with some or all of the following elements:
#' \describe{
#'   \item{\code{alpha}}{The estimated double reduction parameter(s).
#'       In diploids, this value is \code{NULL}.}
#'   \item{\code{r}}{The estimated allele frequency.}
#'   \item{\code{chisq_hwe}}{The chi-square test statistic for testing
#'       against the null of equilibrium.}
#'   \item{\code{df_hwe}}{The degrees of freedom associated with
#'       \code{chisq_hwe}.}
#'   \item{\code{p_hwe}}{The p-value against the null of equilibrium.}
#' }
#'
#' @export
#'
#' @examples
#' thout <- hwefreq(alpha = 0.1, r = 0.3, ploidy = 6)
#' nvec <- c(stats::rmultinom(n = 1, size = 100, prob = thout))
#' hwelike(nvec = nvec)
#'
hwelike <- function(nvec,
                    thresh = 5,
                    effdf = FALSE) {
  ploidy <- length(nvec) - 1
  stopifnot(ploidy %% 2 == 0, ploidy >= 4, ploidy <= 10)
  stopifnot(nvec >= 0)
  stopifnot(length(thresh) == 1, thresh >= 0)
  stopifnot(is.logical(effdf), length(effdf) == 1)
  ibdr <- floor(ploidy / 4)
  omethod <- ifelse(ibdr == 1, "Brent", "L-BFGS-B")

  ## Choose which groups to aggregate ----
  which_keep <- choose_agg(x = nvec, thresh = thresh, like = TRUE)

  numkeep <- sum(which_keep)
  if (numkeep >= ploidy) {
    ## aggregating one group = no aggregation at all.
    which_keep <- rep(TRUE, ploidy + 1)
    numkeep <- ploidy + 1
  }

  ## Return early if too few groups ----
  if (sum(which_keep) - ibdr <= 0) {
    return(
      list(
        alpha = rep(NA_real_, ibdr),
        r = NA_real_,
        chisq_hwe = NA_real_,
        df_hwe = sum(which_keep) - ibdr - 1,
        p_hwe = NA_real_
      )
    )
  }

  ## Find MLE under null ----
  rhat <- sum(0:ploidy * nvec) / (ploidy * sum(nvec))
  minval <- 0.0001
  upper_alpha <- drbounds(ploidy = ploidy)
  oout <- stats::optim(par = rep(minval, ibdr),
                       fn = like_obj2,
                       method = omethod,
                       lower = rep(minval, ibdr),
                       upper = upper_alpha,
                       control = list(fnscale = -1),
                       nvec = nvec,
                       r = rhat,
                       which_keep = which_keep)

  alphahat <- oout$par
  ll_e <- oout$value

  ## MLE under alternative ----
  if (all(which_keep)) {
    q_u <- nvec / sum(nvec)
    ll_u <- stats::dmultinom(x = nvec, prob = q_u, log = TRUE)
  } else {
    ntemp <- c(nvec[which_keep], sum(nvec[!which_keep]))
    q_temp <- ntemp / sum(ntemp)
    ll_u <- stats::dmultinom(x = ntemp, prob = q_temp, log = TRUE)
  }

  ## Test statistic ----
  chisq_hwe <- -2 * (ll_e - ll_u)

  ## Find degrees of freedom ----
  if (all(which_keep)) {
    df_hwe <- ploidy - ibdr - 1
  } else {
    df_hwe <- sum(which_keep) - ibdr - 1 ## unconstrained as sum(which_keep), alpha is ibdr, r is 1.
  }

  TOL <- sqrt(.Machine$double.eps)
  if (effdf) {
    dfadd <- sum((abs(alphahat - minval) < TOL) | (abs(alphahat - upper_alpha) < TOL))
  } else {
    dfadd <- 0
  }
  df_hwe <- df_hwe + dfadd

  ## Run test ----
  if (df_hwe > 0) {
    p_hwe <- stats::pchisq(q = chisq_hwe, df = df_hwe, lower.tail = FALSE)
  } else {
    p_hwe <- NA_real_
  }

  retlist <- list(alpha = alphahat,
                  r = rhat,
                  chisq_hwe = chisq_hwe,
                  df_hwe = df_hwe,
                  p_hwe = p_hwe)

  return(retlist)
}
