#include <Rcpp.h>
using namespace Rcpp;

//' Sample gamete counts from full conditional.
//'
//' @param x A vector of genotype counts. x[i] is the number of individuals
//'     with dosage i.
//' @param p A vector of gamete probabilities. p[i] is the proportion of
//'     gametes with dosage i.
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
IntegerVector samp_gametes(const NumericVector& x,
                           const NumericVector& p) {
  int ploidy = x.length() - 1;
  if (p.length() != (ploidy / 2 + 1)) {
    Rcpp::stop("p should be the same length as ploidy / 2 + 1");
  }
  if (std::abs(Rcpp::sum(p) - 1.0) > (DBL_EPSILON * 100)) {
    Rcpp::stop("p should sum to 1");
  }

  IntegerMatrix A(ploidy / 2 + 1, ploidy / 2 + 1);
  IntegerVector y(ploidy / 2 + 1);

  for (int k = 0; k <= ploidy; k++) {
    int lind = std::max(0, k - ploidy / 2);
    int uind = std::floor((double)k / 2.0);
    int numc = uind - lind + 1;
    IntegerVector ak(numc);
    NumericVector rk(numc);

    // populate rk
    int i = 0;
    for (int j = lind; j <= uind; j++) {
      if (j == k - j) {
        rk(i) = p(j) * p(j);
      } else {
        rk(i) = 2.0 * p(j) * p(k - j);
      }
      i++;
    }
    rk = rk / sum(rk);

    // simulate ak
    R::rmultinom(x(k), rk.begin(), numc, ak.begin());

    // populate A
    i = 0;
    for (int j = lind; j <= uind; j++) {
      A(j, k-j) = ak(i);
      i++;
    }
  }

  // populate y
  for (int k = 0; k <= ploidy / 2; k++) {
    y(k) = 0;
    for (int j = 0; j <= ploidy / 2; j++) {
      if (j < k) {
        y(k) += A(j, k);
      } else if (j == k) {
        y(k) += 2 * A(j, k);
      } else {
        y(k) += A(k, j);
      }
    }
  }

  return y;
}


//' Gibbs sampler under random mating with known genotypes.
//'
//' @param x The vector of genotype counts. x(i) is the number of
//'     individuals that have genotype i.
//' @param alpha Vector of hyperparameters for the gamete frequencies.
//'     Should be length (x.length() - 1) / 2 + 1.
//' @param B The number of sampling iterations.
//' @param T The number of burn-in iterations.
//' @param more A logical. Should we also return posterior draws (\code{TRUE})
//'     or not (\code{FALSE}).
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
Rcpp::List gibbs_known(Rcpp::NumericVector x,
                       Rcpp::NumericVector alpha,
                       int B,
                       int T,
                       bool more) {
  int ploidy = x.length() - 1;
  if (alpha.length() != (ploidy / 2 + 1)) {
    Rcpp::stop("alpha should be the same length as ploidy / 2 + 1");
  }

  Rcpp::List retlist;

  return retlist;
}






