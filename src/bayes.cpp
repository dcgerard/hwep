#include <Rcpp.h>
using namespace Rcpp;


//' Random sample from Dirichlet distribution with n = 1.
//'
//' @param alpha The
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
NumericVector rdirichlet1(NumericVector alpha) {
  int n = alpha.length();
  NumericVector x(n);
  for (int i = 0; i < n; i++) {
    x(i) = R::rgamma(alpha(i), 1.0);
  }
  x = x / Rcpp::sum(x);
  return x;
}


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
IntegerVector samp_gametes(NumericVector x, NumericVector p) {
  int ploidy = x.length() - 1;
  if (p.length() != ploidy / 2 + 1) {
    Rcpp::stop("p should be the same length as ploidy / 2");
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
