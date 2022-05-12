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

//' Multinomial probability mass function
//'
//' @param x The counts.
//' @param p The probabilities.
//' @param lg A logical. Should we return log (true) or not (false).
//'
//' @noRd
//' @author David Gerard
// [[Rcpp::export]]
double dmultinom_cpp(NumericVector x, NumericVector p, bool lg = false) {
  double n = sum(x);
  double ret;
  double TOL = DBL_EPSILON * 100;
  double sval = 0.0;
  for (int i = 0; i < x.length(); i++) {
    if ((p(i) > TOL) || (x(i) > TOL)) {
      sval += x(i) * log(p(i));
    }
  }
  ret = std::lgamma(n + 1.0) - sum(lfactorial(x)) + sval;

  if (!lg) {
    ret = std::exp(ret);
  }
  return ret;
}
