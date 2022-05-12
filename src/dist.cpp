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

//' Dirichlet probability density function.
//'
//' @param x The observed vector of proportions. Should sum to 1.
//' @param alpha The concentation parameters, should all be greater than 0.
//' @param lg A logical. Should we return the log pdf or not?
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
double ddirichlet(NumericVector x, NumericVector alpha, bool lg = false) {
  if (Rcpp::is_false(Rcpp::all(alpha > 0.0))) {
    Rcpp::stop("alpha vector should be all positive");
  }
  if (std::abs(Rcpp::sum(x) - 1.0) > (100 * DBL_EPSILON)) {
    Rcpp::stop("x should sum to 1");
  }
  if (x.length() != alpha.length()) {
    Rcpp::stop("x and alpha should have same length");
  }

  double ll = std::lgamma(Rcpp::sum(alpha)) -
    Rcpp::sum(Rcpp::lgamma(alpha)) +
    Rcpp::sum((alpha - 1.0) * Rcpp::log(x));

  if (!lg) {
    ll = std::exp(ll);
  }

  return ll;
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
