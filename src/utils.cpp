#include <Rcpp.h>
using namespace Rcpp;

//' Log-sum-exponential trick using just two doubles.
//'
//' @param x A double.
//' @param y Another double.
//'
//' @return The log of the sum of the exponential of x and y.
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
double log_sum_exp_2_cpp(double x, double y) {
  double z = std::max(x, y);
  double finalval;
  if (z == R_NegInf) {
    finalval = R_NegInf;
  } else {
    finalval = std::log(std::exp(x - z) + std::exp(y - z)) + z;
  }
  return finalval;
}

//' Log-sum-exponential trick.
//'
//' @param x A vector to log-sum-exp.
//'
//' @return The log of the sum of the exponential
//'     of the elements in \code{x}.
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
double log_sum_exp_cpp(const NumericVector& x) {
  double max_x = Rcpp::max(x);
  double lse; // the log-sum-exp
  // if all -Inf, need to treat this special to avoid -Inf + Inf.
  if (max_x == R_NegInf) {
    lse = R_NegInf;
  } else {
    lse = max_x + std::log(Rcpp::sum(Rcpp::exp(x - max_x)));
  }
  return lse;
}

//' Basic convolution using naive algorithm
//'
//' Same as stats::conv(x, rev(y), type = "open"), but using the naive
//' algorithm. This should be faster for very small vectors (but will
//' be horrible for larger vectors).
//'
//' @param x The first vector
//' @param y The second vector
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
NumericVector conv_cpp(NumericVector x, NumericVector y) {
  NumericVector z(x.length() + y.length() - 1);
  z.fill(0.0);

  for (int k = 0; k < z.length(); k++) {
    // maximum of k and (y.length() - 1)
    int lind;
    if (k > (y.length() - 1)) {
      lind = k - y.length() + 1;
    } else {
      lind = 0;
    }

    // minimum of k and (x.length() - 1)
    int uind;
    if (k < (x.length() - 1)) {
      uind = k;
    } else {
      uind = x.length() - 1;
    }

    for (int i = lind; i <= uind; i++) {
      z(k) += x(i) * y(k - i);
    }
  }

  return z;
}
