#include <Rcpp.h>
using namespace Rcpp;

// Functions from utils.cpp -------
NumericVector conv_cpp(NumericVector x, NumericVector y);
double log_sum_exp_2_cpp(double x, double y);
double log_sum_exp_cpp(const NumericVector& x);

// Functions from dist.cpp --------
NumericVector rdirichlet1(NumericVector alpha);
double ddirichlet(NumericVector x, NumericVector alpha, bool lg);
double dmultinom_cpp(NumericVector x, NumericVector p, bool lg);

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
//' @param lg Should we return the log marginal likelihood (true) or not
//'     (false).
//'
//' @return A list with some or all of the following elements
//' \itemize{
//'   \item{\code{mx}: The estimate of the marginal likelihood}
//'   \item{\code{p_tilde}: The value of p used to evaluate the posterior density}.
//'   \item{\code{p}: The samples of the gamete frequencies}
//'   \item{\code{post}: The samples of the full conditionals of p_tilde.}
//' }
//'
//' @author David Gerard
//'
//' @export
// [[Rcpp::export]]
Rcpp::List gibbs_known(Rcpp::NumericVector x,
                       Rcpp::NumericVector alpha,
                       int B = 10000,
                       int T = 100,
                       bool more = false,
                       bool lg = false) {
  int ploidy = x.length() - 1;
  if (alpha.length() != (ploidy / 2 + 1)) {
    Rcpp::stop("alpha should be the same length as ploidy / 2 + 1");
  }

  NumericVector y(alpha.length());
  NumericVector p = rdirichlet1(alpha);
  NumericVector p_tilde = p;
  NumericVector q = conv_cpp(p_tilde, p_tilde);
  double logpost = dmultinom_cpp(x, q, true) +  ddirichlet(p_tilde, alpha, true);
  double logpihat = R_NegInf;

  // build more output ----
  int nsamp;
  if (more) {
    nsamp = B;
  } else {
    nsamp = 0;
  }
  NumericMatrix pmat(nsamp, ploidy / 2 + 1);
  NumericVector postvec(nsamp);

  for (int i = 0; i < T + B; i++) {
    // sample y
    y = as<NumericVector>(samp_gametes(x, p));

    // sample p
    p = rdirichlet1(y + alpha);

    // Check if we are in burnin
    if (i < T) {
      q = conv_cpp(p, p);
      double lp2 = dmultinom_cpp(x, q, true) +  ddirichlet(p, alpha, true);
      if (lp2 > logpost) {
        p_tilde = p;
        logpost = lp2;
      }
    } else {
      double ptilde_post = ddirichlet(p_tilde, y + alpha, true);
      logpihat = log_sum_exp_2_cpp(logpihat, ptilde_post);

      // include if more
      if (more) {
        NumericMatrix::Row crow = pmat(i - T, _);
        crow = p;
        postvec(i - T) = ptilde_post;
      }

    }
  }
  logpihat -= log((double)B);

  Rcpp::List retlist;
  if (lg) {
    retlist["mx"] = logpost - logpihat;
  } else {
    retlist["mx"] = exp(logpost - logpihat);
  }

  // include if more
  if (more) {
    retlist["p_tilde"] = p_tilde;
    retlist["p"] = pmat;
    retlist["post"] = postvec;
  }

  return retlist;
}

//' Calculate marginal likelihood under alternative using genotype likelihoods.
//'
//' Calculates
//' \deqn{
//' \log \prod_i \sum_k l_{ik}q_{ik}
//' }
//' where q = beta / sum(beta). It can return the exponentiated version of this.
//'
//' @param gl genotype log-liklihoods. Rows index individuals, columns
//'     index genotypes.
//' @param beta The concentration parameters.
//' @param lg Should we return the log (true) or the not (false)?
//'
//' @author David Gerard
//'
//' @noRd
//'
// [[Rcpp::export]]
double plq(NumericMatrix& gl, NumericVector beta, bool lg = false) {
  NumericVector lq = Rcpp::log(beta / Rcpp::sum(beta));
  if (gl.ncol() != beta.length()) {
    Rcpp::stop("Number of columns of gl should equal length of beta");
  }
  double lval = 0.0;

  for (int i = 0; i < gl.nrow(); i++) {
    NumericMatrix::Row lvec = gl(i, _);
    lval += log_sum_exp_cpp(lvec + lq);
  }

  if (!lg) {
    lval = std::exp(lval);
  }

  return lval;
}


//' Sample genotypes from posteriors using genotype likelihoods and genotype priors
//'
//' @param gl The matrix of genotype log-likelihoods. Rows index individuals
//'     and columns index genotypes.
//' @param q The vector of genotype priors (not log-priors).
//'
//' @author David Gerard
//'
//' @noRd
// [[Rcpp::export]]
NumericVector sample_z(NumericMatrix& gl, NumericVector& q) {
  int n = gl.nrow();
  int ploidy = gl.ncol() - 1;
  IntegerVector ivec = Rcpp::seq(0, ploidy);
  NumericVector z(n);
  NumericVector probs(ploidy + 1);
  for (int j = 0; j < n; j++) {
    probs = gl(j, _);
    probs = probs + Rcpp::log(q);
    probs = Rcpp::exp(probs - log_sum_exp_cpp(probs));
    z(j) = Rcpp::sample(ivec, 1, false, probs)(0);
  }
  return z;
}

//' Gibbs sampler under random mating using genotype log-likelihoods.
//'
//' @param gl The matrix of genotype log-likelihoods. The columns index the
//'     dosages and the rows index the individuals. \code{gl[i,j]} is the
//'     genotype log-likelihood for individual i at dosage j. It is assumed
//'     that natural log is used.
//' @param alpha Vector of hyperparameters for the gamete frequencies.
//'     Should be length (x.length() - 1) / 2 + 1.
//' @param B The number of sampling iterations.
//' @param T The number of burn-in iterations.
//' @param more A logical. Should we also return posterior draws (\code{TRUE})
//'     or not (\code{FALSE}).
//' @param lg Should we return the log marginal likelihood (true) or not
//'     (false).
//'
//' @return A list with some or all of the following elements
//' \itemize{
//'   \item{\code{mx}: The estimate of the marginal likelihood}
//'   \item{\code{p_tilde}: The value of p used to evaluate the posterior density}.
//'   \item{\code{p}: The samples of the gamete frequencies}
//'   \item{\code{z}: The samples of the individual genotypes}
//'   \item{\code{post}: The samples of the full conditionals of p_tilde.}
//' }
//'
//' @author David Gerard
//'
//' @examples
//' set.seed(1)
//' ploidy <- 8
//'
//' ## Simulate under the null
//' p <- stats::runif(ploidy / 2 + 1)
//' p <- p / sum(p)
//' q <- stats::convolve(p, rev(p), type = "open")
//' nvec <- c(stats::rmultinom(n = 1, size = 100, prob = q))
//' gl <- simgl(nvec)
//'
//' gibbs_gl(gl = gl, alpha = rep(1, ploidy / 2 + 1), lg = TRUE)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List gibbs_gl(Rcpp::NumericMatrix& gl,
                    Rcpp::NumericVector alpha,
                    int B = 10000,
                    int T = 100,
                    bool more = false,
                    bool lg = false) {
  int ploidy = gl.ncol() - 1;
  int n = gl.nrow();
  if (alpha.length() != (ploidy / 2 + 1)) {
    Rcpp::stop("alpha should be the same length as ploidy / 2 + 1");
  }

  NumericVector y(alpha.length()); // latent gamete counts
  NumericVector x(ploidy + 1); // latent genotype counts
  IntegerVector z(n); // latent genotypes
  NumericVector p = rdirichlet1(alpha); // gamete frequencies
  NumericVector q = conv_cpp(p, p);
  NumericVector p_tilde = p;
  NumericVector q_tilde = conv_cpp(p_tilde, p_tilde);
  double logpost = ddirichlet(p, alpha, true) + plq(gl, q, true);
  double logpihat = R_NegInf; // estimate of posterior

  // build more output ----
  int nsamp;
  if (more) {
    nsamp = B;
  } else {
    nsamp = 0;
  }
  NumericMatrix pmat(nsamp, ploidy / 2 + 1);
  NumericMatrix zmat(nsamp, n);
  NumericVector postvec(nsamp);

  for (int i = 0; i < T + B; i++) {
    // sample z
    z = sample_z(gl, q);

    // calculate x
    x.fill(0.0);
    for (int j = 0; j < n; j++) {
      x(z(j)) = x(z(j)) + 1.0;
    }

    // sample y
    y = as<NumericVector>(samp_gametes(x, p));

    // sample p
    p = rdirichlet1(y + alpha);

    // calculate q
    q = conv_cpp(p, p);

    // Check if we are in the burnin
    if (i < T) {
      double lp2 = ddirichlet(p, alpha, true) + plq(gl, q, true);
      if (lp2 > logpost) {
        logpost = lp2;
        p_tilde = p;
        q_tilde = q;
      }
    } else {
      double ptilde_post = ddirichlet(p_tilde, y + alpha, true);
      logpihat = log_sum_exp_2_cpp(logpihat, ptilde_post);

      // include if more
      if (more) {
        NumericMatrix::Row crow = pmat(i - T, _);
        crow = p;

        NumericMatrix::Row czrow = zmat(i - T, _);
        czrow = z;

        postvec(i - T) = ptilde_post;
      }
    }
  }
  logpihat -= std::log((double)B);

  Rcpp::List retlist;
  if (lg) {
    retlist["mx"] = logpost - logpihat;
  } else {
    retlist["mx"] = exp(logpost - logpihat);
  }

  // include if more
  if (more) {
    retlist["p_tilde"] = p_tilde;
    retlist["p"] = pmat;
    retlist["z"] = zmat;
    retlist["post"] = postvec;
  }

  return retlist;
}


//' Gibbs sampler under the alternative of non-random mating using genotype
//' log-likelihoods.
//'
//' @inheritParams gibbs_gl
//' @param beta The concentration hyperparameter for the genotype frequencies.
//'
//' @return A list with some or all of the following elements
//' \itemize{
//'   \item{\code{mx}: The estimate of the marginal likelihood}
//' }
//'
//' @author David Gerard
//'
//' @examples
//' set.seed(1)
//' ploidy <- 8
//'
//' ## Simulate under the alternative
//' q <- stats::runif(ploidy + 1)
//' q <- q / sum(q)
//' nvec <- c(stats::rmultinom(n = 1, size = 100, prob = q))
//' gl <- simgl(nvec)
//'
//' gibbs_gl_alt(gl = gl, beta = rep(1, ploidy + 1), lg = TRUE)
//'
//' @export
// [[Rcpp::export]]
Rcpp::List gibbs_gl_alt(Rcpp::NumericMatrix& gl,
                        Rcpp::NumericVector beta,
                        int B = 10000,
                        int T = 100,
                        bool more = false,
                        bool lg = false) {
  int ploidy = gl.ncol() - 1;
  int n = gl.nrow();
  if (beta.length() != (ploidy + 1)) {
    Rcpp::stop("beta should be the same length as ploidy + 1");
  }

  NumericVector x(ploidy + 1); // latent genotype counts
  IntegerVector z(n); // latent genotypes
  NumericVector q = rdirichlet1(beta); // genotype frequencies
  NumericVector q_tilde = q;
  double logpost = ddirichlet(q, beta, true) + plq(gl, q, true);
  double logpihat = R_NegInf; // estimate of posterior

  // build more output ----
  int nsamp;
  if (more) {
    nsamp = B;
  } else {
    nsamp = 0;
  }
  NumericMatrix qmat(nsamp, ploidy + 1);
  NumericMatrix zmat(nsamp, n);
  NumericVector postvec(nsamp);
  NumericMatrix xmat(nsamp, ploidy + 1);

  for (int i = 0; i < T + B; i++) {
    // sample z
    z = sample_z(gl, q);

    // calculate x
    x.fill(0.0);
    for (int j = 0; j < n; j++) {
      x(z(j)) = x(z(j)) + 1.0;
    }

    // Sample q
    q = rdirichlet1(x + beta);

    // Check if we are in the burnin
    if (i < T) {
      double lp2 = ddirichlet(q, beta, true) + plq(gl, q, true);
      if (lp2 > logpost) {
        logpost = lp2;
        q_tilde = q;
      }
    } else {
      double qtilde_post = ddirichlet(q_tilde, x + beta, true);
      logpihat = log_sum_exp_2_cpp(logpihat, qtilde_post);

      // include if more
      if (more) {
        NumericMatrix::Row crow = qmat(i - T, _);
        crow = q;

        NumericMatrix::Row czrow = zmat(i - T, _);
        czrow = z;

        NumericMatrix::Row cxrow = xmat(i - T, _);
        cxrow = x;

        postvec(i - T) = qtilde_post;
      }
    }
  }
  logpihat -= std::log((double)B);

  Rcpp::List retlist;
  if (lg) {
    retlist["mx"] = logpost - logpihat;
  } else {
    retlist["mx"] = exp(logpost - logpihat);
  }

  // include if more
  if (more) {
    retlist["q_tilde"] = q_tilde;
    retlist["q"] = qmat;
    retlist["z"] = zmat;
    retlist["x"] = xmat;
    retlist["post"] = postvec;
  }

  return retlist;
}

