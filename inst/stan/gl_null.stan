data {
  int<lower=0> N; // sample size
  int<lower=2> K; // ploidy
  int<lower=2> khalf; // K/2+1 so stan does not complain about integer division
  matrix[N, K+1] gl; // genotype log-likelihoods
  vector[khalf] alpha; // concentration hyperparameters
}

parameters {
  simplex[khalf] p; // gamete frequencies
}

transformed parameters {
  // convolve p to get q.
  simplex[K+1] q; // genotype frequencies
  for (k in 1:(K+1)) {
    int iup = min(k - 1, khalf - 1);
    int ilo = max(0, k - khalf);
    q[k] = 0.0;
    for (i in ilo:iup) {
      q[k] += p[i + 1] * p[k - i];
    }
  }
}

model {
  target += dirichlet_lpdf(p | alpha);
  for (i in 1:N) {
    vector[K+1] lq;
    for (j in 1:(K+1)) {
      lq[j] = gl[i, j] + log(q[j]);
    }
    target += log_sum_exp(lq);
  }
}
