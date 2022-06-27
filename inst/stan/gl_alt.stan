data {
  int<lower=0> N; // sample size
  int<lower=2> K; // ploidy
  matrix[N, K+1] gl; // genotype log-likelihoods
  vector[K + 1] beta; // concentration hyperparameters
}

parameters {
  simplex[K + 1] q; // genotype frequencies
}

model {
  target += dirichlet_lpdf(q |beta);
  for (i in 1:N) {
    vector[K+1] lq;
    for (j in 1:(K+1)) {
      lq[j] = gl[i, j] + log(q[j]);
    }
    target += log_sum_exp(lq);
  }
}

