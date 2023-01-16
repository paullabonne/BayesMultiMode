// see https://mc-stan.org/docs/stan-users-guide/summing-out-the-responsibility-parameter.html
// mixture of normal priors from https://link.springer.com/article/10.1007/s11634-021-00461-8

data {
  int<lower=1> K;          // number of mixture components
  int<lower=1> N;          // number of data points
  real y[N];               // observations
  real<lower=0> e0;
  real<lower=0> a0;                 // prior number of components
  real<lower=0> A0;                 // prior number of components
  real b0;                 // prior mean mean
  real<lower=0> B0;                 // prior mean variance
  real<lower=0> c0;                 // prior variance mean
  real<lower=0> g0;                 // prior variance mean
  real<lower=0> G0;                 // prior variance mean
}

parameters {
  simplex[K] theta;          // mixing proportions
  ordered[K] mu;             // locations of mixture components
  vector<lower=0>[K] sigma;  // scales of mixture components
  vector<lower=0>[(e0>0) ? 0 : 1] alpha;                // parameter mixing proportions
  vector<lower=0>[(G0>0) ? K : 0] C0; // see https://discourse.mc-stan.org/t/if-else-statement-inside-parameter-block/13937/3
}

model {
  vector[K] log_theta;

  if (G0>0){
    C0 ~ gamma(g0, G0);
    sigma^2 ~ inv_gamma(c0, C0);
  } else {
    sigma^2 ~ inv_gamma(c0, g0);
  }
  mu ~ normal(b0, B0);

  if(e0>0){
    theta ~ dirichlet(rep_vector(e0, K));
  } else {
    alpha[1] ~ gamma(a0, A0);
    theta ~ dirichlet(rep_vector(alpha[1], K));
  }

  log_theta = log(theta);  // cache log calculation

  for (n in 1:N) {
    vector[K] lps = log_theta;
    for (k in 1:K) {
      lps[k] += normal_lpdf(y[n] | mu[k], sigma[k]);
    }
    target += log_sum_exp(lps);
  }
}
