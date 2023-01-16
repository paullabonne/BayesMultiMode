// https://groups.google.com/g/stan-users/c/W4hS1mtUtxg?pli=1
// the prior on the degrees of freedom comes from :
// Juárez and Steel (2010) (Model-based clustering of non-Gaussian panel data based on skew-t distributions. Journal of Business & Economic Statistics 28, 52–66.)

data {
  int<lower=1> K;          // number of mixture components
  int<lower=1> N;          // number of data points
  real y[N];               // observations
  real<lower=0> e0;                 // prior mean variance
  real<lower=0> a0;                 // prior number of components
  real<lower=0> A0;                 // prior number of components
  real b0;                 // prior mean mean
  real<lower=0> B0;                 // prior mean variance
  real<lower=0> c0;                 // prior variance mean
  real<lower=0> g0;                 // prior variance mean
  real<lower=0> G0;                 // prior variance mean
  real<lower=0> n0;                 // prior variance mean
  real<lower=0> N0;                 // prior variance mean
}

parameters {
  simplex[K] theta;          // mixing proportions
  ordered[K] mu;             // locations of mixture components
  vector<lower=0>[K] sigmaSQ;  // scales of mixture components
  vector<lower=1>[K] nu;  // scales of mixture components
  vector<lower=0>[(e0>0) ? 0 : 1] alpha;                // parameter mixing proportions
  vector<lower=0>[(G0>0) ? K : 0] C0; // see https://discourse.mc-stan.org/t/if-else-statement-inside-parameter-block/13937/3
}

transformed parameters {
  vector<lower=0>[K] sigma;  // scales of mixture components
  
  for (i in 1:K) {
    sigma[i] = sqrt(sigmaSQ[i]);
  }
}

model {
      vector[K] log_theta;

  if (G0>0){
      C0 ~ gamma(g0, G0);
      sigmaSQ ~ inv_gamma(c0, C0);
  } else {
      sigmaSQ ~ inv_gamma(c0, g0);
  }

  mu ~ normal(b0, B0);

  nu ~ gamma(n0, N0);

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
      lps[k] += student_t_lpdf(y[n] | nu[k], mu[k], sigma[k]);
    }
    target += log_sum_exp(lps);
  }
}
