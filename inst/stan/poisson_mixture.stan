
data {
  int<lower=1> K;          // number of mixture components
  int<lower=1> N;          // number of data points
  int y[N];                // observations
  real<lower=0> e0;        // prior - number of components
  real<lower=0> a0;        // prior - number of components
  real<lower=0> A0;        // prior - number of components
  real l0;                 // prior - lambda
  real<lower=0> L0;        // prior - lambda
}

parameters {
  simplex[K] theta;                       // mixing proportions
  positive_ordered[K] lambda;             // poisson parameters
  vector<lower=0>[(e0>0) ? 0 : 1] alpha;  // parameter mixing proportions
}

model {
  vector[K] log_theta = log(theta);  
  
  if(e0>0){
    theta ~ dirichlet(rep_vector(e0, K));
  } else {
    alpha[1] ~ gamma(a0, A0);
    theta ~ dirichlet(rep_vector(alpha[1], K));
  }
  
  lambda ~ gamma(l0, L0);
  
  for (n in 1:N) {
    vector[K] lps = log_theta;
    for (k in 1:K) {
      lps[k] += poisson_lpmf(y[n]| lambda[k]);
    }
    target += log_sum_exp(lps);
  }
}
