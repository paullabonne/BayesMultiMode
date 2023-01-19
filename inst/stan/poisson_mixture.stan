// see https://mc-stan.org/docs/stan-users-guide/summing-out-the-responsibility-parameter.html

data {
  int<lower=1> K;          // number of mixture components
  int<lower=1> N;          // number of data points
  int y[N];                // observations
  real<lower=0> e0;                 // prior variance mean 
  real<lower=0> a0;         // prior number of components
  real<lower=0> A0;         // prior number of components
  real l0;                 // prior mean mean 1.1
  real<lower=0> L0;        // prior mean variance 1.1/mean(yy)
}

parameters {
  simplex[K] theta;          // mixing proportions
  vector<lower=0>[K] lambda;            
  vector<lower=0>[(e0>0) ? 0 : 1] alpha;                // parameter mixing proportions
}

model {
    vector[K] log_theta = log(theta);  // cache log calculation
    
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
