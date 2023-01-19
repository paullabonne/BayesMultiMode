// This code treats kappas like the mixing parameters :
// they are marginalised and a prior is put on the resulting proportions.
// Stan cannot sample discrete parameters
// It is possible to use a uniform prior on the kappas when marginalising but that does not work well.

// https://mc-stan.org/docs/stan-users-guide/change-point.html for an example of marginalisation with stan
// Gives example for two models including a Gaussian mixture. See also https://arxiv.org/abs/2204.06313 

data {
  int<lower=1> K;          // number of mixture components
  int<lower=1> N;          // number of data points
  int y[N];                // observations
  real<lower=0> e0;                 // prior variance mean 
  real<lower=0> a0;         // prior number of components
  real<lower=0> A0;         // prior number of components
  real l0;                 // prior mean mean 1.1
  real<lower=0> L0;        // prior mean variance 1.1/mean(yy)
  real<lower=0> e0_kappa;                 // prior variance mean 
  real<lower=0> d0;         // prior kappa
  real<lower=0> D0;         //prior kappa
}

transformed data {
  int max_y = max(y);
}

parameters {
  simplex[K] theta;          // mixing proportions
  positive_ordered[K] lambda;             
  vector<lower=0>[(e0>0) ? 0 : 1] alpha;                // parameter mixing proportions
  simplex[max_y+1] kappa[K];             // locations of mixture components
  vector<lower=0>[(e0_kappa>0) ? 0 : 1] alpha_kappa;                // parameter mixing proportions
}

// This block is necessary to integrate, or marginalise, kappa. This is the only way to deal with discrete parameters in 
// Stan; it is not possible to sample discrete parameters. 
transformed parameters {
  vector[K] log_theta = log(theta);     // cache calculations on log mixture proportions.
  vector[N] ll;
  
  for (n in 1:N) {
    // for each observation
    
    vector[K] lps;
    
    for (k in 1:K) {
      // for each mixture component
      vector[max_y+1] log_kappa = log(kappa[k]);
      
      // vector[range] log_kappa = log(kappa[k]); // cache calculations on log kappas proportions.
      vector[max_y+1] lp;
      
      for (j in 0:max_y) {
        // for each kappa parameter
        
        if((y[n] - j) > -1){
          lp[j+1] = log_theta[k] + log_kappa[j+1] + poisson_lpmf(y[n] - j | lambda[k]);  //log(kappa[1+kap]) 
        } else {
          lp[j+1] = log(0);
        }
        
      }
      
      lps[k] = log_sum_exp(lp);
    }
    
    ll[n] = log_sum_exp(lps);
  }
}

model {
  lambda ~ gamma(l0, L0);
  
  if(e0>0){
    theta ~ dirichlet(rep_vector(e0, K));
  } else {
    alpha[1] ~ gamma(a0, A0);
    theta ~ dirichlet(rep_vector(alpha[1], K));
  }
  
  if(e0_kappa>0){
    for (i in 1:K) {
      kappa[i] ~ dirichlet(rep_vector(e0_kappa, max_y+1));
    } 
  } else {
    alpha_kappa[1] ~ gamma(d0, D0);
    for (i in 1:K) {
      kappa[i] ~ dirichlet(rep_vector(alpha_kappa[1], max_y+1));
    } 
  }
  
  target += sum(ll);
}
