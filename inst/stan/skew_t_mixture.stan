// https://groups.google.com/g/stan-users/c/W4hS1mtUtxg?pli=1
// the prior on the degrees of freedom comes from :
// Juárez and Steel (2010) (Model-based clustering of non-Gaussian panel data based on skew-t distributions. Journal of Business & Economic Statistics 28, 52–66.)


functions {
  // page 102-103 of The Skew-Normal and Related Families
  real skew_t_lpdf(real y, real mu, real sigma, real xi, real nu) {
    real z = (y - mu) / sigma;
    real w = sqrt((nu + 1) / (nu + z^2));
    real log_t_f = student_t_lpdf(z | nu, 0, 1);
    real log_t_F = student_t_lcdf(xi*z*w | nu+1, 0, 1);

    return log2() + -log(sigma) + log_t_f + log_t_F;
  }
}

data {
  int<lower=1> K;          // number of mixture components
  int<lower=1> N;          // number of data points
  real y[N];               // observations
  real<lower=0> e0;         // prior number of components
  real<lower=0> a0;         // prior number of components
  real<lower=0> A0;         // prior number of components
  real b0;                 // prior mean mean
  real<lower=0> B0;        // prior mean variance
  real<lower=0> c0;        // prior variance mean
  real<lower=0> g0;     // prior variance mean
  real<lower=0> G0;     // prior variance mean
  real<lower=0> h0;     // prior variance mean
  real<lower=0> H0;     // prior variance mean
  real<lower=0> n0;     // prior variance mean
  real<lower=0> N0;     // prior variance mean
}

parameters {
  simplex[K] theta;          // mixing proportions
  ordered[K] mu;             // locations of mixture components ordered[K]
  vector<lower=0>[K] sigma;  // scales of mixture components
  real xi[K];
  vector<lower=1>[K] nu;
  vector<lower=0>[(G0>0) ? K : 0] C0; // see https://discourse.mc-stan.org/t/if-else-statement-inside-parameter-block/13937/3
  vector<lower=0>[(e0>0) ? 0 : 1] alpha;                // parameter mixing proportions
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

  xi ~ normal(h0, H0);

  if(e0>0){
    theta ~ dirichlet(rep_vector(e0, K));
  } else {
    alpha[1] ~ gamma(a0, A0);
    theta ~ dirichlet(rep_vector(alpha[1], K));
  }
    log_theta = log(theta);  // cache log calculation

  nu ~ gamma(n0, N0);

  for (n in 1:N) {
    vector[K] lps = log_theta;
    for (k in 1:K) {
      lps[k] += skew_t_lpdf(y[n] | mu[k], sigma[k], xi[k], nu[k]);
    }
    target += log_sum_exp(lps);
  }
}
