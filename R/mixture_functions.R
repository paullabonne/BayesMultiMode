#' @importFrom sn dst
#' @importFrom sn dsn

# Vectorise the dst function for vector of nus
dst_vec <- function(x, xi, omega, nu){
  n = length(xi)
  output = rep(NA, n)
  
  for (i in 1:n) {
    output[i] = dst(x, xi = xi[i], omega = omega[i], nu = nu[i])
  }
  
  return(output)
}

# Mixture of normals
Gaussian_mixture <- function(x, p, mu, sigma) {
  mixture = 0
  
  for (i in 1:length(p)) {
    mixture = mixture + p[i] * dnorm(x,
                                     mean = mu[i],
                                     sd = sigma[i])
  }
  
  return(mixture)
}

# Mixture of skew normals
SN_mixture <- function(x, p, mu, sigma, xi) {
  mixture = 0
  
  for (i in 1:length(p)) {
    mixture = mixture + p[i] * dsn(x,
                                   xi = mu[i],
                                   omega = sigma[i],
                                   alpha = xi[i])
  }
  
  return(mixture)
}

# Mixture of student's t-distributions
ST_mixture <- function(x, p, mu, sigma, nu) {
  mixture = 0
  
  for (i in 1:length(p)) {
    mixture = mixture + p[i] * dst(x,
                                   xi = mu[i],
                                   omega = sigma[i],
                                   nu = nu[i])
  }
  
  return(mixture)
}


# currently not used
# Mixture of skewed student's t-distributions
SkT_mixture <- function(x, p, mu, sigma, xi, nu) {
  mixture = 0
  
  for (i in 1:length(p)) {
    mixture = mixture + p[i] * dst(x,
                                   xi = mu[i],
                                   omega = sigma[i],
                                   alpha = xi[i],
                                   nu = nu[i]
    )
  }
  
  return(mixture)
}

# wrapper
dist_mixture <- function(x, dist, pars, K) {
  
  if (dist == "gaussian") {
    output = Gaussian_mixture(x, pars[, 1], pars[, 2], pars[, 3])
  }
  
  if (dist == "student") {
    output = ST_mixture(x, pars[, 1], pars[, 2], pars[, 3], pars[, 4])
  }
  
  if (dist == "skew_normal") {
    output = SN_mixture(x, pars[, 1], pars[, 2], pars[, 3], pars[, 4])
  }
  
  if (dist == "skew_t") {
    output = SkT_mixture(x, pars[, 1], pars[, 2], pars[, 3], pars[, 4], pars[, 5])
  }
  
  return(output)
  
}