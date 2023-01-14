#' @importFrom sn dst
#' @importFrom sn dsn
#' @importFrom stats dnorm
#' 
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
normal_mix <- function(x, p, mu, sigma) {
  mixture = 0
  
  for (i in 1:length(p)) {
    mixture = mixture + p[i] * dnorm(x,
                                     mean = mu[i],
                                     sd = sigma[i])
  }
  
  return(mixture)
}

# Mixture of skew normals
skew_norm_mix <- function(x, p, mu, sigma, xi) {
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
student_mix <- function(x, p, mu, sigma, nu) {
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
skew_t_mix <- function(x, p, mu, sigma, xi, nu) {
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

# Mixture of shifted Poisson distributions
shift_pois_mix <- function(x, p, lambda, kappa) {
  mixture = 0
  for (i in 1:length(kappa)) {
    mixture = mixture + p[i] * dpois(x - kappa[i], lambda[i], log = FALSE)
  }
  return(mixture)
}

# Mixture 
pdf_func_mix <- function(x, pars) {
  mixture = 0
  for (i in 1:nrow(pars)) {
    mixture = mixture + pars[i, ] * pdf_func(x, pars[i, -1])
  }
  return(mixture)
}

# wrapper
dist_mixture <- function(x, dist, pars, pdf_func = NULL) {
 
  if (!is.null(pdf_func)) {
    output = pdf_func_mix(x, pars)
  } else {
    if (dist == "normal") {
      output = normal_mix(x, pars[, "theta"], pars[, "mu"], pars[, "sigma"])
    }
    
    if (dist == "student") {
      output = student_mix(x, pars[, "theta"], pars[, "mu"], pars[, "sigma"], pars[, "nu"])
    }
    
    if (dist == "skew_normal") {
      output = skew_norm_mix(x, pars[, "theta"], pars[, "mu"], pars[, "sigma"], pars[, "xi"])
    }
    
    if (dist == "shift_pois_mix") {
      output = shift_pois_mix(x, pars[, "theta"], pars[, "lambda"], pars[, "lambda"])
    }
  }
  
  return(output)
  
}

# wrapper
dist_pdf <- function(x, dist, pars, pdf_func = NULL) {
  
  if (!is.null(pdf_func)) {
    pdf_func = match.fun(pdf_func) #solves NOTE "pdf_func is undefined"
    
    dist_pdf = pdf_func(x, pars)
    
  } else {
    
    if (dist == "normal") {
      output = dnorm(x, pars[, "mu"], pars[, "sigma"])
    }
    
    if (dist == "student") {
      output = dst_vec(x, pars[, "mu"], pars[, "sigma"], pars[, "nu"])
    }
    
    if (dist == "skew_normal") {
      output = dsn(x, pars[, "mu"], pars[, "sigma"], pars[, "xi"])
    }
    
    if (dist == "shift_pois_mix") {
      output = dpois(x - pars[, "kappa"], pars[, "lambda"])
    }
  }
  
  return(output)
  
}