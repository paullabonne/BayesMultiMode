#' @importFrom sn dst
#' @importFrom sn dsn
#' @importFrom stats dnorm
#' 
# Vectorise the dst function for vector of nus
dst_vec <- function(x, xi, omega, alpha = NULL, nu){
  if(is.null(alpha)) {
    alpha = rep(0, length(xi))
  }
  n = length(xi)
  output = rep(NA, n)
  
  for (i in 1:n) {
    output[i] = dst(x, xi = xi[i], omega = omega[i], alpha = alpha[i], nu = nu[i])
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

# Mixture of shifted Poisson
pois_mix <- function(x, p, lambda, kappa) {
  mixture = 0
  for (i in 1:length(kappa)) {
    mixture = mixture + p[i] * dpois(x, lambda[i], log = FALSE)
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

# Vectorised version of pdf_func
pdf_func_vec <- function(pdf_func) {
  pdf_func = match.fun(pdf_func) #solves NOTE "pdf_func is undefined"
  
  func <- function(x, pars){
    n = nrow(pars)
    output = rep(NA, n)
    
    for (i in 1:n) {
      output[i] = pdf_func(x, pars[i, ])
    }
    
    return(output)
  }
  
  return(func)
}

# Mixture of pdf_func
pdf_func_mix <- function(x, pars, pdf_func) {
  pdf_func = match.fun(pdf_func) #solves NOTE "pdf_func is undefined"
  
  mixture = 0
  for (i in 1:nrow(pars)) {
    mixture = mixture + pars[i, 1] * pdf_func(x, pars[i, -1, drop = F])
  }
  return(mixture)
}

# wrapper
dist_mixture <- function(x, dist, pars, pdf_func = NULL) {
  if (!is.null(pdf_func)) {
    output = pdf_func_mix(x, pars, pdf_func)
    
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
    
    if (dist == "skew_t") {
      output = skew_t_mix(x, pars[, "theta"], pars[, "mu"], pars[, "sigma"], pars[, "xi"], pars[, "nu"])
    }
    
    if (dist == "pois_mix") {
      output = pois_mix(x, pars[, "theta"], pars[, "lambda"])
    }
    
    if (dist == "shift_pois_mix") {
      output = shift_pois_mix(x, pars[, "theta"], pars[, "lambda"], pars[, "kappa"])
    }
  }
  
  return(output)
  
}

# wrapper
dist_pdf <- function(x, dist, pars, pdf_func = NULL) {
  
  if (!is.null(pdf_func)) {
    pdf_func = match.fun(pdf_func) #solves NOTE "pdf_func is undefined"
    output = pdf_func(x, pars)
    
  } else {
    if (dist == "normal") {
      output = dnorm(x, pars[, "mu"], pars[, "sigma"])
    }
    
    if (dist == "student") {
      output = dst_vec(x, xi = pars[, "mu"], omega = pars[, "sigma"], nu = pars[, "nu"])
    }
    
    if (dist == "skew_normal") {
      output = dsn(x, pars[, "mu"], pars[, "sigma"], pars[, "xi"])
    }
    
    if (dist == "skew_t") {
      output = dst_vec(x, xi = pars[, "mu"], omega = pars[, "sigma"],
                       alpha = pars[, "xi"], nu = pars[, "nu"])
    }
    
    if (dist == "shifted_poisson_bis") {
      output = dpois(x - pars[, "kappa"], pars[, "lambda"])
    }
    
    if (dist == "poisson") {
      output = dpois(x, pars[, "lambda"])
    }
  }
  
  return(output)
  
}