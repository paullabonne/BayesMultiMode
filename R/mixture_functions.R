#' @importFrom sn dst
#' @importFrom sn dsn
#' @importFrom stats dnorm
#' 

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
skew_norm_mix <- function(x, p, xi, omega, alpha) {
  mixture = 0
  
  for (i in 1:length(p)) {
    mixture = mixture + p[i] * dsn(x,
                                   xi[i],
                                   omega[i],
                                   alpha[i])
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

# Mixture of pdf_func
pdf_func_mix <- function(x, pars, pdf_func) {
  pdf_func = match.fun(pdf_func) #solves NOTE "pdf_func is undefined"
  
  mixture = 0
  for (i in 1:nrow(pars)) {
    mixture = mixture + pars[i, 1] * pdf_func(x, pars[i, -1])
  }
  return(mixture)
}

# wrapper
dist_mixture <- function(x, dist, pars, pdf_func = NULL) {
  if (!is.null(pdf_func)) {
    output = pdf_func_mix(x, pars, pdf_func)
    
  } else {
    if (dist == "normal") {
      output = normal_mix(x, pars[, "eta"], pars[, "mu"], pars[, "sigma"])
    }
    
    if (dist == "skew_normal") {
      output = skew_norm_mix(x, pars[, "eta"], pars[, "xi"], pars[, "omega"], pars[, "alpha"])
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
      output = dnorm(x, pars["mu"], pars["sigma"])
    }
    
    if (dist == "skew_normal") {
      output = dsn(x, pars["xi"], pars["omega"], pars["alpha"])
    }
    
    if (dist == "shifted_poisson") {
      output = dpois(x - pars["kappa"], pars["lambda"])
    }
    
    if (dist == "poisson") {
      output = dpois(x, pars)
    }
  }
  
  return(output)
  
}