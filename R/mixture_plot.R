#' Plot mixture
#' 
#' @importFrom posterior as_draws_matrix
#' @import ggplot2
#' 
#' @export

mixture_plot <- function(y, K, fit, dist, max_size = 200, col = "magenta", tol = 1e-4) {
  mcmc_output = as_draws_matrix(fit)
  
  ## plot the data
  g = ggplot(data.frame(y = y), aes(y)) +
    geom_histogram(aes(y = ..density..),
                   colour = "white",
                   bins = 70)
  
  ## plot the mixture for each draw
  for (i in sample(1:(min(nrow(mcmc_output), max_size)))) {
    pars = mcmc_output[i, grep("theta", colnames(mcmc_output))]
    pars = rbind(pars, mcmc_output[i, grep("mu", colnames(mcmc_output))])
    pars = rbind(pars, sqrt(mcmc_output[i, grep("sigma", colnames(mcmc_output))]))
    
    if (dist == "student") {
      pars = rbind(pars, mcmc_output[i, grep("nu", colnames(mcmc_output))])
    }
    
    if (dist == "skew_normal") {
      pars = rbind(pars, mcmc_output[i, grep("xi", colnames(mcmc_output))])
    }
    
    if (dist == "skew_t") {
      pars = rbind(pars, mcmc_output[i, grep("xi", colnames(mcmc_output))])
      pars = rbind(pars, mcmc_output[i, grep("nu", colnames(mcmc_output))])
    }
    
    pars = t(pars)
    pars = pars[pars[, 1]>tol, ]
    
    g = g +
      geom_function(fun = dist_mixture,
                    args = list(dist = dist,
                                pars = pars),
                    alpha = 0.1,
                    colour = col)
  }
  
  g
}
