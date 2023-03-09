#' Mode-finding EM algorithm (MEM)
#' 
#' From Li, Jia, Surajit Ray, and Bruce G. Lindsay.
#' "A nonparametric statistical approach to clustering via mode identification."
#' Journal of Machine Learning Research 8 (2007): 1687-1723.
#' 
#' @param mcmc a vector of estimated mixture parameters.
#' @param dist String indicating the distribution of the mixture components.
#' Currently supports "normal", "student" and "skew_normal".
#' @param data Numeric vector of observations.
#' @param pars_names Names of the variables mcmc draws variables
#' @param pdf_func Pdf or pmf of the mixture components associated with the mcmc draws (if estimation not in-house)
#' @param tol_x Tolerance for distance in-between modes. Default is sd(data)/10. If two modes are closer than tol_x, only the first estimated mode is kept.
#' @param show_plot Show the data and estimated modes.
#' 
#' @return A vector estimated modes.
#' 
#' @importFrom sn dst
#' @importFrom sn dsn
#' @importFrom graphics abline curve
#' @importFrom stats optim na.omit var
#' @importFrom assertthat assert_that
#' @importFrom assertthat is.string

#' 
#' @export

MEM <- function(mcmc, dist = "NA", pars_names, data, pdf_func = NULL, tol_x = sd(data)/10, show_plot = FALSE) {
  ## input checks
  fail = "inputs to the Mode-finding EM algorithm are corrupted"
  assert_that(is.vector(mcmc) & length(mcmc) >= 3,
              msg = paste0("mcmc should be a vector of length >= 3", fail))
  assert_that(is.string(dist),
              msg = paste0("dist should be a string", fail))
  assert_that(is.vector(data) & length(data) > 0,
              msg = "data should be a vector of length > 0")
  assert_that(is.vector(tol_x) & tol_x > 0, msg = paste0("tol_x should be a positive scalar", fail))
  assert_that(is.logical(show_plot), msg = paste0("show_plot should be TRUE or FALSE", fail))
  ##
  
  ##
  names_mcmc = str_to_lower(names(mcmc))
  names_mcmc = str_extract(names_mcmc, "[a-z]+")
  names_mcmc = unique(names_mcmc)
  
  assert_that(sum(pars_names %in% names_mcmc)==length(pars_names),
              msg = paste0("the name of the parameters provided by pars_names and those of the mcmc vector do not match; ", fail))
  
  if (dist %in% c("skew_normal")) {
    assert_that(length(pars_names) == 4,
                msg = paste0("the number of elements in pars_names does not match with dist; ", fail)) 
  }
  ##
  
  pars = c()
  for (i in 1:length(pars_names)) {
    pars = cbind(pars, mcmc[grep(pars_names[i], names(mcmc))])
  }
  
  colnames(pars) <- pars_names

  est_mode = rep(NA, nrow(pars))

  nK = nrow(pars)
  post_prob = rep(NA, nK)
  
  # vectorising function
  if (!is.null(pdf_func)) {
    pdf_func <- pdf_func_vec(pdf_func)
  }

  # remove empty components (a feature of some MCMC methods)
  pars = na.omit(pars)
  
  for (j in 1:nK) {
    x = pars[j,2]
    
    delta = 1
    
    while (delta > 1e-8) {
      # E-step
      f_mix = dist_mixture(x, dist, pars, pdf_func)
      
      for (k in 1:nK){ 
        post_prob[k] = pars[k, 1] * dist_pdf(x, dist, pars[k, -1, drop = F], pdf_func)/f_mix 
      }
     
      # M-step
      Min = optim(par = x, Q_func, method = "L-BFGS-B",
                  dist = dist, 
                  post_prob = post_prob,
                  pars = pars[, -1, drop = F],
                  pdf_func = pdf_func,
                  control = list(fnscale = -1))
      
      x1 = Min$par
      
      # check convergence and increment
      delta = abs(x - x1)
      x = x1
    }
 
    ## check that the mode is not too close to other modes
    not_duplicate = TRUE
    
    if (any(is.finite(est_mode))) {
      diff = abs(x-est_mode)
      diff = diff[!is.na(diff)]
      if (any(diff<tol_x)) {
        not_duplicate = FALSE
      }
    }
    
    if (x <= max(data) & x >= min(data) & not_duplicate){
      est_mode[j] = x
    }
  }
  
  if (show_plot) {
    curve(dist_mixture(x, dist, pars), from = min(data), to =  max(data))
    for (x in est_mode) {
      abline(v = x) 
    } 
  }
  
  return(est_mode)
}

#' @keywords internal
Q_func = function(x, dist, post_prob, pars, pdf_func){
  
  pdf = dist_pdf(x, dist, pars, pdf_func = pdf_func)
  pdf[pdf==0] = 1e-10 #otherwise the log operation below through infs
  Q = sum(post_prob * log(pdf))
  
  if(is.na(Q)|!is.finite(Q)){
    Q = -1e6
  }
  
  return(Q)
}