# Function to count the number of modes from an empirical pdf py evaluated at points y
# Authors: N. Basturk, J. Cross, P. Labonne
# Created: 20.03.2015
# inputs: 
#	y  : [vector (N)] of evaluation points
#	py : [vector (N)] of pdf at each evaluation point
# 	mode.sel: [string] optional argument from following options:  'leftmode', 'midmode', 'rightmode'. Defines how the mode is calculated in case 'py' gives adjecent modes.  'leftmode' gives the smallest value, 'rightmode' gives the largest value,  'midmode' gives the mid point (within y) in adjecent modes. Default value is  'leftmode'.
# output: list containing following elements
#	y.peaks	 : [vector (M)] vector of modes
#	py.peaks : [vector (M)] vector of pdf values for each mode
#	num.modes: [integer > 0] number of modes, M.

disc_modes <- function(y, py, mode.sel = NULL){
  mode.sel.options <- c('leftmode', 'midmode', 'rightmode')
  
  if(!is.vector(y)) stop("y should be a vector")
  if(any(y!=sort(y))) stop("y should be ordered")
  if(length(y)!=length(py)) stop("length of y and py should match")
  if(is.null(mode.sel)) mode.sel <- 'leftmode'
  if(all(mode.sel != mode.sel.options))	stop("Mode selection is not properly defined for 'count_modes.R'")
  
  n <- length(y)
  py.left <- c(-Inf,py[1:(n-1)]) # density on the left of every point
  py.right <- c(py[2:n],-Inf)    # density on the right of every point
  
  cond.left <- which(py.left < py)     # condition for increasing values from left
  
  cond.r <- c()# condition for peak from left
  cond.l <- c()# condition for peak from right
  for(i in 1: length(cond.left)){
    i.left <- cond.left[i]
    i.right <-  i.left		
    while(py[i.right] == py.right[i.right]){  # remove non-modes (steps) if any subsequent 
      i.right = i.right + 1
    }
    if(py[i.right] > py.right[i.right]){
      cond.r <- c(cond.r, i.right)
      cond.l <- c(cond.l, i.left)
    }
  }
  
  i.mode <- which(mode.sel  == mode.sel.options)	
  fn.sub <- function(x,...){
    tmp <- cbind(x[cond.l], x[cond.r])
    ret <- rbind(tmp[,1], rowMeans(tmp), tmp[,2])
    return(ret)
  }
  y.peaks <- fn.sub(y)[i.mode,]
  py.peaks <- fn.sub(py)[i.mode,]
  num.modes <- length(y.peaks)
  
  return(list(y.peaks=y.peaks,py.peaks=py.peaks,num.modes=num.modes))
}


#' @keywords internal
### COUNT NUMBER OF MODES
fn.sub.mixpois <- function(mcmc, y, which.r, pars_names, tol_p, dist, pdf_func = NULL){
  
  if (!is.null(pdf_func)) {
    pdf_func <- pdf_func_vec(pdf_func)
  }
  
  ## input checks
  fail = "inputs to the discrete mode-finding algorithm are corrupted"
  assert_that(is.vector(mcmc) & length(mcmc) >= 3,
              msg = paste0("mcmc should be a vector of length >= 3", fail))
  assert_that(is.string(dist),
              msg = paste0("dist should be a string", fail))
  assert_that(is.vector(y) & length(y) > 0,
              msg = "y should be a vector of length > 0")
  assert_that(is.vector(tol_p) & tol_p > 0, msg = paste0("tol_p should be a positive scalar", fail))
  ##
  
  ##
  names_mcmc = str_to_lower(names(mcmc))
  names_mcmc = str_extract(names_mcmc, "[a-z]+")
  names_mcmc = unique(names_mcmc)
  
  assert_that(sum(pars_names %in% names_mcmc)==length(pars_names),
              msg = paste0("the name of the parameters provided by pars_names and those of the mcmc vector do not match; ", fail))
  
  if (dist %in% c("shifted_poisson")) {
    assert_that(length(pars_names) == 3,
                msg = paste0("the number of elements in pars_names does not match with dist; ", fail))
  }
  if (dist %in% c("poisson")) {
    assert_that(length(pars_names) == 2,
                msg = paste0("the number of elements in pars_names does not match with dist; ", fail)) 
  }
  ##
  
  pars = c()
  for (i in 1:length(pars_names)) {
    pars = cbind(pars, mcmc[grep(pars_names[i], names(mcmc))])
  }
  
  colnames(pars) <- pars_names
  ##
  
  # mcmc = mcmc[!is.na(mcmc)]
  Khat = nrow(pars)
  
  # keep = which(pars[,1] > tol_p)
  # pars = pars[keep, , drop = F]
  
  ### Getting individual component densities
  pdf = matrix(0, nrow=length(y), ncol=Khat) 
  for(k in 1:nrow(pars)){
    pdf[,k] = pars[k,1] * dist_pdf(y, dist, pars[k, -1, drop = F], pdf_func)
  }
  
  ### summing up to get the mixture
  py <- rowSums(pdf, na.rm = T)
  
  ### Finding the modes
  out <- disc_modes(y, py, mode.sel = "leftmode")
  r <- num.modes <- out$num.modes
  if(which.r == 2){
    # mode locations
    r = rep(NA, length(y))
    r[1:length(out$y.peaks)] = out$y.peaks
  }
  if(which.r == 3){
    r = rep(NA, length(y))
    #density at modes
    r[1:length(out$py.peaks)] = out$py.peaks
  }
  if(which.r == 4){
    r = rep(0,length(y))
    
    # account for flat modes
    out_right <- disc_modes(y,py, mode.sel = "rightmode")
    for(i in 1:length(out$y.peaks)){
      r[y==out$y.peaks[i]] = 1
      
      if(out$y.peaks[i]!=out_right$y.peaks[i]){
        r[y==out_right$y.peaks[i]] = 1
        
        diff = out_right$y.peaks[i]-out$y.peaks[i]
        if(diff>1){
          for(j in 1:(diff-1))
            r[y==out$y.peaks[i]+j] = 1
        }
      }
    }
    #
  }
  return(r)
}
