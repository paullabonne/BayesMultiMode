# Function to count the number of modes from an empirical pdf py evaluated at points y
# Author: N. Basturk
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
fn.sub.mixpois<-function(theta.draws_i, y, which.r, Khat = Khat){
  theta.draws_i = theta.draws_i[!is.na(theta.draws_i)]

  theta <- cbind(theta.draws_i[1:Khat],
                 theta.draws_i[(Khat+1):(2*Khat)],
                 theta.draws_i[(2*Khat+1):(3*Khat)])
  
  p <- theta[,1]
  kappa <- theta[,2]
  lambda <-  theta[,3]
  
  ### Getting individual component densities
  pdf.J = matrix(nrow=length(y),ncol=Khat) 
  for(j in 1:Khat){
    pdf.J[,j] = dpois((y-kappa[j]),lambda[j]) * p[j]
  }
  
  ### summing up to get the mixture
  py <- rowSums(pdf.J)
  
  ### Finding the modes
  out <- disc_modes(y,py, mode.sel = "leftmode")
  r <- num.modes <- out$num.modes
  if(which.r == 2){
    # mode locations
    r = rep(NA, Khat*2)
    r[1:length(out$y.peaks)] = out$y.peaks
  }
  if(which.r == 3){
    r = rep(NA, Khat*2)
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
