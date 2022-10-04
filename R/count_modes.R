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

count_modes <- function(y, py, mode.sel = NULL){
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
