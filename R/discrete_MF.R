#' Mode-finding algorithm for mixture of discrete distributions
#' 
#' @param mcmc Vector of estimated mixture parameters.
#' @param data Vector of observations used for estimating the mixture.
#' @param pars_names Names of of the mcmc variables.
#' @param dist String indicating the distribution of the mixture components.
#' Currently supports "poisson" and "shifted_poisson". Default is "NA".
#' @param pdf_func Pmf of the mixture components associated with the mcmc draws
#' (if mcmc estimation has not been carried out with BayesMultiMode); default is null.
#' @param type Type of modes, either unique or all (the latter includes flat modes); default is "all".
#' @param show_plot If true show the data and estimated modes; default is false.
#' 
#' @return A vector estimated modes.
#' 
#' \insertRef{schaap_genome-wide_2013}{BayesMultiMode}\cr
#' 
#' @importFrom assertthat assert_that
#' @importFrom assertthat is.string
#' 
#' @export
discrete_MF <- function(mcmc, data, pars_names, dist = "NA",
                        pdf_func = NULL, type = "all", show_plot = FALSE){
  
  if (!is.null(pdf_func)) {
    pdf_func <- pdf_func_vec(pdf_func)
  }
  
  x = min(data):max(data)
  
  ## input checks
  assert_that(is.vector(mcmc),
              msg = "mcmc should be a vector")
  assert_that(is.string(dist),
              msg = "dist should be a string")
  assert_that(is.vector(x) & length(x) > 0,
              msg = "x should be a vector of length > 0")
  ##
  
  ##
  names_mcmc = str_to_lower(names(mcmc))
  names_mcmc = str_extract(names_mcmc, "[a-z]+")
  names_mcmc = unique(names_mcmc)
  
  assert_that(sum(pars_names %in% names_mcmc)==length(pars_names),
              msg = "the name of the parameters provided by pars_names and those of the mcmc vector do not match")
  
  if (dist %in% c("shifted_poisson")) {
    assert_that(length(pars_names) == 3,
                msg = "the number of elements in pars_names does not match with dist")
  }
  if (dist %in% c("poisson")) {
    assert_that(length(pars_names) == 2,
                msg = "the number of elements in pars_names does not match with dist")
  }
  ##
  
  pars = c()
  for (i in 1:length(pars_names)) {
    pars = cbind(pars, mcmc[grep(pars_names[i], names(mcmc))])
  }
  
  colnames(pars) <- pars_names
  ##
  
  Khat = nrow(pars)
  
  ### Getting individual component densities
  pdf_k = matrix(0, nrow=length(x), ncol=Khat) 
  for(k in 1:nrow(pars)){
    pdf_k[,k] = pars[k,1] * dist_pdf(x, dist, pars[k, -1, drop = F], pdf_func)
  }
  
  ### summing up to get the mixture
  py <- rowSums(pdf_k, na.rm = T)

  # change in the pdf
  d_py = py[-1] - py[-length(py)]
  
  # where does the pdf decrease ?
  x_decrease = x[d_py<0]
  
  # Only keep the points where the pdf starts to decrease; these are modes
  d2_py = c(0, x_decrease[-1] - x_decrease[-length(x_decrease)])
  x_decrease = x_decrease[which(d2_py!=1)]
  
  # get pdf at these modes 
  pdf_modes = py[x %in% x_decrease]
  
  # get all the points at these peaks (it's unlikely but there might be flat modes) 
  loc_modes = x[which(py %in% pdf_modes)]
  
  if (length(loc_modes) != length(x_decrease)) {
    cat("Warning : Some of these modes are flat.")
  }
  
  if (show_plot) {
    plot(x, py, type = "l", xlab = "", ylab = "")
    for (mode_i in loc_modes) {
      abline(v = mode_i)
    }
  }
  
  output = rep(NA, length(x))
  
  if (type == "unique") {
    output[1:length(loc_modes)] = x_decrease
  }
  
  if (type == "all") {
    output[1:length(loc_modes)] = loc_modes
  }

  return(output)
}