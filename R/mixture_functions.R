#' @importFrom sn dst
#' @importFrom sn dsn
#' @importFrom stats dnorm
#' 

# Mixture of pdf_func
#' @keywords internal
pdf_func_mix <- function(x, pars, pdf_func) {
  pdf_func = match.fun(pdf_func) #solves NOTE "pdf_func is undefined"
  
  mixture = 0
  for (i in 1:nrow(pars)) {
    mixture = mixture + pars[i, "eta"] * pdf_func(x, pars[i,])
  }
  return(mixture)
}

#' @keywords internal
test_and_export <- function(p, pdf_func, dist, pars_names, dist_type, par_type, loc) {
  
  list_func = list()
  
  # that pdf_func can be computed when provided
  if (!is.null(pdf_func)) {
    assert_that(!is.na(pdf_func(1, vec_to_mat(p, pars_names)[1,])),
                msg = "running pdf_func with pars provided returns NA") 
    assert_that(!is.na(dist_type),
                msg = "dist_type must be provided when argument pdf_func is used")
    assert_that(dist_type %in% c("continuous", "discrete"),
                msg = "dist_type should be either discrete or continuous")
  }
  
  msg_0 = paste0("variable names in ", par_type, " should be ")
  
  if (!is.na(dist)) {
    
    if (dist == "poisson"){
      assert_that(sum(pars_names %in% c("eta", "lambda"))==2,
                  msg = paste0(msg_0, "eta and lambda when dist = poisson"))
      pdf_func <- function(x, pars) dpois(x, pars["lambda"])
    }
    
    if (dist == "shifted_poisson"){
      assert_that(sum(pars_names %in% c("eta", "kappa", "lambda"))==3,
                  msg = paste0(msg_0, "eta and lambda when dist = shifted_poisson"))
      pdf_func <- function(x, pars) dpois(x - pars["kappa"], pars["lambda"])
    }
    
    if (dist == "normal"){
      assert_that(sum(pars_names %in% c("eta", "mu", "sigma"))==3,
                  msg = paste0(msg_0, "eta, mu and sigma when dist = normal"))
      pdf_func <- function(x, pars) dnorm(x, pars["mu"], pars["sigma"])
    }
    
    if (dist == "skew_normal"){
      assert_that(sum(pars_names %in% c("eta", "xi", "omega", "alpha"))==4,
                  msg = paste0(msg_0, "eta, xi, omega and alpha when dist = skew_normal"))
      pdf_func <- function(x, pars) dsn(x, pars["xi"], pars["omega"], pars["alpha"])
      loc = "xi"
    }
    
    if (dist %in% c("normal", "skew_normal")) {
      dist_type = "continuous"
    } else if (dist %in% c("poisson", "shifted_poisson")) {
      dist_type = "discrete"
    } else {
      stop("Unsupported distribution; dist should be either normal, skew_normal, poisson or shifted_poisson")
    } 
  }
  
  if ((is.na(dist) | !(dist %in% c("normal", "skew_normal"))) & dist_type != "discrete") {
    assert_that(!is.na(loc) & is.string(loc),
                msg = "loc must be given when using a continuous distribution other than the normal distribution")
    assert_that(loc %in% pars_names[pars_names!="eta"],
                msg = paste0("loc must a parameter included in ", par_type, " other than eta"))
  }
  
  list_func$dist_type = dist_type
  list_func$pdf_func = pdf_func
  list_func$loc = loc
  
  return(list_func)
}

#' @keywords internal
vec_to_mat <- function(pars, pars_names) {
  pars_mat = c()
  for (i in 1:length(pars_names)) {
    pars_mat = cbind(pars_mat, pars[grep(pars_names[i], names(pars))])
  }
  colnames(pars_mat) <- pars_names
  
  return(pars_mat)
}