#' Plots SFM MCMC output.
#' 
#' Show plots of the MCMC estimation.
#' @param sfm_mcmc a list. Output of `sfm_mcmc_spmix()` containing the parameter draws from the posterior distribution at each MCMC iteration.
#' @param S (a number between 0 and 1) The first S*M draws will be discarded as a burn-in. M is the total number of MCMC iterations.
#' @returns 
#' A list showing :
#' \itemize{
#'        \item 1: Trace plots of pre-processed draws;
#'        \item 2: The posterior distribution of each parameter after post-processing;
#'        \item 3: The posterior probability of the number of components after post-processing;
#'        \item 4: Estimated mixture density after post-processing at each iteration and histogram of the data.
#' }
#' 
#' @import dplyr
#' @import tidyr
#' @importFrom magrittr %>% %<>%
#' @import ggplot2
#' @importFrom ggh4x facet_nested_wrap
#' @importFrom stringr str_remove
#' 
#' @examples
#' # Example with simulated data ================================================
#' #set seed for random number generation
#' set.seed(1) 
#' 
#' # Set the parameters for drawing from a two-component shifted Poisson mixture
#' p1 = 0.3
#' p2 = 1-p1
#' kap1 = 3
#' kap2 = 0
#' lam1 = 1
#' lam2 = 0.5
#' length_data = 70
#' 
#' # Generate data
#' y <- c(rpois(length_data*p1, lam1)+kap1, rpois(length_data*p2, lam2)+kap2)
#' 
#' # Set parameters for the SFM MCMC estimation
#' M = 1000 # Number of MCMC iterations 
#' Jmax = 4 # Maximum number of mixture components
#' 
#' # Estimation with SFM MCMC
#' sfm_mcmc = sfm_mcmc_spmix(y=y,Jmax=Jmax, M=M)
#' 
#' # Generating plots
#' 
#' # Proportion of draws burned in
#' S = 0.5
#' 
#' plots_mcmc(sfm_mcmc, S=S)
#' 
#' # Example with DNA data =====================================================
#' \donttest{
#' y = d4z4
#' M = 5000 # Number of MCMC iterations 
#' Jmax = 10 # Maximum number of mixture components
#' 
#' # Estimation with SFM MCMC
#' sfm_mcmc = sfm_mcmc_spmix(y=y,Jmax=Jmax, M=M)
#' 
#' # Generating plots
#' 
#' # Proportion of draws burned in
#' S = 0.5
#' 
#' plots_mcmc(sfm_mcmc, S=S)
#' }
#' @export

plots_mcmc <- function(sfm_mcmc, S){
  component <- value <- nb <- x <- NULL
  
  if(S<=0 || S>1){
    stop("S should be number greater than 0 and inferior or equal to 1")
  }
  
  # list to store graphs
  graphs = list()
  
  # number of iterations
  M = sfm_mcmc$M
  
  # trace plots
  if(sfm_mcmc$mixt=="shifted_poisson"){
    df_lam = as_tibble(sfm_mcmc$p_draws) %>%
      mutate(draw = rep(1:M)) %>%
      gather(-draw, key=component, value = value) %>%
      mutate(variable = "Lambda")
    
    df_kap = as_tibble(sfm_mcmc$kappa_draws) %>%
      mutate(draw = rep(1:M)) %>%
      gather(-draw, key=component, value = value) %>%
      mutate(variable = "Kappa")
    
    df_p = as_tibble(sfm_mcmc$lambda_draws) %>%
      mutate(draw = rep(1:M)) %>%
      gather(-draw, key=component, value = value) %>%
      mutate(variable = "Mixture weights")
    
    df_parameters = rbind(df_lam,
                          df_kap,
                          df_p) %>%
      mutate(component = str_remove(component,"V"))
  }
  
  g_pre = ggplot(df_parameters) +
    ggtitle("Trace of preprocessed draws") +
    xlim(0,M) +
    geom_line(aes(draw, value, colour=component)) +
    facet_wrap(~variable, scales = "free_y",nrow = 3) +
    xlab(NULL) + ylab(NULL) +
    theme_gg +
    theme(legend.position="none")
  
  g_parameters = ggplot(filter(df_parameters,draw>M*S)) +
    ggtitle("Posterior distribution of each parameter") +
    geom_histogram(aes(x=value, fill=component), colour="white") +
    facet_nested_wrap(~variable + component, scales="free",ncol=length(unique(df_parameters$component))) +
    xlab(NULL) + ylab(NULL) +
    theme_gg +
    theme(legend.position="none")
  
  graphs[[1]] <- g_pre
  
  graphs[[2]] <- g_parameters
  
  # Histogram of K0
  # Here we count the number of components per draw which are non-zero.
  tsnj = sfm_mcmc$snj[(M*S+1):M,]
  K0 = rep(NA, nrow(tsnj))
  for (i in 1:nrow(tsnj)){
    K0[i] = sum(tsnj[i,]!=0)
  }
  
  possible_nb_components = unique(K0)
  post_prob_nb_components = rep(NA,length(possible_nb_components))
  for (i in 1:length(possible_nb_components)){
    post_prob_nb_components[i] = length(K0[K0==possible_nb_components[i]])/length(K0)
  }

  df_g0 = tibble(nb = possible_nb_components,
                 value = post_prob_nb_components)
  
  g2 = ggplot(data=df_g0, aes(x=nb, y=value)) +
    ggtitle("Posterior prob. nb. of components") +
    ylim(0, 1) +
    scale_x_continuous(breaks=possible_nb_components) +
    xlab("") + ylab("") +
    theme_gg +
    geom_bar(stat="identity")
  
  graphs[[3]] <- g2
  
  # generating data frame with density of y
  y = sfm_mcmc$y
  d_y = rep(NA,length(unique(y)))
  for (i in 1:length(d_y)){
    d_y[i] = length(y[y==unique(y)[i]])/length(y)
  }
  unique(y)
  
  df_y_temp = tibble(density = d_y,
                     x = unique(y))
  
  df_y = tibble(x = seq(min(y),max(y),1)) %>%
    left_join(df_y_temp, by=c("x"="x")) %>%
    mutate(density = ifelse(is.na(density),NA,density))
  
  # Uncertainty around the mixture
  
  # burn-in and discarding empty components
  post = post_sfm_mcmc(sfm_mcmc, S)
  theta_draw_slim = post$theta_draws_slim
  J_ne = post$J_ne
  
  mixture_uncertainty = matrix(NA,length(df_y$x),nrow(theta_draw_slim))
  
  if(sfm_mcmc$mixt=="shifted_poisson"){
    p_slim = theta_draw_slim[,1:J_ne]
    kappa_slim = theta_draw_slim[, (J_ne+1):(2*J_ne)]
    lambda_slim = theta_draw_slim[, (2*J_ne+1):(3*J_ne)]
    
    for (draw in 1:nrow(p_slim)){
      pdf.J = matrix(0, nrow=length(df_y$x),ncol=ncol(p_slim))
      for(j in 1:ncol(p_slim)){
        if(!is.na(p_slim[draw,j])){
          pdf.J[,j] = dpois((df_y$x-kappa_slim[draw,j]),lambda_slim[draw,j]) * p_slim[draw,j]
        }
      }
      # summing up to get the mixture
      mixture_uncertainty[,draw] <- rowSums(pdf.J)
    }
  }

  # 
  df_y = tibble(x = seq(min(y),max(y),1)) %>%
    left_join(df_y_temp, by=c("x"="x")) %>%
    mutate(density = ifelse(is.na(density),NA,density)) %>%
    cbind(mixture_uncertainty)
  
  df_y %<>%
    gather(-x,-density,key="component",value="value")
  
  if (M*S>1000){
    alpha=0.01
  } else {
    alpha=0.1
  }
  
  graphs[[4]] = ggplot(df_y, aes(x=x)) + 
    theme_gg +
    ggtitle("Estimated mixture density at each iteration and histogram of the data") +
    theme(legend.position="none") +
    xlab("Repeat units") + ylab("Density") +
    geom_col(data = filter(df_y,component=="1"),aes(y=density,fill="grey33"),colour="white",alpha=1) +
    geom_line(aes(y=value,colour=component),alpha=alpha) +
    scale_colour_manual(values=rep("#FF6347",length(unique(df_y$component)))) +
    scale_fill_manual(name = "",
                    values = c("grey33"), # Color specification
                    labels = c("Data density"))
  
  return(graphs)
}
