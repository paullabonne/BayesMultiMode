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
  
  return(graphs)
}
