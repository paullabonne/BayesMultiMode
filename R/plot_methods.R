#' Plot method for \code{BayesMixture} objects
#' 
#' Plot an estimated mixture for a given number of draws with a frequency distribution of the data.
#' 
#' @param x An object of class \code{BayesMixture}.
#' @param max_size The number of MCMC draws to plot.
#' @param transparency transparency of the density lines. Default is 0.1. Should be greater than 0 and below or equal to 1.
#' @param ... Not used.
#' 
#' @importFrom posterior as_draws_matrix
#' @importFrom posterior as_draws_df
#' @importFrom scales hue_pal
#' @importFrom ggpubr ggarrange
#' @importFrom assertthat assert_that
#' @importFrom dplyr tibble
#' @importFrom dplyr mutate
#' @importFrom dplyr as_tibble
#' @importFrom dplyr filter
#' @importFrom dplyr left_join
#' @importFrom tidyr gather
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @import ggplot2
#' 
#' @export

plot.BayesMixture <- function(x, max_size = 250, 
                              transparency = 0.1, ...) {
  density <- component <- value <- NULL
  
  assert_that(inherits(x, "BayesMixture"), msg = "input should be an object of class BayesMixture")
  assert_that(is.scalar(transparency) & transparency >= 0 & transparency <= 1,
              msg = "transparency should be a scalar between zero and one")
  
  mcmc = x$mcmc
  dist = x$dist
  y = x$data
  pars_names = x$pars_names

  if (x$dist_type == "continuous") {
    ## plot the data
    g = ggplot(data.frame(y = y), aes(y)) +
      theme_gg +
      theme(legend.position="none") +
      xlab("") + ylab("Density") +
      geom_histogram(aes(y = after_stat(density)),
                     fill="grey33",
                     colour = "white")
    
    ## plot the mixture for each draw
    for (i in sample(nrow(mcmc),min(nrow(mcmc), max_size))) {
      mcmc_i = mcmc[i, ]
      pars = c()
      for (j in 1:length(pars_names)) {
        pars = cbind(pars, mcmc_i[grep(pars_names[j], names(mcmc_i))])
      }
      
      colnames(pars) <- pars_names
      pars = na.omit(pars)
     
      g = g +
        geom_function(fun = dist_mixture,
                      args = list(dist = dist,
                                  pars = pars),
                      alpha = transparency,
                      colour = "#FF6347")
    } 
  }
  
  if (x$dist_type == "discrete") {
    ####### Discrete distribution
    d_y = rep(NA,length(unique(y)))
    for (i in 1:length(d_y)){
      d_y[i] = length(y[y==unique(y)[i]])/length(y)
    }
    
    df_y_temp = tibble(density = d_y,
                       x = unique(y))
    
    x_all = seq(min(y),max(y),1)
    
    mixture_uncertainty = matrix(NA, length(x_all), nrow(mcmc))
    
    if(x$dist %in% c("poisson", "shifted_poisson")){
      eta = mcmc[, grep("eta", colnames(mcmc))]
      lambda = mcmc[, grep("lambda", colnames(mcmc))]
      
      if (x$dist == "shifted_poisson") {
        kappa = mcmc[, grep("kappa", colnames(mcmc))]
      } else {
        kappa = matrix(0,nrow(lambda), ncol(lambda))
      }
      
      for (draw in sample(1:(min(nrow(mcmc), max_size)))) {
        ##
        
        pdf = matrix(0, nrow=length(x_all),ncol=ncol(eta))
        for(j in 1:ncol(eta)){
          if(!is.na(eta[draw,j])){
            pdf[,j] = dpois((x_all-kappa[draw, j, drop = T]), lambda[draw, j, drop = T]) * eta[draw, j, drop = T]
          }
        } 
        
        # summing up to get the mixture
        mixture_uncertainty[,draw] <- rowSums(pdf)
      }
    }
    
    # 
    df_y = tibble(x = seq(min(y),max(y),1)) %>%
      left_join(df_y_temp, by=c("x"="x")) %>%
      mutate(density = ifelse(is.na(density),NA,density)) %>%
      cbind(mixture_uncertainty)
    
    df_y %<>%
      gather(-x,-density,key="component",value="value")
    
    g = ggplot(df_y, aes(x=x)) + 
      theme_gg +
      theme(legend.position="none") +
      xlab("") + ylab("Probability") +
      geom_col(data = filter(df_y,component=="1"),aes(y=density,fill="grey33"),colour="white",alpha=1) +
      geom_line(aes(y=value,colour=component),alpha= transparency) +
      scale_colour_manual(values=rep("#FF6347",length(unique(df_y$component)))) +
      scale_fill_manual(name = "",
                        values = c("grey33"), # Color specification
                        labels = c("Data density"))
    
  }
  
  g
}


#' Plot method for \code{BayesMode} objects
#' 
#' @param x An object of class \code{BayesMode}.
#' @param graphs which plot to show ? Default is all three c("p1", "number", "loc").
#' @param ... Not used.
#' 
#' @importFrom posterior as_draws_matrix
#' @importFrom ggpubr ggarrange
#' @import ggplot2
#' 
#' @export
plot.BayesMode <- function(x, graphs = c("p1", "number", "loc"), ...) {
  Pb <- value <- location_at_modes <- probs_modes <- unique_modes <- prob_nb_modes <- NULL
  
  stopifnot(inherits(x, "BayesMode"))
  assert_that(is.vector(graphs) & is.character(graphs),
              msg = "graphs should be a character vector")
  assert_that(sum(graphs %in% c("p1", "number", "loc"))>=1,
              msg = "graphs should include at least p1, number or loc")
  
  modes = x$modes
  p1 = x$p1
  tb_nb_modes = x$tb_nb_modes
  table_location = x$table_location
  
  df_g0 = tibble(Pb = "Pb",
                 value = (1-p1))
  
  g0 = ggplot(data=df_g0, aes(x=Pb, y=value)) +
    ggtitle("Nb. modes > 1") +
    theme_gg +
    ylim(0, 1) +
    xlab("") + ylab("Posterior probability") +
    geom_bar(stat="identity")
  
  df_g1 = as_tibble(t(table_location))
  
  g1 = ggplot(data=df_g1, aes(x=location_at_modes, y=probs_modes)) +
    theme_gg +
    # scale_x_continuous(breaks=df_g1$possible_nb_modes) +
    ggtitle("Mode locations") +
    xlab("") + ylab("Posterior probability") +
    geom_bar(stat="identity")
  
  if (x$dist_type == "continuous") {
    g1 = g1 + ylim(0, max(df_g1$probs_modes))
  } else {
    g1 = g1 + ylim(0, 1)
  }
  
  df_g2 = as_tibble(t(tb_nb_modes))
  
  g2= ggplot(data=df_g2, aes(x=unique_modes, y=prob_nb_modes)) +
    theme_gg +
    scale_x_continuous(breaks=df_g2$unique_modes) +
    ggtitle("Number of modes") +
    ylim(0, 1) +
    xlab("") + ylab("Posterior probability") +
    geom_bar(stat="identity")
  
  # selecting which graphs to show
  plot_list = list()
  i = 0
  
  widths_p = rep(NA, length(graphs))
  
  if("p1" %in% graphs) {
    i = i + 1
    plot_list[[i]] <- g0
    widths_p[i] = 0.7
  }
  
  if("number" %in% graphs) {
    i = i + 1
    plot_list[[i]] <- g2
    widths_p[i] = 1
  }
  
  if("loc" %in% graphs) {
    i = i + 1
    plot_list[[i]] <- g1
    widths_p[i] = 1
  }
  
  if (i > 1) {
    g <- ggarrange(plotlist = plot_list,
                   ncol = length(graphs), nrow = 1, widths = widths_p)  
  } else {
    g <- plot_list[[i]] 
  }
  
  g
}