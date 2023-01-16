#' Plot mixture
#' 
#' @param x ...
#' @param max_size ...
#' @param tol_p ...
#' @param colour ...
#' @param transparency ...
#' @param ... ...
#' 
#' @importFrom posterior as_draws_matrix
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

plot.BayesMixture <- function(x, max_size = 200, tol_p = 1e-3,
                              colour = "magenta", transparency = 0.1, ...) {
  component <- value <- NULL
  
  assert_that(inherits(x, "BayesMixture"), msg = "input should be an object of class BayesMixture")
  
  mcmc = x$mcmc
  dist = x$dist
  y = x$data
  pars_names = x$pars_names
  
  if (x$dist_type == "continuous") {
    ## plot the data
    g = ggplot(data.frame(y = y), aes(y)) +
      geom_histogram(aes_string(y = "..density.."),
                     colour = "white",
                     bins = 70)
    
    ## plot the mixture for each draw
    for (i in sample(1:(min(nrow(mcmc), max_size)))) {
      mcmc_i = mcmc[i, ]
      pars = c()
      for (i in 1:length(pars_names)) {
        pars = cbind(pars, mcmc_i[grep(pars_names[i], colnames(mcmc_i))])
      }

      colnames(pars) <- pars_names
      
      est_mode = rep(NA, nrow(pars))
      
      pars = pars[pars[,1] > tol_p, ]

      g = g +
        geom_function(fun = dist_mixture,
                      args = list(dist = dist,
                                  pars = pars),
                      alpha = transparency,
                      colour = colour)
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
    
    if(x$dist=="shifted_poisson"){
      theta = mcmc[, grep("theta", colnames(mcmc))]
      kappa = mcmc[, grep("kappa", colnames(mcmc))]
      lambda = mcmc[, grep("lambda", colnames(mcmc))]
      
      for (draw in sample(1:(min(nrow(mcmc), max_size)))) {
        pdf.J = matrix(0, nrow=length(x_all),ncol=ncol(theta))
        for(j in 1:ncol(theta)){
          if(!is.na(theta[draw,j])){
            pdf.J[,j] = dpois((x_all-kappa[draw, j, drop = T]),lambda[draw, j, drop = T]) * theta[draw, j, drop = T]
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
    
    g = ggplot(df_y, aes(x=x)) + 
      theme_gg +
      ggtitle("Estimated mixture density at each iteration and histogram of the data") +
      theme(legend.position="none") +
      xlab("Repeat units") + ylab("Density") +
      geom_col(data = filter(df_y,component=="1"),aes(y=density,fill="grey33"),colour="white",alpha=1) +
      geom_line(aes(y=value,colour=component),alpha= transparency) +
      scale_colour_manual(values=rep("#FF6347",length(unique(df_y$component)))) +
      scale_fill_manual(name = "",
                        values = c("grey33"), # Color specification
                        labels = c("Data density"))
    
  }
  
  g
}


#' Plot modes
#' @param x ...
#' @param colour ...
#' @param ... ...
#' 
#' @importFrom posterior as_draws_matrix
#' @importFrom ggpubr ggarrange
#' @import ggplot2
#' 
#' @export
plot.BayesMode <- function(x, colour = "magenta", ...) {
  Pb <- value <- location_at_modes <- probs_modes <- unique_modes <- prob_nb_modes <- NULL
  
  stopifnot(inherits(x, "BayesMode"))
  
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
  
  g <- ggarrange(g0, g2, g1,
                 ncol = 3, nrow = 1, widths = c(0.7,1, 1))
  
  g
}