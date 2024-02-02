#' Plot method for `bayes_mixture` objects
#' 
#' Plot an estimated mixture for a given number of draws with a frequency distribution of the data.
#' 
#' @param x An object of class `bayes_mixture`.
#' @param draws The number of MCMC draws to plot.
#' @param bins (for continuous mixtures) Number of bins for the histogram of
#' the data. Passed to `geom_histogram()`.
#' @param alpha transparency of the density lines. Default is 0.1. Should be greater than 0 and below or equal to 1.
#' @param ... Not used.
#' 
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
plot.bayes_mixture <- function(x, draws = 250,
                               bins = 30,
                               alpha = 0.1, ...) {
  density <- component <- value <- NULL
  
  assert_that(is.scalar(alpha) & alpha >= 0 & alpha <= 1,
              msg = "alpha should be a scalar between zero and one")
  
  mcmc = x$mcmc
  pdf_func = x$pdf_func
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
                     colour = "white",
                     bins = bins)
    
    ## plot the mixture for each draw
    for (i in sample(nrow(mcmc),min(nrow(mcmc), draws))) {
      pars = vec_to_mat(mcmc[i, , drop = T], pars_names)
      pars = na.omit(pars)
      
      g = g +
        geom_function(fun = pdf_func_mix,
                      args = list(pdf_func = pdf_func,
                                  pars = pars),
                      alpha = alpha,
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
    
    mixture_uncertainty = matrix(NA, length(x_all), draws)
    j = 1
    for (i in sample(nrow(mcmc),min(nrow(mcmc), draws))) {
      ##
      pars = vec_to_mat(mcmc[i, ], pars_names)
      pars = na.omit(pars)
      
      mixture_uncertainty[,j] = pdf_func_mix(x_all, pars, pdf_func)
      j = j+1
    }
    
    # 
    df_y = tibble(x = seq(min(y),max(y),1)) %>%
      left_join(df_y_temp, by=c("x"="x")) %>%
      mutate(density = ifelse(is.na(density),0,density)) %>%
      cbind(mixture_uncertainty)
    
    df_y %<>%
      gather(-x,-density,key="component",value="value")
    
    g = ggplot(df_y, aes(x=x)) + 
      theme_gg +
      theme(legend.position="none") +
      xlab("") + ylab("Probability") +
      geom_col(data = filter(df_y,component=="1"),aes(y=density,fill="grey33"),colour="white",alpha=1) +
      geom_line(aes(y=value,colour=component),alpha= alpha) +
      scale_colour_manual(values=rep("#FF6347",length(unique(df_y$component)))) +
      scale_fill_manual(name = "",
                        values = c("grey33"), # Color specification
                        labels = c("Data density"))
    
  }
  
  g
}


#' Plot method for `bayes_mode` objects
#' 
#' @param x An object of class `bayes_mode`.
#' @param graphs which plot to show ? Default is all three c("p1", "number", "loc").
#' @param draw Plot modes in a given mcmc draw; note that `graphs` is discarded. Default is `NULL`.
#' @param ... Not used.
#' 
#' @importFrom ggpubr ggarrange
#' @importFrom assertthat is.scalar
#' @import ggplot2
#' 
#' @export
plot.bayes_mode <- function(x, graphs = c("p1", "number", "loc"), draw = NULL, ...) {
  Pb <- value <- `posterior probability` <- `number of modes` <- `mode location` <- NULL
  
  if (!is.null(draw)) {
    assert_that(is.scalar(draw),
                msg = paste0("draw should be a scale greater than zero and ",
                             "inferior to the number of mcmc draws used in",
                             deparse(substitute(x))))
    
    BayesMix = x$BayesMix
  
    mix = mixture(BayesMix$mcmc[draw,],
            dist = BayesMix$dist,
            pdf_func = BayesMix$pdf_func,
            dist_type = BayesMix$dist_type,
            loc = BayesMix$loc,
            range = x$range)
    
    modes = mix_mode(mix)
    
    plot(modes)
  } else {
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
    
    g1 = ggplot(data=df_g1, aes(x = `mode location`, y = `posterior probability`)) +
      theme_gg +
      # scale_x_continuous(breaks=df_g1$possible_nb_modes) +
      ggtitle("Mode locations") +
      xlab("") + ylab("Posterior probability") +
      geom_bar(stat="identity")
    
    if (x$dist_type == "continuous") {
      g1 = g1 + ylim(0, max(df_g1$`posterior probability`))
    } else {
      g1 = g1 + ylim(0, 1)
    }
    
    df_g2 = as_tibble(t(tb_nb_modes))
    
    g2= ggplot(data=df_g2, aes(x = `number of modes`, y = `posterior probability`)) +
      theme_gg +
      scale_x_continuous(breaks=df_g2$`number of modes`) +
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
  
}

#' Plot method for `mixture` objects
#' 
#' @param x An object of class `mixture`.
#' @param from the lower limit of the range over which the function will be plotted.
#' @param to the upper limit of the range over which the function will be plotted.
#' @param ... Not used.
#' 
#' @importFrom graphics curve
#' 
#' @export
plot.mixture <- function(x, from = NULL, to = NULL, ...) {
  assert_that(is.null(from)|(is.numeric(from) && is.finite(to)),
              is.null(to)|(is.numeric(to) && is.finite(to)),
              msg = "arguments from and to must be numeric and finite")
  
  pars = x$pars
  mode_est = x$mode_estimates
  pdf_func = x$pdf_func
  dist = x$dist
  
  if (dist %in% c("normal", "skew_normal")) {
    par_names = str_extract(names(pars), "[a-z]+")
    mu = pars[par_names %in% c("mu", "xi")]
    sigma = pars[par_names %in% c("sigma", "omega")]
    
    # calculate min and max for x axis
    if (is.null(from)) {
      min_mu = min(mu)
      min_sigma = sigma[mu == min_mu]
      min_x = min_mu - 4*min_sigma
    } else {
      min_x = from
    }
    
    if (is.null(to)) {
      max_mu = max(mu)
      max_sigma = sigma[mu == max_mu]
      max_x = max_mu + 4*max_sigma
    } else {
      max_x = to
    }
    
    pars = vec_to_mat(pars, par_names)
    curve(pdf_func_mix(x, pars, pdf_func), from = min_x,
          to = max_x, xlab = "", ylab = "")
    
  } else if (x$dist_type == "continuous") {
    assert_that(!is.null(from) & !is.null(to),
                msg = "arguments from and to must be provided when plotting
                a continuous distribution other than the normal and skew
                normal provided by BayesMultiMode.")
    
    assert_that(is.finite(from) && is.finite(to),
                msg = "from and to must be finite")
    
    par_names = str_extract(names(pars), "[a-z]+")
    pars = vec_to_mat(pars, par_names)
    curve(pdf_func_mix(x, pars, pdf_func), from = from,
          to = to, xlab = "", ylab = "")
    
  } else if (x$dist_type  == "discrete") {
    assert_that(!is.null(from) & !is.null(to),
                msg = "arguments from and to must be provided when plotting
                a discrete mixture.")
    assert_that(is.finite(from) && is.finite(to),
                msg = "from and to must be finite")
    xx = round(from):round(to)
    par_names = str_extract(names(pars), "[a-z]+")
    pars_mat = vec_to_mat(pars, par_names)
    py = pdf_func_mix(xx, pars_mat, pdf_func)
    
    plot(xx, py, type = "h", xlab = "", ylab = "", lwd = 4,
         xlim = c(from, to))
  }
}

#' Plot method for `mix_mode` objects
#' 
#' @param x An object of class `mix_mode`.
#' @param from the lower limit of the range over which the function will be plotted.
#' @param to the upper limit of the range over which the function will be plotted.
#' @param ... Not used.
#' 
#' @importFrom graphics curve abline
#' 
#' @export
plot.mix_mode <- function(x, from = NULL, to = NULL, ...) {
  mix = mixture(x$pars, dist = x$dist,
                pdf_func = x$pdf_func,
                dist_type = x$dist_type,
                range = x$range)
  
  plot(mix, from = from, to = to)
  for (m in x$mode_estimates) {
    abline(v = m, col = "red")
  }
}

### ggplot theme
#' @keywords internal
theme_gg <-  ggplot2::theme_bw()+ ggplot2::theme(strip.background=element_blank(),
                                                 strip.text=element_text(size=11),
                                                 title=element_text(size=11),
                                                 panel.border=element_blank(),
                                                 # panel.grid.major=element_blank(),
                                                 #panel.grid.minor=element_blank(),
                                                 #legend.key=element_rect(colour="white"),
                                                 legend.position="bottom",
                                                 legend.box.margin=margin(-10,0,-10,0),
                                                 legend.title=element_blank(),
                                                 legend.text=element_text(size=11),
                                                 axis.text=element_text(size=11),
                                                 axis.line.y = element_line(colour = 'grey', size=0.5, linetype="solid"),
                                                 axis.line.x = element_line(colour = 'grey', size=0.5, linetype="solid"),
                                                 plot.title=element_text(hjust=0.5, size=12, face="bold"))

