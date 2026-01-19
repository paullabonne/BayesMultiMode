options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
library("BayesMultiMode")


###################################################
### galaxy data
###################################################

# Loading the galaxy data
y = galaxy

set.seed(123)
library("multimode") # for conventional multimodality tests

# small function to run multimodality tests
tests_mode <- function(y) {
  
  tests = c(modetest(y, method = "SI")$p.value,
            modetest(y, method = "HY")$p.value,
            modetest(y, method = "FM")$p.value,
            modetest(y, method = "HH")$p.value,
            modetest(y, method = "CH")$p.value,
            modetest(y, method = "ACR")$p.value)
  
  names(tests) = c("SI", "HY", "FM", "HH", "CH", "ACR")
  
  return(tests)
}

# run multimodality tests on the galaxy series
tests_mode(y)


## Bayesian analysis
# sfm mcmc on the galaxy data
mix_mcmc = bayes_fit(data = y,
                     K = 10,
                     dist = "normal",
                     nb_iter = 2000,
                     burnin = 1000,
                     print = FALSE)


# plotting the estimated mixture in 100 draws
plot(mix_mcmc, draws = 100, alpha = 0.2)


# bayesian mode inference using mcmc results
mode_mcmc = bayes_mode(mix_mcmc)

summary(mode_mcmc)


# graphical output of the Bayesian mode inference
plot(mode_mcmc, draws = 100, alpha = 0.2)


# illustrating the effect of changing the rounding precision for mode locations (with argument rd)
library("ggpubr") 
library("ggplot2") 

mode_mcmc_a = bayes_mode(mix_mcmc, rd = 1)
mode_mcmc_b = bayes_mode(mix_mcmc, rd = 0)

plot_a = plot(mode_mcmc_a, graphs = "loc") +
  ggtitle("Rounding at the first decimal")

plot_b = plot(mode_mcmc_b, graphs = "loc") +
  ggtitle("Rounding at the unit level")

ggarrange(plot_a, plot_b,
          ncol = 1, nrow = 2, widths = c(0.5, 0.5))


# showing the results from the mode finding algorithm at one mcmc iteration
par(mar = c(2, 2, 2, 2))
mode_mcmc = plot(mode_mcmc, draw = 1)


###################################################
### world economic growth data
###################################################
library("pwt10") # world gdp per capita data
library("dplyr")

# small function for running the Bayesian mode inference repeatedly
estimation_growth <- function(start, end) {
  
  # retrieving economic growth data
  y = pwt10.0 %>%
    select(year, country, rgdpe, pop) %>%
    filter(year %in% rep(start:end, 1)) %>%
    group_by(country) %>%
    summarise(rgdpe = as.numeric(mean(rgdpe / pop, na.rm = T) / 1000)) %>%
    na.omit() %>%
    select(rgdpe) %>%
    unlist()
  
  # omitting the 5 largest economies
  y = sort(y, decreasing = T)
  y = y[-c(1:5)]
  
  # mcmc estimation
  mix = bayes_fit(y, dist = "normal", K = 10, print = F)
  
  # plotting the estimated mixture in 100 draws
  plot_mix = plot(mix, draws = 100, alpha = 0.1) + ylab(NULL) +
    ggtitle(paste0(start, "'s"))
  
  # mode inference
  modes = bayes_mode(mix, rd = 0)

  # conventional tests
  tests = tests_mode(y)
  
  return(list(mix = mix,
              plot_mix = plot_mix,
              modes = modes,
              tests = tests))
}

# mcmc for each decade
res_60s = estimation_growth(1960, 1969)
res_70s = estimation_growth(1970, 1979)
res_80s = estimation_growth(1980, 1989)
res_90s = estimation_growth(1990, 1999)
res_00s = estimation_growth(2000, 2009)
res_10s = estimation_growth(2010, 2019)

# plot
ggarrange(res_60s$plot,
          res_70s$plot,
          res_80s$plot,
          res_90s$plot,
          res_00s$plot,
          res_10s$plot,
          common.legend = T,
          ncol = 1, nrow = 6)


# retrieving posterior probabilities of location points being modes
t60 = as_tibble(t(res_60s$modes$p_mode_loc)) %>%
  mutate(decade = "1960's")

t70 = as_tibble(t(res_70s$modes$p_mode_loc)) %>%
  mutate(decade = "1970's")

t80 = as_tibble(t(res_80s$modes$p_mode_loc)) %>%
  mutate(decade = "1980's")

t90 = as_tibble(t(res_90s$modes$p_mode_loc)) %>%
  mutate(decade = "1990's")

t00 = as_tibble(t(res_00s$modes$p_mode_loc)) %>%
  mutate(decade = "2000's")

t10 = as_tibble(t(res_10s$modes$p_mode_loc)) %>%
  mutate(decade = "2010's")

locations = rbind(t60, t70, t80, t90, t00, t10)

# plotting
ggplot(data=locations, aes(x = `mode location`, y = `posterior probability`)) +
  facet_wrap(~ decade, nrow = 6) +
  theme_minimal() +
  xlab("") + ylab("Posterior probability") +
  geom_bar(stat="identity")


# showing conventional tests and posterior probability of unimodality for each decade
library("xtable")
tests = rbind(c(res_60s$tests, res_60s$modes$p1),
              c(res_70s$tests, res_70s$modes$p1),
              c(res_80s$tests, res_80s$modes$p1),
              c(res_90s$tests, res_90s$modes$p1),
              c(res_00s$tests, res_00s$modes$p1),
              c(res_10s$tests, res_10s$modes$p1))

rownames(tests) = c("1960's", "1970's", "1980's",
                    "1990's", "2000's", "2010's")
colnames(tests)[7] = "BayesMultiMode"

print(xtable(tests, caption = "P-values from conventional tests of multimodality for the Penn World Tables data (columns 1-6). Null hypothesis of unimodality with alternative hypothesis of at least two modes. The 7th colmn shows the posterior probability of multimodality from BayesMultiMode.", label = "tb:tests"), caption.placement = "bottom", include.rownames=TRUE)


# posterior probabilities attached to the number of modes for each decade
library("tidyr")

nb_modes = mutate(as_tibble(t(res_60s$modes$p_nb_modes)),
                  decade = "1960's") %>%
  rbind(mutate(as_tibble(t(res_70s$modes$p_nb_modes)),
               decade = "1970's")) %>%
  rbind(mutate(as_tibble(t(res_80s$modes$p_nb_modes)),
               decade = "1980's")) %>%
  rbind(mutate(as_tibble(t(res_90s$modes$p_nb_modes)),
               decade = "1990's")) %>%
  rbind(mutate(as_tibble(t(res_00s$modes$p_nb_modes)),
               decade = "2000's")) %>%
  rbind(mutate(as_tibble(t(res_10s$modes$p_nb_modes)),
               decade = "2010's")) %>%
  spread(key = "number of modes", value = "posterior probability") %>%
  select(-decade) %>%
  as.matrix()

nb_modes[is.na(nb_modes)] = 0
rownames(nb_modes) = c("1960's", "1970's", "1980's",
                    "1990's", "2000's", "2010's")

print(xtable(nb_modes, caption = "Posterior probabilities attached to the number of modes for the Penn World Tables data.", label = "tb:nb_modes"), caption.placement = "bottom", include.rownames=TRUE)


###################################################
### DNA data
###################################################

# conventional tests
table_tests = rbind(tests_mode(d4z4),
                    tests_mode(ct47))

row.names(table_tests) = c("d4z4", "ct47")

table_tests


# mcmc estimation
mcmc_d4z4 = bayes_fit(data = d4z4,
                      K = 20,
                      dist = "shifted_poisson",
                      nb_iter = 2000,
                      burnin = 1000,
                      print = FALSE)

mcmc_ct47 = bayes_fit(data = ct47,
                      K = 10,
                      dist = "shifted_poisson",
                      nb_iter = 2000,
                      burnin = 1000,
                      print = FALSE)


# plotting estimated mixtures at 100 draws
p1 = plot(mcmc_d4z4, draws = 100, alpha = 0.2) +
  ggtitle("d4z4")

p2 = plot(mcmc_ct47, draws = 100, alpha = 0.2) +
  ggtitle("ct47")

ggarrange(p1, p2,
          ncol = 2, nrow = 1, widths = c(0.5, 0.5))


# plotting Bayesian mode inference results
p1 = plot(bayes_mode(mcmc_d4z4)) +
  ggtitle("d4z4")

p2 = plot(bayes_mode(mcmc_ct47)) +
  ggtitle("ct47")

ggarrange(p1, p2, labels = c("d4z4","ct47"),
          ncol = 1, nrow = 2, widths = c(0.5, 0.5))


###################################################
### cyclone data
###################################################
library("BNPmix") # for external mcmc estimation

# retrieving the cyclone data (maximum intensity the Eastern North Pacific basin since 1981)
y = cyclone %>%
  filter(BASIN == "EP",
         SEASON > "1981") %>%
  select(max_wind) %>%
  unlist()

# Bayesian MCMC mixture estimation using the BNPmix package
PY_result = PYdensity(y,
                      mcmc = list(niter = 2000,
                                  nburn = 1000,
                                  print_message = FALSE),
                      output = list(out_param = TRUE))


# transforming the results to have one variable per column and draws as rows
mcmc_py = list()

for (i in 1:length(PY_result$p)) {
  k = length(PY_result$p[[i]][, 1])
  
  draw = c(PY_result$p[[i]][, 1],
           PY_result$mean[[i]][, 1],
           sqrt(PY_result$sigma2[[i]][, 1]),
           i)
  
  names(draw)[1 : k] = paste0("eta", 1 : k)
  names(draw)[(k + 1):(2*k)] = paste0("mu", 1 : k)
  names(draw)[(2 * k + 1) : (3 * k)] = paste0("omega", 1 : k)
  names(draw)[3 * k + 1] = "draw"
  
  mcmc_py[[i]] = draw
}

mcmc_py = as.matrix(bind_rows(mcmc_py))


# creating an object of class BayesMixture
py_BayesMix = bayes_mixture(mcmc = mcmc_py,
                            data = y,
                            burnin = 0,
                            # the burnin has already been discarded
                            dist = "normal",
                            vars_to_keep = c("eta", "mu", "omega"),
                            vars_to_rename = c("sigma" = "omega"))


# plotting the estimated mixture in 100 draws
plot(py_BayesMix, draws = 100, alpha = 0.2)


# Bayesian mode inference
bayesmode = bayes_mode(py_BayesMix)


# Plotting the results of the Bayesian mode inference
plot(bayesmode)


# posterior probability of multimodality
1 - bayesmode$p1


