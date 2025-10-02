# load package
devtools::load_all()

# set seed for
set.seed(123)

N = 1000
r = 2
p = 0.8
Y = rnbinom(1e6, size = r, prob = 1-p)
mu = r * (1 - (1-p)) / (1-p)
mu - mean(Y)

p = mu/(mu + r)

r = (1-rtilde)/rtilde
rtilde = apply(bayesmix$mcmc, 2, mean)[c("rtilde1", "rtilde2")]
r = (1-rtilde)/rtilde
round(apply(bayesmix$mcmc, 2, mean),2)
set.seed(123)

N = 1000
J = 2
p = c(0.5, 0.5)

true_Mu = c(5, 1)
true_R = c(2, 2)
true_r_tilde = 1 / (1 + true_R)
true_P = true_Mu / (true_Mu + true_R)
P_R = 1 - true_P

Y1 = rnbinom(round(p[1] * N), size = true_R[1], prob = P_R[1])
Y2 = rnbinom(round(p[2] * N), size = true_R[2], prob = P_R[2])

sum(dnbinom(Y1, size = true_R[1], prob = P_R[1], log=TRUE))
Y = c(10 + Y1, 0+Y2)
hist(Y, breaks = 50, main = "Histogram of Simulated Data", xlab = "Counts")

## ============================================================

# estimation
burnin = 4000
bayesmix = bayes_fit(data = Y,
                     K = 2,
                     dist = "shifted_poisson",
                     nb_iter = 5000,
                     burnin = burnin,
                     print = T)

plot(bayesmix, draws = 200)
bayesplot::mcmc_trace(bayesmix$mcmc[, c("kappa1", "kappa2")])
bayesmode = bayes_mode(bayesmix)
bayesmode$p_mode_loc
plot(bayesmode)

average_ll = mean(bayesmix$loglik[burnin:5000])
sd_ll = sd(bayesmix$loglik[burnin:5000])

### Illustrates the benefit of the shift with simulation
### shifted NB vs shifted Poisson
### More tricky: Simulation to show that our NB works better
### measure of performance of the MCMC
### adaptive MH step
