#' Bayesian estimation of mixture distributions
#' 
#' Fit a univariate mixture using sparse finite mixture (SFM) or telescoping MCMC.
#' 
#' @param data Vector of observations.
#' @param K Maximum number of count components when `zero_inflated = TRUE`;
#' otherwise, number of SFM slots or upper bound on the sampled telescoping count.
#' @param dist String indicating the distribution of the mixture components;
#' currently supports `"normal"`, `"skew_normal"`, `"poisson"`, `"shifted_poisson"`,
#' `"neg_binomial"` and `"shifted_neg_binomial"`.
#' @param priors List of priors; default is an empty list which implies the following priors:\cr
#' `a0 = 1`,\cr `A0 = 200`,\cr `b0 = median(y)`,\cr `B0 = (max(y) - min(y))^2` (normal),\cr
#' `D_xi = 1`,\cr `D_psi =1`, (skew normal: `B0 = diag(D_xi,D_psi)`), \cr `c0 = 2.5`,\cr
#' `l0 = 1.1` (poisson),\cr `l0 = 5` (shifted poisson),\cr `L0 = 1.1/median(y)`,\cr `L0 = l0 - 1` (shifted poisson),\cr
#' `g0 = 0.5`,\cr `G0 = 100*g0/c0/B0` (normal),\cr 
#' `G0 = g0/(0.5*var(y))` (skew normal).
#' @param nb_iter Number of MCMC iterations; default is `2000`.
#' @param burnin Number of MCMC iterations used as burnin; default is `nb_iter/2`.
#' @param print Showing MCMC progression ? Default is `TRUE`.
#' @param zero_inflated Add a structural-zero component to a Poisson or negative-binomial
#' mixture (shifted or unshifted)? Default is `FALSE`. Its probability has a
#' Beta prior with positive shape parameters `priors$zero_a` and
#' `priors$zero_b`, both defaulting to `1`. Only `sampler = "sfm"` supports
#' zero inflation; the telescoping sampler does not yet support it.
#' @param sampler `"sfm"` (default) fixes the total number of components.
#' `"telescoping"` samples that count for every supported `dist` and requires
#' `zero_inflated = FALSE`.
#' @param K_init Initial total count for telescoping, from 1 to `K`. May exceed the number of distinct observations.
#' @param weight_model Telescoping Dirichlet model: `"dynamic"` uses concentration divided by the current count per component; `"static"` uses concentration per component.
#' @param K_prior Telescoping count prior: `"bnb"` gives `K - 1` a beta-negative-binomial(1, 4, 3) prior; `"uniform"` gives every count equal mass. Both priors are truncated to the integers from 1 to the supplied cap.
#' @param concentration `"default"` follows the prior helpers in the authors' `telescope` package.
#' With dynamic weights it fixes concentration at 1; use `"F_6_3"` to sample it
#' from the F(6, 3) prior suggested for clustering. With static weights,
#' `"default"` and `"G_1_20"` both use a Gamma prior with shape 1 and rate 20.
#' A positive number fixes concentration for either model; `NULL` samples it
#' from the corresponding prior. SFM uses `a0` and `A0` instead.
#' 
#' @return A list of class \code{bayes_mixture} containing:
#'  \item{data}{Same as argument.}
#'  \item{mcmc}{Matrix of MCMC draws where the rows corresponding to burnin have been discarded;}
#'  \item{mcmc_all}{Matrix of MCMC draws.}
#'  \item{loglik}{Log likelihood at each MCMC draw.}
#'  \item{K}{Storage width, including the structural-zero component if enabled. For telescoping this is the cap, not a draw of the count.}
#'  \item{state_all, state}{For telescoping, matrices of total count (`K`), occupied count (`Kplus`) and concentration, before and after burn-in. Neither count is a count of modes.}
#'  \item{cap_mass}{For telescoping, the fraction of retained draws at the cap. A large value warrants checking a larger cap.}
#'  \item{zero_inflated}{Present and `TRUE` when structural zeros are enabled.}
#'  \item{count_K}{When zero inflation is enabled, the requested number of count components.}
#'  \item{dist}{Same as argument.}
#'  \item{pdf_func}{The pdf/pmf of the mixture components.}
#'  \item{dist_type}{Type of the distribution, i.e. continuous or discrete.}
#'  \item{pars_names}{Names of the mixture components' parameters.}
#'  \item{loc}{Name of the location parameter of the mixture components.}
#'  \item{nb_var}{Number of variables/parameters in the mixture distribution.}
#' 
#' @details
#' The telescoping sampler follows Algorithms 2 and 3 of
#' \doi{10.1214/21-BA1294}. Each draw includes weights and parameters for
#' occupied and empty components. Slots beyond its sampled count contain `NA`.
#' For constant Gaussian data the range-based variance prior is undefined;
#' telescoping uses `B0 = 1` unless a positive variance is supplied.
#' For constant or single-observation skew-normal data, telescoping uses a
#' variance of 1 to set the defaults for `C0` and `G0`.
#' For shifted and unshifted negative binomial components the mean prior uses the maximum
#' of the full data, not the cluster maximum. `priors$sig` sets the initial
#' proposal scale (default 0.1); adaptation stops at the end of burn-in.
#'
#' 
#' Let \eqn{y_i}, \eqn{i=1,\dots,n} denote observations.
#' A general mixture of \eqn{K} distributions from the same 
#' parametric family is given by:
#' \deqn{y_i \sim \sum_{k=1}^{K}\pi_k p(\cdot|\theta_k)}
#' with \eqn{\sum_{k=1}^{K}\pi_k=1} and \eqn{\pi_k\geq 0}, \eqn{k=1, ...,K}.
#' \cr\cr
#' The exact number of components does not have to be known *a priori*
#' when using an SFM MCMC approach. Rather, an upper bound is specified for the
#' number of components and the weights of superfluous components are shrunk
#' towards zero during estimation. Following \insertCite{malsiner-walli_model-based_2016;textual}{BayesMultiMode}
#' a symmetric Dirichlet prior is used for the mixture weights:
#' \deqn{\pi_k \sim \text{Dirichlet}(e_0,\dots,e_0),}
#' where a Gamma hyperprior is used on the concentration parameter \eqn{e_0}:\cr\cr
#' \deqn{e_0 \sim \text{Gamma}\left(a_0, A_0\right).}
#' 
#' **Mixture of Normal distributions**
#' 
#' Normal components take the form:
#' \deqn{p(y_i|\mu_k,\sigma_k) = \frac{1}{\sqrt{2 \pi} \
#'   \sigma_k} \exp\left( - \, \frac{1}{2}            \left(  \frac{y_i -
#'       \mu_k}{\sigma_k} \right)^2     \right).}
#' 
#' Independent conjugate priors are used for \eqn{\mu_k} and \eqn{\sigma^2_k}
#' (see for instance Malsiner-Walli et al. 2016):
#' \deqn{\mu_k \sim \text{Normal}( \text{b}_0, \text{B}_0),}
#' \deqn{\sigma^{-2}_k \sim \text{Gamma}( \text{c}_0, \text{C}_0),}
#' \deqn{C_0 \sim \text{Gamma}( \text{g}_0, \text{G}_0).}
#' 
#' 
#' **Mixture of skew-Normal distributions**
#' 
#' We use the skew-Normal of \insertCite{azzalini_1985;textual}{BayesMultiMode} which takes the form:
#' \deqn{p(y_i| \xi_k,\omega_k,\alpha_k) = \frac{1}{\omega_k\sqrt{2\pi}} \ \exp\left( - \,
#' \frac{1}{2}            \left(  \frac{y_i - \xi_k}{\omega_k} \right)^2\right) \
#' \left(1 + \text{erf}\left( \alpha_k\left(\frac{y_i - \xi_k}{\omega_k\sqrt{2}}\right)\right)\right),}
#' where \eqn{\xi_k} is a location parameter, \eqn{\omega_k} a scale parameter and \eqn{\alpha_k}
#' the shape parameter introducing skewness. For Bayesian estimation, we adopt the approach of
#' \insertCite{fruhwirth-schnatter_bayesian_2010;textual}{BayesMultiMode} and use the following reparameterised random-effect model:
#' \deqn{z_i \sim TN_{[0,\infty)}(0, 1),}
#' \deqn{y_i|(S_i = k) = \xi_k + \psi_k z_i + \epsilon_i, \quad \epsilon_i \sim N(0, \sigma^2_k),}
#' where the parameters of the skew-Normal are recovered with
#' \deqn{\alpha_k = \frac{\psi_k}{\sigma_k}, \qquad \omega^2_k = \sigma^2_k + \psi^2_k.}
#' By defining a regressor \eqn{x_i = (1, z_i)'}, the skew-Normal mixture can be seen as
#' random effect model and sampled using standard techniques. Thus we use priors similar to
#' the Normal mixture model:
#' \deqn{(\xi_k, \psi_k)' \sim \text{Normal}(\text{b}_0, \text{B}_0),}
#' \deqn{\sigma^{-2}_k \sim \text{Gamma}(\text{c}_0, \text{C}_0),}
#' \deqn{\text{C}_0 \sim \text{Gamma}( \text{g}_0, \text{G}_0).}
#' We set \deqn{\text{b}_0 = (\text{median}(y), 0)'} and \deqn{\text{B}_0 = \text{diag}(\text{D}\_\text{xi}, \text{D}\_\text{psi})} with D_xi = D_psi = 1.
#' 
#' 
#' **Mixture of Poisson distributions**
#' 
#' Poisson components take the form:
#' \deqn{p(y_i|\lambda_k) = \frac{1}{y_i!} \, \lambda^{y_i}_k \,\exp(-\lambda_k).}
#' The prior for \eqn{\lambda_k} follows from \insertCite{viallefont2002bayesian;textual}{BayesMultiMode}:
#' \deqn{\lambda_k \sim \text{Gamma}(\text{l}_0,\text{L}_0).}
#' For `dist = "poisson"` with zero inflation, the default rate prior uses
#' `L0 = 1.1/median(data[data > 0])`, or `1.1` for all-zero data.
#' Telescoping uses the same fallback when the sample median is zero.
#'
#' **Zero-inflated count mixtures**
#'
#' With `zero_inflated = TRUE`, the model is
#' \deqn{p(y) = \pi_0 1\{y=0\} + (1-\pi_0)\sum_{k=1}^{K}\eta_k f_k(y).}
#' The count weights sum to one and retain the sparse Dirichlet prior.
#' The structural-zero probability has an independent Beta prior.
#' Zeros may be allocated to either the zero atom or a count component;
#' positive observations may only be allocated to count components. A shifted
#' count component can generate zero only when its shift is zero.
#' In the returned draws, component 1 is the zero atom and `eta1` is
#' \eqn{\pi_0}; components 2 through `K + 1` have unconditional weights
#' \eqn{(1-\pi_0)\eta_k}. The atom has `lambda1 = 0` for Poisson, or
#' `mu1 = 0` and fixed `rtilde1 = 0.25` for negative binomial.
#' In either shifted family, the atom also has `kappa1 = 0`.
#' The fixed negative-binomial dispersion is only a representation of the
#' zero atom and is not estimated. Negative-binomial fits require at least
#' one positive observation to define their existing mean prior.
#' Plotting and mode inference include the zero atom; use `range = c(1, max(data))`
#' in [bayes_mode()] to restrict reported modes to positive counts.
#' 
#' 
#' **Mixture of shifted-Poisson distributions**
#' 
#' Shifted-Poisson components take the form
#' \deqn{p(y_i |\lambda_k, \kappa_k) = \frac{1}{(y_i - \kappa_k)!} \,
#' \lambda^{(y_i - \kappa_k)!}_k \,\exp(-\lambda_k)}
#' where \eqn{\kappa_k} is a location or shift parameter with uniform prior, see \insertCite{Cross2024;textual}{BayesMultiMode}.
#' 
#' @references
#' \insertAllCited{}
#' 
#' @importFrom assertthat assert_that
#' @importFrom assertthat is.scalar
#' @importFrom assertthat is.string
#' 
#' @examples
#' fit <- bayes_fit(galaxy, K = 30, dist = "normal", sampler = "telescoping",
#'                  K_init = 3, nb_iter = 100, burnin = 50, print = FALSE)
#' table(fit$state[, "K"], fit$state[, "Kplus"])
#' fit$cap_mass
#'
#' # Example with galaxy data ================================================
#' set.seed(123) 
#' 
#' # retrieve galaxy data
#' y = galaxy
#'
#' # estimation
#' bayesmix = bayes_fit(data = y,
#'                            K = 5, #not many to run the example rapidly
#'                            dist = "normal",
#'                            nb_iter = 500, #not many to run the example rapidly
#'                            burnin = 100)
#'                            
#' # plot estimated mixture
#' # plot(bayesmix, max_size = 200)
#' 
#' # Changing priors ================================================
#' set.seed(123) 
#' 
#' # retrieve galaxy data
#' y = galaxy
#'
#' # estimation
#' K = 5
#' bayesmix = bayes_fit(data = y,
#'                            K = K, #not many to run the example rapidly
#'                            dist = "normal",
#'                            priors = list(a0 = 10,
#'                                          A0 = 10*K),
#'                            nb_iter = 500, #not many to run the example rapidly
#'                            burnin = 100)
#'                            
#' # plot estimated mixture
#' # plot(bayesmix, max_size = 200)
#' 
#' # Example with DNA data =====================================================
#' \donttest{
#' set.seed(123) 
#' 
#' # retrieve DNA data
#' y = d4z4
#'
#' # estimation
#' bayesmix = bayes_fit(data = y,
#'                            K = 5, #not many to run the example rapidly
#'                            dist = "shifted_poisson",
#'                            nb_iter = 500, #not many to run the example rapidly
#'                            burnin = 100)
#'                            
#' # plot estimated mixture
#' # plot(bayesmix, max_size = 200)
#' }
#' 
#' @export
bayes_fit <- function(data,
                      K,
                      dist,
                      priors = list(),
                      nb_iter = 2000,
                      burnin = nb_iter/2,
                      print = TRUE,
                      zero_inflated = FALSE,
                      sampler = c("sfm", "telescoping"),
                      K_init = 1L,
                      weight_model = c("dynamic", "static"),
                      K_prior = c("bnb", "uniform"),
                      concentration = "default") {
  sampler <- match.arg(sampler)
  if (sampler == "telescoping") {
    weight_model <- match.arg(weight_model)
    K_prior <- match.arg(K_prior)
    assert_that(!isTRUE(zero_inflated),
          msg = "Telescoping does not support zero inflation; use sampler = 'sfm' for zero-inflated count mixtures.")
    assert_that(is.numeric(data), is.vector(data), length(data) > 0,
                msg = "Telescoping requires a nonempty numeric vector.")
    assert_that(is.scalar(K_init), is.numeric(K_init), is.finite(K_init), K_init == trunc(K_init),
                K_init >= 1, K_init <= K, msg = "K_init must be an integer between 1 and K.")
    if (is.character(concentration)) {
      # Follow the authors' telescope::priorOnAlpha_spec() and telescope::priorOnE0_spec() defaults.
      # For clustering, the paper suggests "F_6_3" instead of fixing dynamic alpha at 1.
      sampled_prior <- if (weight_model == "dynamic") "F_6_3" else "G_1_20"
      assert_that(is.string(concentration), concentration %in% c("default", sampled_prior),
                  msg = paste0("For ", weight_model, " weights, use concentration = 'default', '",
                               sampled_prior, "', NULL or a positive number."))
      concentration <- if (weight_model == "dynamic" && concentration == "default") 1 else NULL
    }
    if (!is.null(concentration)) assert_that(is.numeric(concentration), is.scalar(concentration),
      is.finite(concentration), concentration > 0, msg = "concentration must be positive and finite.")
    if (any(c("a0", "A0") %in% names(priors))) stop("a0 and A0 are SFM priors; use concentration for telescoping.")
  } else {
    if (!missing(K_init) || !missing(weight_model) || !missing(K_prior) || !missing(concentration))
      stop("K_init, weight_model, K_prior and concentration are telescoping controls.")
    assert_that(is.vector(data) & length(data) > K,
                msg = "data should be a vector of length greater than K")
  }
  
  assert_that(all(is.finite(data)),
              msg = "data should only include numeric finite values")
  assert_that(is.string(dist) & dist %in% c("normal", "skew_normal", "poisson", "shifted_poisson", "neg_binomial", "shifted_neg_binomial"),
              msg = paste0("Unsupported distribution;\n",
                           "dist should be either\n",
                           "'normal', 'skew_normal', 'poisson', 'shifted_poisson', 'neg_binomial' or 'shifted_neg_binomial'"))
  assert_that(is.scalar(nb_iter), round(nb_iter) == nb_iter, nb_iter > 0, msg = "nb_iter should be a positive integer")
  assert_that(is.scalar(burnin), burnin >= as.integer(sampler == "sfm"), burnin < nb_iter, round(burnin) == burnin,
              msg = "burnin should be an integer smaller than nb_iter (and positive for SFM)")
  assert_that(is.scalar(K), is.finite(K), round(K) == K, K > 0, msg = "K should be a positive integer")
  assert_that(is.logical(print), msg = "print should be either TRUE or FALSE")
  assert_that(is.logical(zero_inflated), is.scalar(zero_inflated), !is.na(zero_inflated),
              msg = "zero_inflated should be either TRUE or FALSE")
  
  if (dist %in% c("poisson", "shifted_poisson", "neg_binomial", "shifted_neg_binomial")) {
    assert_that(!any(!data%%1==0),
                msg = "data must include only integer values when using Poisson, shifted Poisson, negative Binomial or shifted negative Binomial mixtures.")
    assert_that(min(data) > -1,
                msg = "data should not include negative values when using Poisson, shifted Poisson, negative Binomial or shifted negative Binomial mixtures.")
    dist_type = "discrete"
  } else {
    assert_that(!zero_inflated, msg = "zero_inflated is supported only for count mixtures.")
    dist_type = "continuous"
  }
  
  priors = check_priors(priors, dist, data, zero_inflated, telescoping = sampler == "telescoping")

  if (sampler == "telescoping") {
    family <- switch(dist, normal = tele_normal_family(data, priors),
             skew_normal = tele_skew_normal_family(data, priors),
             poisson = tele_poisson_family(data, priors, FALSE),
             shifted_poisson = tele_poisson_family(data, priors, TRUE),
             neg_binomial = tele_nb_family(data, priors, burnin, FALSE),
             shifted_neg_binomial = tele_nb_family(data, priors, burnin, TRUE))
    chain <- gibbs_telescoping(data, K, K_init, nb_iter, family, weight_model, K_prior, concentration, print)
    fit <- bayes_mixture(chain$draws, data, burnin, dist, dist_type = dist_type, loglik = chain$loglik)
    fit$sampler <- sampler
    fit$weight_model <- weight_model
    fit$K_prior <- K_prior
    fit$state_all <- chain$state
    fit$state <- chain$state[seq.int(burnin + 1L, nb_iter), , drop = FALSE]
    fit$cap_mass <- mean(fit$state[, "K"] == K)
    return(fit)
  }
  
  mcmc <- gibbs_SFM(y = data,
                    K = K,
                    nb_iter = nb_iter,
                    nb_burnin = burnin,
                    priors = priors,
                    print = print,
                    dist = dist,
                    zero_inflated = zero_inflated)
  
  # extract loglik from mcmc
  ll_id = which(colnames(mcmc) == "loglik")
  loglik = mcmc[,ll_id]
  mcmc = mcmc[,-ll_id]
  
  BayesMixture = bayes_mixture(mcmc = mcmc,
                               data = data,
                               burnin = burnin,
                               dist = dist,
                               dist_type = dist_type,
                               loglik = loglik)
  if (zero_inflated) {
    BayesMixture$zero_inflated = TRUE
    BayesMixture$count_K = K
  }
  
  return(BayesMixture)
}