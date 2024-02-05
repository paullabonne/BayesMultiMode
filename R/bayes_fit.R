#' Bayesian estimation of mixture distributions
#' 
#' Estimation of a univariate mixture with unknown number of components using a sparse finite mixture Markov chain Monte Carlo (SFM MCMC) algorithm.
#' 
#' @param data Vector of observations.
#' @param K Maximum number of mixture components.
#' @param dist String indicating the distribution of the mixture components;
#' currently supports `"normal"`, `"skew_normal"`, `"poisson"` and `"shifted_poisson"`.
#' @param priors List of priors; default is an empty list which implies the following priors:\cr
#' `a0 = 1`,\cr `A0 = 200`,\cr `b0 = median(y)`,\cr `B0 = (max(y) - min(y))^2` (normal),\cr
#' `D_xi = 1`,\cr `D_psi =1`, (skew normal: `B0 = diag(D_xi,D_psi)`), \cr `c0 = 2.5`,\cr
#' `l0 = 1.1` (poisson),\cr `l0 = 5` (shifted poisson),\cr `L0 = 1.1/median(y)`,\cr `L0 = l0 - 1` (shifted poisson),\cr
#' `g0 = 0.5`,\cr `G0 = 100*g0/c0/B0` (normal),\cr 
#' `G0 = g0/(0.5*var(y))` (skew normal).
#' @param nb_iter Number of MCMC iterations; default is `2000`.
#' @param burnin Number of MCMC iterations used as burnin; default is `nb_iter/2`.
#' @param print Showing MCMC progression ? Default is `TRUE`.
#' 
#' @return A list of class \code{bayes_mixture} containing:
#'  \item{data}{Same as argument.}
#'  \item{mcmc}{Matrix of MCMC draws where the rows corresponding to burnin have been discarded;}
#'  \item{mcmc_all}{Matrix of MCMC draws.}
#'  \item{loglik}{Log likelihood at each MCMC draw.}
#'  \item{K}{Number of components.}
#'  \item{dist}{Same as argument.}
#'  \item{pdf_func}{The pdf/pmf of the mixture components.}
#'  \item{dist_type}{Type of the distribution, i.e. continuous or discrete.}
#'  \item{pars_names}{Names of the mixture components' parameters.}
#'  \item{loc}{Name of the location parameter of the mixture components.}
#'  \item{nb_var}{Number of variables/parameters in the mixture distribution.}
#' 
#' @details
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
#' \deqn{\omega_k = \frac{\psi_k}{\sigma_k}, \qquad \omega^2_k = \sigma^2_k + \psi^2_k.}
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
                      print = TRUE) {
  
  assert_that(is.vector(data) & length(data) > K,
              msg = "data should be a vector of length greater than K")
  assert_that(all(is.finite(data)),
              msg = "data should only include numeric finite values")
  assert_that(is.string(dist) & dist %in% c("normal", "skew_normal", "poisson", "shifted_poisson"),
              msg = paste0("Unsupported distribution;\n",
                           "dist should be either\n",
                           "'normal', 'skew_normal', 'poisson' or 'shifted_poisson'"))
  assert_that(is.scalar(nb_iter), round(nb_iter) == nb_iter, nb_iter > 0, msg = "nb_iter should be a positive integer")
  assert_that(is.scalar(burnin), burnin > 0, burnin < nb_iter, round(burnin) == burnin,
              msg = "nb_iter should be a positive integer lower than burnin")
  assert_that(is.scalar(K), round(K) == K, K > 0, msg = "K should be a positive integer")
  assert_that(is.logical(print), msg = "print should be either TRUE or FALSE")
  
  if (dist %in% c("poisson", "shifted_poisson")) {
    assert_that(!any(!data%%1==0),
                msg = "data must include only integer values when using Poisson or shifted Poisson mixtures.")
    assert_that(min(data) > -1,
                msg = "data should not include negative values when using Poisson or shifted Poisson mixtures.")
    dist_type = "discrete"
  } else {
    dist_type = "continuous"
  }
  
  priors = check_priors(priors, dist, data)
  
  mcmc <- gibbs_SFM(y = data,
                    K = K,
                    nb_iter = nb_iter,
                    priors = priors,
                    print = print,
                    dist = dist)
  
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
  
  return(BayesMixture)
}