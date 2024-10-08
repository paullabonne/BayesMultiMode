% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayes_trace.R
\name{bayes_trace}
\alias{bayes_trace}
\title{Trace plots}
\usage{
bayes_trace(BayesMix, mcmc_vars = NULL, with_burnin = FALSE, ...)
}
\arguments{
\item{BayesMix}{An object of class \code{bayes_mixture}.}

\item{mcmc_vars}{Variables to plot; default is all the variable in the MCMC output.}

\item{with_burnin}{Plot all draws ?}

\item{...}{Additional arguments passed to function \code{\link[bayesplot:MCMC-traces]{bayesplot::mcmc_trace()}}.}
}
\value{
A trace plot.
}
\description{
This is wrapper around the \code{\link[bayesplot:MCMC-traces]{bayesplot::mcmc_trace()}} function from package \code{bayesplot}.
}
\examples{
# Example with galaxy data ================================================
set.seed(123)

# retrieve galaxy data
y <- galaxy

# estimation
bayesmix <- bayes_fit(
  data = y,
  K = 5, # not many to run the example rapidly
  dist = "normal",
  nb_iter = 500, # not many to run the example rapidly
  burnin = 100
)

# trace plot
# bayes_trace(bayesmix)

}
