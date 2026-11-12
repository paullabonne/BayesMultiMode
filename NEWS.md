## Development version
* Poisson and negative-binomial mixtures now support `zero_inflated = TRUE`,
	adding a fixed zero atom with an independent Beta prior controlled by `zero_a`
	and `zero_b`. `K` specifies the count components; returned component 1 is the
	zero atom and `eta1` is its structural-zero probability. Count weights retain
	the sparse-mixture prior. Existing fits are unchanged when the option is off.
* Mode-count probabilities retain their matrix dimensions when every posterior
	draw has the same number of modes, including zero-inflated fits.
* Shifted Poisson mixtures now jointly update the shift and rate of each non-empty
	component at fixed `kappa + lambda`, using the likelihood and existing Gamma rate
	prior. The subsequent conditional Gamma update can change the total component
	mean. Empty-component prior draws and output formats are unchanged; seeded draw
	sequences change.
* Normal, skew-normal, Poisson, and shifted Poisson samplers now use the same
	vectorized allocation routine as the negative binomial samplers. Poisson PMFs
	are evaluated at distinct counts and mapped back to the observations, including
	shift likelihoods. Priors, parameter updates, and output formats are unchanged;
	vectorized allocation changes seeded draw sequences. Continuous densities are
	not deduplicated.
* Negative binomial and shifted negative binomial samplers now evaluate PMFs at
	distinct counts and map values back to the original observations. Allocation
	uses vectorized inverse-CDF sampling instead of per-observation multinomial
	calls. Per-observation allocation state, priors, and output formats are unchanged.
	PMF deduplication preserves seeded calculations; vectorized allocation changes
	the random-number sequence while preserving the categorical probabilities.
* Negative binomial and shifted negative binomial mixtures now tune `rtilde`
	proposal scales every 50 burn-in iterations, targeting 30-40% acceptance with
	diminishing adjustments. Empty components do not contribute to proposal counts,
	and adaptation stops after burn-in. Component sorting preserves proposal scales
	and batch histories. Priors and the joint shift-mean update are unchanged.
* Shifted negative binomial mixtures now jointly update the shift and mean of each
	non-empty component. The update samples all valid integer shifts conditional on
	`kappa + mu` and adjusts `mu` to preserve that sum, using the existing uniform
	priors and fixed full-data bounds. The subsequent mean and shape updates are
	unchanged; the component mean is not fixed across iterations.

## Version 0.7.1
* Range argument not optional any more when plotting mixtures
* More details in summary methods
* Added a df to the mode output showing the mixture density in each draw
* Added message for users to be aware of potential label-switching to summary and plot_trace functions.

## Version 0.7.0
* Major changes around the structure of the package
* bayes_estimation() renamed to bayes_fit()
* new_BayesMixture() renamed to bayes_mixture()
* Added a class "mixture" representing estimated mixtures which is used as input to mode estimation functions; see mixture()
* Added mix_mode() which calls the mode finding algorithms; discrete_MF, fixed_point and MEM have been removed
* Added a class "mix_mode"
* Removed the plotting option inside the mode estimation functions
* Added print, plot and summary methods for all classes.

## Version 0.6.0
* Added examples to new_BayesMixture() and bayes_mode()
* Fixed typo on the tolerance arguments
* Added details to mode and mixture estimation functions
* Improved the documentation

## Version 0.5.1
* Improved robustness of the gibbs sampler when initial classification is bad

## Version 0.5.0 
* Restructured the code to make use of S3 generic functions plot() and summary()
* The package now handles continuous data
* Added a mixtures of normals and skew_normals
* Added galaxy and cyclone data
* Added support for external MCMC output
* Added tests through testthats