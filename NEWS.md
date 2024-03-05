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