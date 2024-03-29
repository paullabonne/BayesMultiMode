% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{cyclone}
\alias{cyclone}
\title{Tropical cyclones lifetime maximum intensity}
\format{
A dataset with three columns showing the identification of the cyclone, its year of occurrence and its lifetime maximum intensity (LMI).
LMI is calculated as the maximum wind speed for each cyclone with unit ks.
}
\source{
https://www.ncei.noaa.gov/products/international-best-track-archive
}
\usage{
cyclone
}
\description{
Dataset constructed using the International Best Track Archive for Climate Stewardship (IBTrACS).
The distribution of tropical cyclones lifetime maximum intensity across the globe is known
to be bimodal which has important implications for climate modelling.
}
\references{
\insertRef{knapp_international_2010}{BayesMultiMode}\cr\cr
\insertRef{knapp_international_2018}{BayesMultiMode}
}
\keyword{datasets}
