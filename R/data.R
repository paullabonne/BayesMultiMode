#' Autosomal macrosatellite repeats d4z4
#'
#' Macrosatellite repeats D4Z4 in the subtelomere of chromosome 4q.\cr
#' Locus (hg18): 4q35.2 \cr
#' Unit (kb): 3.3 \cr
#' Restriction enzyme: EcoRI + HindIII/EcoRI + BlnI/XapI \cr
#' Encoded product : DUX4
#' @format
#' A vector of counts with 410 elements.
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{schaap_genome-wide_2013}{BayesMultiMode}
"d4z4"

#' X chromosomal macrosatellite repeats ct47
#'
#' Repeat units that encode for a cancer testis antigen.\cr
#' Locus (hg18): Xq24 \cr
#' Unit (kb): 4.8 \cr
#' Restriction enzyme: EcoRI \cr
#' Encoded product : cancer testis antigen 47
#'
#' @format
#' A vector of counts with 410 elements.
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{schaap_genome-wide_2013}{BayesMultiMode}
"ct47"

#' Tropical cyclones lifetime maximum intensity
#'
#' Dataset constructed using the International Best Track Archive for Climate Stewardship (IBTrACS).
#' The distribution of tropical cyclones lifetime maximum intensity across the globe is known
#' to be bimodal which has important implications for climate modelling.
#'
#' @format
#' A dataset with three columns showing the identification of the cyclone, its year of occurrence and its lifetime maximum intensity (LMI).
#' LMI is calculated as the maximum wind speed for each cyclone with unit ks.
#'
#' @source
#' https://www.ncei.noaa.gov/products/international-best-track-archive
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{knapp_international_2010}{BayesMultiMode}\cr\cr
#' \insertRef{knapp_international_2018}{BayesMultiMode}
"cyclone"

#' Galaxy series
#'
#' Velocity at which 82 galaxies in the Corona Borealis region are moving away from our galaxy, scaled by 1000.
#'
#' @source
#' https://people.maths.bris.ac.uk/~mapjg/mixdata
#'
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{Richardson_Green_1997_RJMCMC}{BayesMultiMode}
"galaxy"
