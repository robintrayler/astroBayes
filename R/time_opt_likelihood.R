# This function is in development for astrochron: An R Package for Astrochronology
# Copyright (C) 2020 Stephen R. Meyers
# this function has been modified to work with the astroBayes package in development
# by Robin B. Trayler, Mark D. Schmitz and Stephen R. Meyers

# sed_rate should be in m/Ma
# target_frequency should be in Ma, 1/Ma, with cycles labeled


#' Haario et al. (2001)
#' @name time_opt_likelihood

#' @import "tidyverse"
#' @import "dplyr"
#' @import "tibble"
#' @importFrom magrittr "%>%"
#'
#' @return a a proposed value from `chain`
#' @md
#' @export
#'
time_opt_likelihood <- function(cyclostrat_data,
                                sed_rate,
                                target_frequency,
                                f_low = 0.035,
                                f_high = 0.065,
                                roll = 10 ^ 3) {

  # convert sed_rate to m/ka
  sed_rate <- sed_rate * 0.001 # convert to m/ka

  # pull out eccentricity frequencies and convert to ka
  target_E <- target_frequency %>%
    filter(orbital_cycle == 'eccentricity') %>%
    pull(period) * 1000

  # pull out precession frequencies and convert to ka
  target_P <- target_frequency %>%
    filter(orbital_cycle == 'precession') %>%
    pull(period) * 1000

  # put it all togehter
  target_total <- c(target_E, target_P)

  # prepare data array --------------------------------------------------------
  cyclostrat_data = data.frame(cyclostrat_data)
  n_pts <- length(cyclostrat_data[, 1])
  dx <- cyclostrat_data[2, 1] - cyclostrat_data[1, 1]

  # standardize data series ---------------------------------------------------
  cyclostrat_data[2]=cyclostrat_data[2] - colMeans(cyclostrat_data[2])
  cyclostrat_data[2]=cyclostrat_data[2] / sapply(cyclostrat_data[2], sd)

  # check minimum and maximum sedimentation rates.
  # sedmin is now in m/ka, dx is in meters.
  NyqFreq <- sed_rate / (2 * dx)
  RayFreq <- sed_rate / (n_pts * dx)

  # CALIBRATE DEPTH SERIES (m) TO TIME (ka) -----------------------------------
  time_series = cyclostrat_data
  # create new time vector
  # it is the index vector for time
  it <- seq(1, n_pts, by = 1)
  time = (dx / sed_rate) * (it - 1)
  time_series[1] = time

  # bandpass precession or short eccentricity band
  bp = taner(time_series,
             padfac=2,
             flow=f_low,
             fhigh=f_high,
             roll=roll,
             demean=T,
             detrend=F,
             addmean=F,
             genplot=F,
             verbose=F)

  # hilbert transform for instantaneous amplitude
  hil = hilbert(bp,
                padfac=2,
                demean=T,
                detrend=F,
                addmean=F,
                genplot=F,
                verbose=F)
  # standardize hil to unit variance
  hil[2] <- hil[2] - colMeans(hil[2])
  hil[2] <- hil[2] / sd(hil[, 2])

  # execute functions ---------------------------------------------------------
  # for precession modulations
  result <- fit_it(sed_rate = sed_rate,
                   time_series = hil,
                   target_in = target_E,
                   n_pts = n_pts,
                   dx = dx)

  power_out <- fit_it(sed_rate,
                      time_series = time_series,
                      target_in = target_total,
                      n_pts = n_pts,
                      dx = dx)

  # for short eccentricity modulations
  answer <- data.frame(
    sed_rate = result$sed_rate * 1000,
    LL_env   = result$logLL,
    LL_pow   = power_out$logLL) %>%
    mutate(LL_total = LL_env + LL_pow)

  return(answer)
}
