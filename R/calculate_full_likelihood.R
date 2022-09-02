# this is a wrapper function to calculate model likelihoods
# there are two options `malinverno` and `time_opt`
#

calculate_full_likelihood <- function(sed_rate,         # vector of sedimentation rates
                                      segment_edges,    # vector segment boundaries
                                      cyclostrat_data,  # cyclostrat data
                                      tuning_frequency, # tuning frequencies
                                      method = NA) {

  # use the default method if unspecified
  if(is.na(method)) {method = malinverno}
  # preallocate storage -------------------------------------------------------
  LL <- vector(length = length(sed_rate))

  if(method == 'malinverno') {
    LL <- full_malinverno_likelihood(sed_rate = sed_rate,
                                     segment_edges = segment_edges,
                                     cyclostrat_data = cyclostrat_data,
                                     tuning_frequency = tuning_frequency)
  }

  # if(method == 'time_opt') {
  #   # LL[i] <- time_opt_likelihood(cyclostrat_data = current_cyclostrat,
  #   #                              sed_rate = sed_rate[i],
  #   #                              tuning_frequency = tuning_frequency)$LL_total
  # }


  # calculate joint probability
  LL <- sum(LL)

  return(LL)
}

# INPUTS
# rate = vector of sedimentation rates
# segment_edges = sedimentation rate change points
# cyclostrat = cyclostratigraphic record spanning the range of segment_edges
# tuning frequencies = vector of tuning frequencies to use
full_malinverno_likelihood <- function(sed_rate,
                                       segment_edges,
                                       cyclostrat_data,
                                       tuning_frequency) {

  # prepare the data ----------------------------------------------------------
  # calculate the floating age model
  age_model <- c(0, cumsum(diff(segment_edges) / sed_rate))

  # form an interpolation function
  f <- approxfun(x = segment_edges,
                 y = age_model)
  # apply and overwrite interpolation function
  p <- cyclostrat_data %>%
    # transform position into time
    mutate(position = f(position)) %>%
    # reinterpolate using median spacing
    astrochron::linterp(genplot = FALSE,
                        verbose = FALSE) %>%
    # calculate periodogram
    astrochron::periodogram(output = 1,
                            verbose = FALSE,
                            genplot = FALSE,
                            background = 1,
                            f0 = TRUE) %>%
    # calculate probability distribution
    mutate(probability = Power / AR1_Fit,
           probability = probability / sum(probability),
           time_freq   = Frequency) %>%
    # form a new interpolation function
    with(approxfun(x = time_freq,
                   y = probability))

  # calculate joint probability
  LL <- p(tuning_frequency) %>%
    log() %>%
    sum()

  return(LL)
}
