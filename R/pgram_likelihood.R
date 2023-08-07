pgram_likelihood <- function(sed_rate,
                             layer_boundaries,
                             cyclostrat,
                             target_frequency){
  # this function is an implementation of the probability calculations from
  # Malinverno et al. (2010).
  # INPUTS
  # sed_rate = vector of sedimentation rates
  # layer_boundaries = sedimentation rate change points
  # cyclostrat = cyclostratigraphic record spanning the range of layer_boundaries
  # target_frequency = vector of tuning frequencies to use
  # OUTPUTS
  # LL = log-likelihood of `sed_rate`

  # preallocate storage -------------------------------------------------------
  LL <- vector(length = length(sed_rate))

  for(i in 1:length(sed_rate)) {
    # loop through all layers an calculate scaled periodogram -----------------
    f <- cyclostrat %>%
      filter(position > layer_boundaries[i] &
               position < layer_boundaries[i + 1]) %>%
      astrochron::periodogram(output = 1,
                  verbose = FALSE,
                  genplot = FALSE,
                  background = 1,
                  f0 = TRUE,
                  padfac = 10) %>%
      mutate(probability = (Power / AR1_Fit),
             probability = (probability / sum(probability)),
             time_freq   = Frequency * sed_rate[i]) %>%
      with(approxfun(x   = time_freq,
                     y   = probability))

    # calculate probability of tuning frequencies
    LL[i] <- f(target_frequency) %>%
      log() %>%
      sum()
  }
  LL <- LL %>% sum()
  return(LL)
}

