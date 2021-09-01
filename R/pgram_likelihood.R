# INPUTS
# rate = vector of sedimentation rates
# segment_edges = sedimentation rate change points
# cyclostrat = cyclostratigraphic record spanning the range of segment_edges
# tuning frequencies = vector of tuning frequencies to use
pgram_likelihood <- function(sed_rate,
                             segment_edges,
                             cyclostrat,
                             tuning_frequency){
  # preallocate storage -------------------------------------------------------
  LL <- vector(length = length(sed_rate))

  for(i in 1:length(sed_rate)) {
    # loop through all segments an calculate scaled periodogram ---------------
    f <- cyclostrat %>%
      filter(position > segment_edges[i] &
               position < segment_edges[i + 1]) %>%
      periodogram(output = 1,
                  verbose = FALSE,
                  genplot = FALSE,
                  background = 1) %>%
      mutate(probability = Power / AR1_Fit,
             probability = probability / sum(probability),
             time_freq   = Frequency * sed_rate[i]) %>%
      with(approxfun(x = time_freq,
                     y = probability))
    # calculate probability of tuning frequencies
    LL[i] <- f(tuning_frequency) %>% log() %>% sum()
  }
  LL <- LL %>% sum()
  return(LL)
}
