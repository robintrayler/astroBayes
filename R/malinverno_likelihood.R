# this function is an implementation of the probability calculations from
# Malinverno et al. (2010).
# INPUTS
# sed_rate = vector of sedimentation rates
# segment_edges = sedimentation rate change points
# cyclostrat = cyclostratigraphic record spanning the range of segment_edges
# tuning frequencies = dataframe of tuning frequencies to use
# OUTPUTS
# LL = log-likelihood of `sed_rate`

malinverno_likelihood <- function(cyclostrat_data,
                                  sed_rate,
                                  tuning_frequency) {
  # calculate probability distribution and make interpolation function
  f <- cyclostrat_data %>%
    periodogram(output = 1,
                verbose = FALSE,
                genplot = FALSE,
                background = 1,
                f0 = TRUE) %>%
    # calculate PDF
    mutate(probability = (Power / AR1_Fit),
           # normalize area to 1
           probability = (probability / sum(probability)),
           # convert cycles/m to Ma/m
           time_freq   = Frequency * sed_rate) %>%
    # make interpolation function
    with(approxfun(x   = time_freq,
                   y   = probability))

  # calculate probability of tuning frequencies
  LL <- f(tuning_frequency$frequency) %>%
    log() %>%
    sum()
  return(LL)
}


